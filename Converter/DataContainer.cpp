
#include <algorithm>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <stdexcept>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "DataContainer.h"
#include "GridUtility.h"

using namespace amrex;

namespace {
std::string tec_upper(std::string value) {
  for (char& c : value)
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return value;
}

std::string tec_trim(const std::string& value) {
  const auto first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos)
    return {};
  return value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
}

// Comments and assignments in quoted titles must never become header tokens.
std::string tec_strip_comment(const std::string& line) {
  bool quoted = false;
  for (size_t i = 0; i < line.size(); ++i) {
    if (line[i] == '"') {
      size_t backslashes = 0;
      for (size_t j = i; j > 0 && line[j - 1] == '\\'; --j)
        ++backslashes;
      if (backslashes % 2 == 0)
        quoted = !quoted;
    }
    if (line[i] == '#' && !quoted)
      return line.substr(0, i);
  }
  return line;
}

struct TecToken {
  std::string value;
  bool quoted = false;
};

std::vector<TecToken> tec_tokens(const std::string& text) {
  std::vector<TecToken> tokens;
  for (size_t i = 0; i < text.size();) {
    const unsigned char c = text[i];
    if (std::isspace(c) || c == ',') {
      ++i;
      continue;
    }
    TecToken token;
    if (c == '"') {
      token.quoted = true;
      ++i;
      bool closed = false;
      while (i < text.size()) {
        if (text[i] == '"') {
          ++i;
          closed = true;
          break;
        }
        if (text[i] == '\\' && i + 1 < text.size() &&
            (text[i + 1] == '"' || text[i + 1] == '\\'))
          ++i;
        token.value += text[i++];
      }
      if (!closed)
        throw std::runtime_error("unterminated quoted header text");
    } else if (c == '=') {
      token.value = text[i++];
    } else if (c == '(' || c == '[') {
      // Preserve variable-location expressions as one assignment value.
      const size_t start = i;
      int depth = 0;
      do {
        if (text[i] == '(' || text[i] == '[')
          ++depth;
        else if (text[i] == ')' || text[i] == ']')
          --depth;
        ++i;
      } while (i < text.size() && depth > 0);
      if (depth != 0)
        throw std::runtime_error("unbalanced header expression");
      token.value = text.substr(start, i - start);
    } else {
      const size_t start = i;
      while (i < text.size() &&
             !std::isspace(static_cast<unsigned char>(text[i])) &&
             text[i] != ',' && text[i] != '=' && text[i] != '"')
        ++i;
      token.value = text.substr(start, i - start);
    }
    tokens.push_back(std::move(token));
  }
  return tokens;
}

size_t tec_size(const std::string& text, const std::string& description) {
  size_t value = 0;
  if (text.empty())
    throw std::runtime_error("missing " + description);
  for (const unsigned char c : text) {
    if (!std::isdigit(c))
      throw std::runtime_error("invalid integer " + description + ": " + text);
    const size_t digit = c - '0';
    if (value > (std::numeric_limits<size_t>::max() - digit) / 10)
      throw std::runtime_error("integer overflow in " + description);
    value = 10 * value + digit;
  }
  return value;
}

size_t tec_product(size_t count, size_t width, size_t maxSize,
                   const std::string& description) {
  if (width != 0 && count > maxSize / width)
    throw std::runtime_error("size overflow in " + description);
  return count * width;
}

bool tec_next_number(std::istream& stream, std::string& token) {
  while (stream >> token) {
    if (token[0] == '#') {
      stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }
    return true;
  }
  return false;
}
} // namespace

int TECDataContainer::read() {
  pointValues.clear();
  connectivity.clear();
  varNames.clear();
  nCell = nBrick = 0;
  nVar = nDim = 0;
  elementType = ZoneType::Type::UNSET;
  std::fill(std::begin(coordinateColumns), std::end(coordinateColumns), -1);

  try {
    std::ifstream input(filename, std::ios::binary);
    if (!input)
      throw std::runtime_error("cannot open input file");

    std::string variables, zone;
    enum class Record { NONE, VARIABLES, ZONE, METADATA };
    Record record = Record::NONE;
    bool hasVariables = false, hasZone = false, hasData = false;
    std::streampos dataStart;
    std::string line;
    while (input >> std::ws) {
      const auto lineStart = input.tellg();
      const int first = input.peek();
      if ((first >= '0' && first <= '9') || first == '+' || first == '-' ||
          first == '.') {
        dataStart = lineStart;
        hasData = true;
        break;
      }
      if (!std::getline(input, line))
        break;
      line = tec_trim(tec_strip_comment(line));
      if (line.empty())
        continue;
      const auto wordEnd = line.find_first_of(" \t=,");
      const std::string word = tec_upper(line.substr(0, wordEnd));
      if (word == "VARIABLES") {
        if (hasVariables || hasZone)
          throw std::runtime_error("duplicate or misplaced VARIABLES record");
        hasVariables = true;
        record = Record::VARIABLES;
        variables += line.substr(9) + "\n";
      } else if (word == "ZONE") {
        if (hasZone)
          throw std::runtime_error("multiple zones are unsupported");
        hasZone = true;
        record = Record::ZONE;
        zone += line.substr(4) + "\n";
      } else if (word == "TITLE" || word == "AUXDATA" ||
                 word == "DATASETAUXDATA") {
        // Validate quotes even in metadata that is deliberately not retained.
        tec_tokens(line);
        record = Record::METADATA;
      } else if (record == Record::VARIABLES) {
        variables += line + "\n";
      } else if (record == Record::ZONE) {
        zone += line + "\n";
      } else {
        throw std::runtime_error(
            "unsupported header record or non-ASCII payload");
      }
    }
    if (!hasVariables || !hasZone)
      throw std::runtime_error("required VARIABLES or ZONE record is missing");

    const auto variableTokens = tec_tokens(variables);
    size_t firstVariable = 0;
    if (!variableTokens.empty() && variableTokens[0].value == "=")
      firstVariable = 1;
    for (size_t i = firstVariable; i < variableTokens.size(); ++i) {
      const auto& token = variableTokens[i];
      if (token.value.empty() || (!token.quoted && token.value == "="))
        throw std::runtime_error("invalid variable label");
      varNames.push_back(token.value);
    }
    if (varNames.empty() || varNames.size() > static_cast<size_t>(INT_MAX))
      throw std::runtime_error("invalid or overflowing variable count");
    nVar = static_cast<int>(varNames.size());

    std::map<std::string, std::string> properties;
    const auto zoneTokens = tec_tokens(zone);
    for (size_t i = 0; i < zoneTokens.size();) {
      if (i + 2 >= zoneTokens.size() || zoneTokens[i].quoted ||
          zoneTokens[i + 1].value != "=")
        throw std::runtime_error("malformed ZONE assignment");
      const std::string key = tec_upper(zoneTokens[i].value);
      if (!properties.emplace(key, zoneTokens[i + 2].value).second)
        throw std::runtime_error("duplicate ZONE property " + key);
      i += 3;
    }
    for (const auto& property : properties) {
      const auto& key = property.first;
      if (key == "VARLOCATION") {
        const auto location = tec_upper(property.second);
        if (location.find("CELLCENTERED") != std::string::npos ||
            location.find("NODAL") == std::string::npos)
          throw std::runtime_error("only nodal VARLOCATION is supported");
      } else if (key != "T" && key != "N" && key != "NODES" && key != "E" &&
                 key != "ELEMENTS" && key != "F" && key != "DATAPACKING" &&
                 key != "ET" && key != "ZONETYPE" && key != "SOLUTIONTIME" &&
                 key != "STRANDID") {
        throw std::runtime_error(
            "unsupported ZONE property " + key +
            " (ordered grids and sharing are unsupported)");
      }
    }
    auto count = [&](const std::string& legacy, const std::string& modern) {
      if (properties.count(legacy) && properties.count(modern))
        throw std::runtime_error("duplicate count aliases " + legacy + "/" +
                                 modern);
      const auto found =
          properties.find(properties.count(legacy) ? legacy : modern);
      if (found == properties.end())
        throw std::runtime_error("missing finite-element count " + legacy);
      const size_t value = tec_size(found->second, legacy);
      if (value == 0)
        throw std::runtime_error("empty meshes are unsupported: " + legacy +
                                 "=0");
      return value;
    };
    const size_t points = count("N", "NODES");
    const size_t elements = count("E", "ELEMENTS");
    bool hasPacking = false;
    for (const auto& key : { "F", "DATAPACKING" }) {
      if (!properties.count(key))
        continue;
      hasPacking = true;
      const auto packing = tec_upper(properties.at(key));
      if ((std::string(key) == "F" && packing != "FEPOINT") ||
          (std::string(key) == "DATAPACKING" && packing != "POINT"))
        throw std::runtime_error(
            "only finite-element POINT/FEPOINT packing is supported");
    }
    if (!hasPacking)
      throw std::runtime_error("missing POINT/FEPOINT packing declaration");
    for (const auto& key : { "ET", "ZONETYPE" }) {
      if (!properties.count(key))
        continue;
      const auto name = tec_upper(properties.at(key));
      ZoneType::Type type = ZoneType::Type::UNSET;
      if (name == "BRICK" || name == "FEBRICK")
        type = ZoneType::Type::BRICK;
      else if (name == "QUADRILATERAL" || name == "FEQUADRILATERAL")
        type = ZoneType::Type::QUAD;
      else
        throw std::runtime_error("unsupported finite-element topology " + name);
      if (elementType != ZoneType::Type::UNSET && elementType != type)
        throw std::runtime_error("conflicting element topology declarations");
      elementType = type;
    }
    if (elementType == ZoneType::Type::UNSET)
      throw std::runtime_error("missing BRICK or QUADRILATERAL element type");

    for (int i = 0; i < nVar; ++i) {
      const auto label = tec_upper(tec_trim(varNames[i]));
      if (label.empty())
        continue;
      const auto coordinate = std::string("XYZ").find(label[0]);
      if (coordinate == std::string::npos)
        continue;
      if (label.size() > 1 &&
          !std::isspace(static_cast<unsigned char>(label[1])) &&
          label[1] != '[' && label[1] != '(' && label[1] != '{')
        continue;
      if (coordinateColumns[coordinate] != -1)
        throw std::runtime_error("duplicate coordinate variable " + label);
      coordinateColumns[coordinate] = i;
    }
    if (coordinateColumns[0] < 0 || coordinateColumns[1] < 0 ||
        (elementType == ZoneType::Type::BRICK && coordinateColumns[2] < 0))
      throw std::runtime_error(
          "missing required X/Y coordinates (and Z for BRICK)");
    nDim = coordinateColumns[2] >= 0 ? 3 : 2;
    const size_t vertices = elementType == ZoneType::Type::BRICK ? 8 : 4;
    const size_t values =
        tec_product(points, nVar, pointValues.max_size(), "point values");
    const size_t indices = tec_product(elements, vertices,
                                       connectivity.max_size(), "connectivity");
    tec_product(points, 3, pointValues.max_size(), "coordinates");
    if (values > std::numeric_limits<size_t>::max() - indices)
      throw std::runtime_error("size overflow in numeric token count");
    if (!hasData)
      throw std::runtime_error("numeric point data is missing");
    input.clear();
    input.seekg(0, std::ios::end);
    const auto end = input.tellg();
    if (end < dataStart ||
        static_cast<uintmax_t>(end - dataStart) / 2 + 1 < values + indices)
      throw std::runtime_error("truncated data: declared counts exceed the "
                               "available numeric payload");
    input.seekg(dataStart);

    pointValues.resize(values);
    connectivity.resize(indices);
    std::string token;
    for (size_t i = 0; i < values; ++i) {
      if (!tec_next_number(input, token))
        throw std::runtime_error("truncated point data at point " +
                                 std::to_string(i / nVar + 1));
      for (char& c : token)
        if (c == 'D' || c == 'd')
          c = 'E';
      char* last = nullptr;
      const float value = std::strtof(token.c_str(), &last);
      if (last == token.c_str() || last != token.c_str() + token.size() ||
          !std::isfinite(value))
        throw std::runtime_error("invalid ASCII value at point " +
                                 std::to_string(i / nVar + 1));
      pointValues[i] = value;
    }
    for (size_t i = 0; i < indices; ++i) {
      if (!tec_next_number(input, token))
        throw std::runtime_error("truncated connectivity at element " +
                                 std::to_string(i / vertices + 1));
      const size_t index = tec_size(token, "connectivity");
      if (index == 0 || index > points)
        throw std::runtime_error("connectivity out of range at element " +
                                 std::to_string(i / vertices + 1));
      connectivity[i] = index;
    }
    if (tec_next_number(input, token))
      throw std::runtime_error(
          "unexpected trailing payload or multiple zones are unsupported");
    if (input.bad())
      throw std::runtime_error("I/O error reading numeric payload");
    nCell = points;
    nBrick = elements;
    return iSuccess;
  } catch (const std::exception& error) {
    Print() << "Error reading Tecplot ASCII " << filename << ": "
            << error.what() << '\n';
    pointValues.clear();
    connectivity.clear();
    return iFail;
  }
}

void TECDataContainer::get_loc(Vector<float>& vars) {
  vars.resize(3 * nCell);
  for (size_t i = 0; i < nCell; ++i)
    for (int coordinate = 0; coordinate < 3; ++coordinate)
      vars[3 * i + coordinate] =
          coordinateColumns[coordinate] >= 0
              ? pointValues[i * nVar + coordinateColumns[coordinate]]
              : 0.0f;
}

void AMReXDataContainer::read_header(const std::string& headerName, int& nVar,
                                     int& nDim, Real& time, int& finest_level,
                                     RealBox& domain, Box& cellBox,
                                     Vector<std::string>& varNames) {
  BL_PROFILE("AMReXDataContainer::read_header");

  std::ifstream HeaderFile(headerName, std::ifstream::in);
  if (!HeaderFile.is_open()) {
    Abort("Error: cannot open header file: " + headerName);
  }

  HeaderFile.precision(17);

  std::string versionName;

  HeaderFile >> versionName;
  HeaderFile >> nVar;

  varNames.clear();
  varNames.reserve(nVar);
  for (int ivar = 0; ivar < nVar; ++ivar) {
    std::string var;
    HeaderFile >> var;
    varNames.push_back(var);
  }

  HeaderFile >> nDim;
  HeaderFile >> time;
  HeaderFile >> finest_level;

  for (int i = 0; i < nDim; ++i) {
    Real lo;
    HeaderFile >> lo;
    domain.setLo(i, lo);
  }

  for (int i = 0; i < nDim; ++i) {
    Real hi;
    HeaderFile >> hi;
    domain.setHi(i, hi);
  }

  for (int i = 0; i < finest_level; ++i) {
    int refRatio;
    HeaderFile >> refRatio;
  }

  if (finest_level == 0) {
    HeaderFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  // Read the base grid and ignore refined level grids.
  HeaderFile >> cellBox;
  HeaderFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

void AMReXDataContainer::read_header() {
  const std::string headerName = filename + "/Header";

  int finestLev;
  AMReXDataContainer::read_header(headerName, nVar, nDim, time, finestLev,
                                  domain, cellBox, varNames);

  SetFinestLevel(finestLev);
}

int AMReXDataContainer::read() {
  BL_PROFILE("AMReXDataContainer::read");

  Print() << "Reading in " << filename << std::endl;

  read_header();

  nCell = 0;
  nBrick = 0;
  isCellNumbered = false;

  Grid grid(Geom(0), get_amr_info(), nGst);

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    VisMF::Read(mf[iLev],
                filename + "/Level_" + std::to_string(iLev) + "/Cell");

    grid.SetBoxArray(iLev, mf[iLev].boxArray());
    grid.SetDistributionMap(iLev, mf[iLev].DistributionMap());
  }
  grid.SetFinestLevel(n_lev() - 1);

  regrid(grid.boxArray(0), &grid);

  return iSuccess;
}

size_t AMReXDataContainer::loop_cell(bool doStore, Vector<float>& vars,
                                     bool doStoreLoc) {
  BL_PROFILE("AMReXDataContainer::loop_cell");

  if (doStore) {
    vars.clear();
    if (nCell > 0) {
      vars.reserve(nCell *
                   (doStoreLoc ? 3 : (n_lev() > 0 ? mf[0].nComp() : 1)));
    }
  }

  const bool needNumbering = !isCellNumbered;
  size_t iCount = 0;

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    const int ncomp = mf[iLev].nComp();
    const auto& geom = Geom(iLev);
    const auto dx = geom.CellSize();
    const auto probLo = geom.ProbLo();

    for (MFIter mfi(iCell[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto& status = cell_status(iLev)[mfi].array();
      const Array4<Real>& cell = iCell[iLev][mfi].array();
      const Array4<Real>& data = mf[iLev][mfi].array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k) {
#if AMREX_SPACEDIM > 2
        const float z = static_cast<float>(probLo[iz_] + (k + 0.5) * dx[iz_]);
#else
        constexpr float z = 0.0f;
#endif
        for (int j = lo.y; j <= hi.y; ++j) {
          const float y = static_cast<float>(probLo[iy_] + (j + 0.5) * dx[iy_]);
          for (int i = lo.x; i <= hi.x; ++i) {
            if (!bit::is_refined(status(i, j, k))) {
              iCount++;
              if (needNumbering) {
                cell(i, j, k) = static_cast<Real>(iCount);
              }
              if (doStore) {
                if (doStoreLoc) {
                  const float x =
                      static_cast<float>(probLo[ix_] + (i + 0.5) * dx[ix_]);
                  vars.push_back(x);
                  vars.push_back(y);
                  vars.push_back(z);
                } else {
                  for (int iVar = 0; iVar < ncomp; iVar++)
                    vars.push_back(static_cast<float>(data(i, j, k, iVar)));
                }
              }
            }
          }
        }
      }
    }

    if (needNumbering) {
      iCell[iLev].FillBoundary();
    }
  }

  if (needNumbering) {
    for (int iLev = n_lev() - 2; iLev >= 0; --iLev) {
      fill_fine_lev_bny_from_coarse(
          iCell[iLev], iCell[iLev + 1], 0, iCell[iLev].nComp(), ref_ratio[iLev],
          Geom(iLev), Geom(iLev + 1), cell_status(iLev + 1), pc_interp);
    }
    isCellNumbered = true;
  }

  nCell = iCount;
  return iCount;
}

size_t AMReXDataContainer::loop_zone(bool doStore, Vector<size_t>& zones) {
  BL_PROFILE("AMReXDataContainer::loop_zone");

  size_t iBrick = 0;

  if (doStore) {
    zones.clear();
    if (nBrick > 0) {
      zones.reserve(nBrick * 8);
    }
  }

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    for (MFIter mfi(iCell[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& cell = iCell[iLev][mfi].array();
      const auto& status = cell_status(iLev)[mfi].array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      // Loop over valid cells + lower end ghost cells
      for (int k = lo.z - 1; k <= hi.z; ++k) {
        for (int j = lo.y - 1; j <= hi.y; ++j) {
          for (int i = lo.x - 1; i <= hi.x; ++i) {
            const bool inValidBox = (i >= lo.x && j >= lo.y && k >= lo.z);
            if ((iLev > 0 && bit::is_lev_boundary(status(i, j, k))) ||
                (inValidBox && !bit::is_refined(status(i, j, k)))) {
              bool isBrick = true;

              for (int kk = k; kk <= k + 1 && isBrick; ++kk) {
                for (int jj = j; jj <= j + 1 && isBrick; ++jj) {
                  for (int ii = i; ii <= i + 1; ++ii) {
                    if (cell(ii, jj, kk) == 0) {
                      isBrick = false;
                      break;
                    }
                  }
                }
              }

              if (isBrick) {
                iBrick++;
                if (doStore) {
                  zones.push_back(static_cast<size_t>(cell(i, j, k)));
                  zones.push_back(static_cast<size_t>(cell(i + 1, j, k)));
                  zones.push_back(static_cast<size_t>(cell(i + 1, j + 1, k)));
                  zones.push_back(static_cast<size_t>(cell(i, j + 1, k)));
                  zones.push_back(static_cast<size_t>(cell(i, j, k + 1)));
                  zones.push_back(static_cast<size_t>(cell(i + 1, j, k + 1)));
                  zones.push_back(
                      static_cast<size_t>(cell(i + 1, j + 1, k + 1)));
                  zones.push_back(static_cast<size_t>(cell(i, j + 1, k + 1)));
                }
              }
            }
          }
        }
      }
    }
  }

  nBrick = iBrick;
  return iBrick;
}

void AMReXDataContainer::smooth(int nSmooth) {
  BL_PROFILE("AMReXDataContainer::smooth");

  constexpr Real coef = 0.5;
  constexpr Real weightSelf = 1.0 - coef;
  constexpr Real weightNei = coef / 2.0;

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    const int ng = 1;
    MultiFab mfOld(mf[iLev].boxArray(), mf[iLev].DistributionMap(),
                   mf[iLev].nComp(), ng);

    for (int iSmooth = 0; iSmooth < nSmooth; iSmooth++) {
      auto smooth_dir = [&](int iDir) {
        MultiFab::Copy(mfOld, mf[iLev], 0, 0, mf[iLev].nComp(), 0);
        mfOld.FillBoundary(Geom(iLev).periodicity());

        int dIdx[3] = { 0, 0, 0 };
        dIdx[iDir] = 1;
        const int di = dIdx[ix_];
        const int dj = dIdx[iy_];
        const int dk = dIdx[iz_];

        const int ncomp = mf[iLev].nComp();

        for (MFIter mfi(mf[iLev]); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();

          const auto& status = cell_status(iLev)[mfi].array();
          Array4<Real> const& arr = mf[iLev][mfi].array();
          Array4<Real> const& tmp = mfOld[mfi].array();

          const auto lo = lbound(bx);
          const auto hi = ubound(bx);

          for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
              for (int i = lo.x; i <= hi.x; ++i) {
                if (bit::is_lev_edge(status(i, j, k))) {
                  continue;
                }

                auto is_bnd_or_refined = [&](int ci, int cj, int ck) {
                  const auto st = status(ci, cj, ck);
                  return bit::is_lev_boundary(st) || bit::is_refined(st);
                };

                if (is_bnd_or_refined(i - di, j - dj, k - dk) ||
                    is_bnd_or_refined(i + di, j + dj, k + dk)) {
                  continue;
                }

                for (int iVar = 0; iVar < ncomp; ++iVar) {
                  const Real neiSum = tmp(i - di, j - dj, k - dk, iVar) +
                                      tmp(i + di, j + dj, k + dk, iVar);
                  arr(i, j, k, iVar) =
                      weightSelf * arr(i, j, k, iVar) + weightNei * neiSum;
                }
              }
            }
          }
        }

        mf[iLev].FillBoundary(Geom(iLev).periodicity());
      };

      smooth_dir(ix_);
      smooth_dir(iy_);
      smooth_dir(iz_);
    }
  }
}

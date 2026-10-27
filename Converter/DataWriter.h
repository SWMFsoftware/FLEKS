#ifndef _DATAWRITER_H_
#define _DATAWRITER_H_

#include <algorithm>
#include <cassert>
#include <cctype>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <type_traits>
#include <utility>

#include "DataContainer.h"
#include "Parallel.h"
#include "ZLibCompressor.h"

// Whether std::to_chars(first, last, T, chars_format, precision) exists.
// Library version checks are not enough: e.g. libstdc++ 11 hides the
// floating-point overloads from compilers such as Intel icpc.
template <typename T, typename = void>
struct has_float_to_chars : std::false_type {};
template <typename T>
struct has_float_to_chars<
    T, std::void_t<decltype(std::to_chars(
           std::declval<char*>(), std::declval<char*>(), std::declval<T>(),
           std::chars_format::general, 0))>> : std::true_type {};

template <typename T> char* format_float_impl(char* p, T v, int digits) {
  if constexpr (has_float_to_chars<T>::value)
    return std::to_chars(p, p + 32, v, std::chars_format::general, digits).ptr;
  else
    return p + std::snprintf(p, 32, "%.*g", digits, static_cast<double>(v));
}

// Write v at p as printf("%.*g", digits, v) does, and return the end of the
// text. At most 32 characters are written.
inline char* format_float(char* p, float v, int digits) {
  return format_float_impl(p, v, digits);
}

// Write n bytes in pieces. A single write of 2 GiB or more fails on macOS
// (write(2) rejects sizes above INT_MAX), which large meshes reach.
inline void write_bytes(std::ostream& os, const void* data, uint64_t n) {
  constexpr uint64_t maxChunk = uint64_t(1) << 30;
  const char* p = static_cast<const char*>(data);
  while (n > 0 && os) {
    const uint64_t chunk = std::min(n, maxChunk);
    os.write(p, static_cast<std::streamsize>(chunk));
    p += chunk;
    n -= chunk;
  }
}

// Prints the wall-clock time of consecutive phases of a task, e.g.
// "[timing] TEC: write cells: 1.234 s". phase() ends the running phase and
// starts the next one; the destructor ends the last one.
class PhaseTimer {
public:
  explicit PhaseTimer(std::string taskIn) : task(std::move(taskIn)) {}

  ~PhaseTimer() { end_phase(); }

  void phase(const std::string& name) {
    end_phase();
    current = name;
    start = std::chrono::steady_clock::now();
  }

  // End the running phase without starting a new one.
  void stop() { end_phase(); }

  PhaseTimer(const PhaseTimer&) = delete;
  PhaseTimer& operator=(const PhaseTimer&) = delete;

private:
  void end_phase() {
    if (current.empty())
      return;
    const std::chrono::duration<double> dt =
        std::chrono::steady_clock::now() - start;
    const auto oldFlags = std::cout.flags();
    const auto oldPrecision = std::cout.precision(3);
    std::cout << "[timing] " << task << ": " << current << ": " << std::fixed
              << dt.count() << " s" << std::endl;
    std::cout.precision(oldPrecision);
    std::cout.flags(oldFlags);
    current.clear();
  }

  std::string task;
  std::string current;
  std::chrono::steady_clock::time_point start;
};

class DataWriter {
public:
  DataWriter(DataContainer* dcIn, const std::string& filenameIn) {
    dc = dcIn;
    filename = filenameIn;
    sourceFilename = filenameIn;
  };

  virtual ~DataWriter() {};

  virtual int write() = 0;

  std::string type_string() { return fileTypeString.at(fType); }

  virtual void print() {
    std::cout << "========DataWriter========\n";
    std::cout << "Write data to: " << filename << "\n";
    std::cout << "File type: " << type_string() << "\n";
    std::cout << "================================" << std::endl;
  }

protected:
  bool can_write() const { return can_write_to(filename); }

  bool can_write_to(const std::filesystem::path& output) const {
    const auto reject = [&](const std::string& reason) {
      std::cerr << "Error: cannot write output " << output << ": " << reason
                << "\n";
      return false;
    };
    std::error_code error;
    // weakly_canonical cannot resolve a symlink whose target does not exist.
    // Opening such an output could create a file inside the source directory.
    const auto status = std::filesystem::symlink_status(output, error);
    if (error && error != std::errc::no_such_file_or_directory)
      return reject("cannot inspect path: " + error.message());
    error.clear();
    const bool exists = std::filesystem::exists(output, error);
    if (error)
      return reject("cannot inspect target: " + error.message());
    if (std::filesystem::is_symlink(status) && !exists)
      return reject("dangling symlink");

    const auto src = std::filesystem::weakly_canonical(sourceFilename, error);
    if (error)
      return reject("cannot resolve input path: " + error.message());
    const auto dst = std::filesystem::weakly_canonical(output, error);
    if (error)
      return reject("cannot resolve output path: " + error.message());
    // Compare path components, so directory inputs cannot contain the output
    // that -D would subsequently remove. Canonical paths also resolve symlinks.
    const auto mismatch =
        std::mismatch(src.begin(), src.end(), dst.begin(), dst.end());
    if (mismatch.first == src.end())
      return reject("aliases or lies inside the input");
    if (exists && std::filesystem::equivalent(src, dst, error))
      return reject("aliases the input file");
    if (error)
      return reject("cannot check file identity: " + error.message());
    return true;
  }

  DataContainer* dc;
  std::string filename;
  std::string sourceFilename;
  FileType fType;
  std::ofstream outFile;

  ZoneType zoneType;
};

class TECWriter : public DataWriter {
public:
  TECWriter(DataContainer* dcIn, const std::string& filenameIn)
      : DataWriter(dcIn, filenameIn) {
    fType = FileType::TEC;
    filename = filenameIn + ".dat";
  };

  ~TECWriter() {};

  int write() override {

    if (!can_write())
      return iFail;

    zoneType.set_type(dc->zone_type());
    if (zoneType.n_vertex() <= 0) {
      std::cerr << "Error: unsupported zone topology.\n";
      return iFail;
    }

    PhaseTimer timer("TEC");

    timer.phase("count cells and zones");
    size_t nCell = dc->count_cell();
    size_t nBrick = dc->count_zone();

    timer.phase("get cell data");
    amrex::Vector<float> vars;
    dc->get_cell(vars);

    timer.phase("get zones");
    amrex::Vector<size_t> zones;
    dc->get_zones(zones);

    const int nVar = dc->n_var();
    const int nVertex = zoneType.n_vertex();
    if (vars.size() != nCell * nVar || zones.size() != nBrick * nVertex) {
      std::cerr << "Error: inconsistent cell or zone data.\n";
      return iFail;
    }

    timer.phase("write header");
    outFile.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
    if (!outFile) {
      std::cerr << "Error: cannot open output: " << filename << "\n";
      return iFail;
    }
    outFile.precision(std::numeric_limits<float>::max_digits10);

    //-----------Write header---------------
    outFile << "TITLE = " << std::quoted(filename) << "\n";

    outFile << "VARIABLES = ";

    auto varNames = dc->var_names();
    for (int i = 0; i < varNames.size(); ++i) {
      outFile << std::quoted(varNames[i]);
      if (i != varNames.size() - 1) {
        outFile << ',' << " ";
      }
    }
    outFile << "\n";

    outFile << "ZONE "
            << " N=" << nCell << ", E=" << nBrick
            << ", F=FEPOINT, ET=" << zoneType.tec_string() << "\n";
    //-----------------------------------------

    // Format the numbers into text. format_float() gives the same text as the
    // stream used for the header (%.9g), but is much faster. Blocks of rows
    // are formatted by several threads and written in their original order.
    constexpr int floatDigits = std::numeric_limits<float>::max_digits10;

    // formatRow(i, p) writes row i at p and returns the new end pointer. A row
    // of n numbers takes at most rowMax(n) characters.
    auto rowMax = [](int n) { return 33 * static_cast<size_t>(n) + 1; };
    auto write_rows = [&](size_t nRow, size_t rowBytes, auto&& formatRow) {
      if (nRow == 0)
        return;
      // About 4 MiB of text per block bounds the buffer memory per thread.
      constexpr size_t blockBytes = size_t(4) << 20;
      const size_t rowsPerBlock =
          std::min(nRow, std::max<size_t>(1, blockBytes / rowBytes));
      const size_t nBlock = (nRow + rowsPerBlock - 1) / rowsPerBlock;
      const size_t nThread = std::min(converter_num_threads(), nBlock);
      std::vector<std::vector<char> > bufs(nThread);
      std::vector<size_t> sizes(nThread);
      for (auto& b : bufs)
        b.resize(rowsPerBlock * rowBytes);

      for (size_t first = 0; first < nBlock; first += nThread) {
        const size_t nActive = std::min(nThread, nBlock - first);
        run_parallel(nActive, [&](size_t t) {
          const size_t iBegin = (first + t) * rowsPerBlock;
          const size_t iEnd = std::min(nRow, iBegin + rowsPerBlock);
          char* const begin = bufs[t].data();
          char* p = begin;
          for (size_t i = iBegin; i < iEnd; ++i)
            p = formatRow(i, p);
          sizes[t] = p - begin;
        });
        for (size_t t = 0; t < nActive; ++t)
          outFile.write(bufs[t].data(), sizes[t]);
      }
    };

    // Write cell data
    timer.phase("write cells");
    write_rows(nCell, rowMax(nVar), [&](size_t i, char* p) {
      const float* row = vars.data() + i * nVar;
      for (int j = 0; j < nVar; ++j) {
        p = format_float(p, row[j], floatDigits);
        *p++ = ' ';
      }
      *p++ = '\n';
      return p;
    });

    // Write zone data
    timer.phase("write zones");
    write_rows(nBrick, rowMax(nVertex), [&](size_t i, char* p) {
      const size_t* row = zones.data() + i * nVertex;
      for (int j = 0; j < nVertex; ++j) {
        p = std::to_chars(p, p + 32, row[j]).ptr;
        *p++ = ' ';
      }
      *p++ = '\n';
      return p;
    });

    timer.phase("close file");
    if (outFile.is_open()) {
      outFile.close();
    }

    if (!outFile) {
      std::cerr << "Error: writing output failed: " << filename << "\n";
      return iFail;
    }
    return iSuccess;
  }
};

class VTKWriter : public DataWriter {
public:
  VTKWriter(DataContainer* dcIn, const std::string& filenameIn)
      : DataWriter(dcIn, filenameIn) {
    fType = FileType::VTK;
    filename = filenameIn + ".vtk";

    saveBinary = true;
  };

  ~VTKWriter() {};

  int write() override {
    if (!can_write())
      return iFail;
    zoneType.set_type(dc->zone_type());

    PhaseTimer timer("VTK");

    timer.phase("count cells and zones");
    size_t nCell = dc->count_cell();
    size_t nBrick = dc->count_zone();
    const int nVertex = zoneType.n_vertex();
    const size_t maxInt = std::numeric_limits<int>::max();
    if (nVertex <= 0 || nCell > maxInt / 3 ||
        nBrick > maxInt / static_cast<size_t>(nVertex + 1)) {
      std::cerr << "Error: mesh exceeds legacy VTK integer limits.\n";
      return iFail;
    }

    timer.phase("get cell data");
    amrex::Vector<float> vars;
    dc->get_cell(vars);

    timer.phase("get locations");
    amrex::Vector<float> xyz;
    // If nDim == 2, the coordinates of the third dimension are 0.
    dc->get_loc(xyz);
    assert(xyz.size() == 3 * nCell);

    //=== Brick data ===
    timer.phase("get zones");
    amrex::Vector<size_t> zones;
    amrex::Vector<int> bricksInt;
    dc->get_zones(zones);
    bricksInt.resize(zones.size());
    {
      bool isValid = true;
      for (size_t n = 0; n < zones.size(); ++n) {
        const size_t index = zones[n];
        isValid &= (index != 0 && index <= nCell);
        bricksInt[n] = static_cast<int>(index - 1);
      }
      if (!isValid) {
        std::cerr << "Error: invalid mesh connectivity.\n";
        return iFail;
      }
    }
    // The size_t copy is no longer needed.
    zones = amrex::Vector<size_t>();
    int* brickData = bricksInt.data();
    //===================

    amrex::Vector<int> brickType(nBrick, zoneType.vtk_index());

    // All variables are scalars, so vardim is 1.
    amrex::Vector<int> vardim(dc->n_var(), 1);

    // Cell-based: 0  (This is connectivity cell, NOT simulation cell)
    // Point-based: 1
    amrex::Vector<int> centering(dc->n_var(), 1);

    amrex::Vector<amrex::Vector<char> > varnameStorage(dc->n_var());
    amrex::Vector<char*> varnames(dc->n_var());
    const auto varNames = dc->var_names();
    std::set<std::string> usedNames;
    for (int i = 0; i < dc->n_var(); ++i) {
      std::string base = varNames[i];
      for (auto& ch : base)
        if (std::isspace(static_cast<unsigned char>(ch)))
          ch = '_';
      if (base.empty())
        base = "variable";
      std::string name = base;
      for (size_t suffix = 2; !usedNames.insert(name).second; ++suffix)
        name = base + "_" + std::to_string(suffix);
      if (name != varNames[i])
        std::cout << "VTK variable: \"" << varNames[i] << "\" -> \"" << name
                  << "\"\n";
      varnameStorage[i].assign(name.begin(), name.end());
      varnameStorage[i].push_back('\0');
      varnames[i] = varnameStorage[i].data();
    }

    timer.phase("transpose variables");
    const int nVarOut = dc->n_var();
    amrex::Vector<amrex::Vector<float> > valueStorage(
        nVarOut, amrex::Vector<float>(nCell));
    amrex::Vector<float*> v(nVarOut);
    for (int i = 0; i < nVarOut; ++i) {
      v[i] = valueStorage[i].data();
    }
    // Walk the interleaved data once instead of once per variable.
    for (size_t j = 0; j < nCell; ++j) {
      const float* row = vars.data() + j * nVarOut;
      for (int i = 0; i < nVarOut; ++i) {
        v[i][j] = row[i];
      }
    }
    // The interleaved copy is no longer needed.
    vars = amrex::Vector<float>();

    timer.phase("write file");
    if (!write_unstructured_mesh(
            filename.c_str(), saveBinary, static_cast<int>(nCell), xyz.data(),
            static_cast<int>(nBrick), brickType.data(), brickData, dc->n_var(),
            vardim.data(), centering.data(), varnames.data(), v.data())) {
      std::cerr << "Error: writing output failed: " << filename << "\n";
      return iFail;
    }

    return iSuccess;
  }

private:
  bool saveBinary;
};

class VTMWriter : public DataWriter {
public:
  VTMWriter(DataContainer* dcIn, const std::string& filenameIn,
            bool useCompressionIn = false)
      : DataWriter(dcIn, filenameIn), useCompression(useCompressionIn) {
    auto srcPath = std::filesystem::path(sourceFilename).lexically_normal();
    // A trailing separator (including one left by /.) has an empty filename.
    while (srcPath.has_relative_path() && srcPath.filename().empty())
      srcPath = srcPath.parent_path();
    const std::string name = srcPath.filename().string();
    const std::string dirName = name + "_vtm";
    const std::filesystem::path parentDir = srcPath.parent_path();
    outputDir = parentDir.empty() ? std::filesystem::path(dirName)
                                  : parentDir / dirName;
    filename = (outputDir / (name + ".vtm")).string();
    pieceFileName = name + "_0.vtu";
    fType = FileType::VTM;
  }
  ~VTMWriter() override = default;

  void print() override {
    std::cout << "========DataWriter========\n";
    std::cout << "Write data to: " << filename << "\n";
    std::cout << "File type: " << type_string()
              << (useCompression ? " (compressed)" : "") << "\n";
    std::cout << "================================" << std::endl;
  }

  int write() override {
    zoneType.set_type(dc->zone_type());
    if (zoneType.n_vertex() <= 0) {
      std::cerr << "Error: unsupported zone topology.\n";
      return iFail;
    }

    PhaseTimer timer("VTM");

    timer.phase("count cells and zones");
    const size_t nCell = dc->count_cell();
    const size_t nBrick = dc->count_zone();
    const int nVertex = zoneType.n_vertex();

    const std::filesystem::path pieceFilePath = outputDir / pieceFileName;
    if (!can_write_to(outputDir) || !can_write_to(filename) ||
        !can_write_to(pieceFilePath))
      return iFail;

    std::error_code ec;
    std::filesystem::create_directories(outputDir, ec);
    if (ec) {
      std::cerr << "Error: cannot create directory " << outputDir << ": "
                << ec.message() << "\n";
      return iFail;
    }

    // Points
    timer.phase("get locations");
    amrex::Vector<float> xyz;
    dc->get_loc(xyz);
    assert(xyz.size() == 3 * nCell);

    // Connectivity
    timer.phase("get zones");
    amrex::Vector<size_t> zones;
    dc->get_zones(zones);
    amrex::Vector<int64_t> connectivity(zones.size());
    {
      bool isValid = true;
      for (size_t n = 0; n < zones.size(); ++n) {
        const size_t index = zones[n];
        isValid &= (index != 0 && index <= nCell);
        connectivity[n] = static_cast<int64_t>(index - 1);
      }
      if (!isValid) {
        std::cerr << "Error: invalid mesh connectivity.\n";
        return iFail;
      }
    }
    // The size_t copy is no longer needed.
    zones = amrex::Vector<size_t>();

    // Offsets and cell types
    amrex::Vector<int64_t> offsets(nBrick);
    for (size_t i = 0; i < nBrick; ++i) {
      offsets[i] = static_cast<int64_t>((i + 1) * nVertex);
    }
    const uint8_t vtkType = static_cast<uint8_t>(zoneType.vtk_index());
    amrex::Vector<uint8_t> cellTypes(nBrick, vtkType);

    // Variables
    const int nVar = dc->n_var();
    const auto varNames = dc->var_names();
    std::vector<std::string> cleanNames(nVar);
    std::set<std::string> usedNames;
    for (int i = 0; i < nVar; ++i) {
      std::string base = varNames[i];
      for (auto& ch : base)
        if (std::isspace(static_cast<unsigned char>(ch)))
          ch = '_';
      if (base.empty())
        base = "variable";
      std::string name = base;
      for (size_t suffix = 2; !usedNames.insert(name).second; ++suffix)
        name = base + "_" + std::to_string(suffix);
      if (name != varNames[i])
        std::cout << "VTK variable: \"" << varNames[i] << "\" -> \"" << name
                  << "\"\n";
      cleanNames[i] = name;
    }

    timer.phase("get cell data");
    amrex::Vector<float> vars;
    dc->get_cell(vars);

    // Determine machine endianness
    const uint16_t endianTest = 1;
    const char* endianStr =
        (*reinterpret_cast<const uint8_t*>(&endianTest) == 1) ? "LittleEndian"
                                                              : "BigEndian";

    // Compute appended raw binary block offsets
    std::vector<uint64_t> varOffsets(nVar);
    uint64_t pointsOffset = 0;
    uint64_t connOffset = 0;
    uint64_t offsetsOffset = 0;
    uint64_t typesOffset = 0;

    const uint64_t pointsBytes =
        static_cast<uint64_t>(3 * nCell) * sizeof(float);
    const uint64_t connBytes =
        static_cast<uint64_t>(connectivity.size()) * sizeof(int64_t);
    const uint64_t offsetsBytes =
        static_cast<uint64_t>(offsets.size()) * sizeof(int64_t);
    const uint64_t typesBytes =
        static_cast<uint64_t>(cellTypes.size()) * sizeof(uint8_t);

    std::vector<std::vector<uint8_t> > varBuffers;
    std::vector<uint8_t> pointsBuf;
    std::vector<uint8_t> connBuf;
    std::vector<uint8_t> offsetsBuf;
    std::vector<uint8_t> typesBuf;

    timer.phase(useCompression ? "compress" : "compute offsets");
    if (useCompression) {
      ZLibCompressor comp;
      varBuffers.resize(nVar);
      uint64_t currentOffset = 0;
      amrex::Vector<float> varBuf(nCell);
      for (int i = 0; i < nVar; ++i) {
        for (size_t j = 0; j < nCell; ++j) {
          varBuf[j] = vars[j * nVar + i];
        }
        varOffsets[i] = currentOffset;
        if (!comp.compress_vtk_block(
                varBuf.data(), static_cast<uint64_t>(nCell) * sizeof(float),
                varBuffers[i])) {
          std::cerr << "Error: failed to compress variable " << cleanNames[i]
                    << "\n";
          return iFail;
        }
        currentOffset += varBuffers[i].size();
      }

      pointsOffset = currentOffset;
      if (!comp.compress_vtk_block(xyz.data(), pointsBytes, pointsBuf)) {
        std::cerr << "Error: failed to compress points.\n";
        return iFail;
      }
      currentOffset += pointsBuf.size();

      connOffset = currentOffset;
      if (!comp.compress_vtk_block(connectivity.data(), connBytes, connBuf)) {
        std::cerr << "Error: failed to compress connectivity.\n";
        return iFail;
      }
      currentOffset += connBuf.size();

      offsetsOffset = currentOffset;
      if (!comp.compress_vtk_block(offsets.data(), offsetsBytes, offsetsBuf)) {
        std::cerr << "Error: failed to compress offsets.\n";
        return iFail;
      }
      currentOffset += offsetsBuf.size();

      typesOffset = currentOffset;
      if (!comp.compress_vtk_block(cellTypes.data(), typesBytes, typesBuf)) {
        std::cerr << "Error: failed to compress cell types.\n";
        return iFail;
      }
      currentOffset += typesBuf.size();
    } else {
      const uint64_t headerBytes = sizeof(uint64_t);
      uint64_t currentOffset = 0;
      for (int i = 0; i < nVar; ++i) {
        varOffsets[i] = currentOffset;
        const uint64_t byteSize = static_cast<uint64_t>(nCell) * sizeof(float);
        currentOffset += headerBytes + byteSize;
      }
      pointsOffset = currentOffset;
      currentOffset += headerBytes + pointsBytes;
      connOffset = currentOffset;
      currentOffset += headerBytes + connBytes;
      offsetsOffset = currentOffset;
      currentOffset += headerBytes + offsetsBytes;
      typesOffset = currentOffset;
      currentOffset += headerBytes + typesBytes;
    }

    // 1. Write the .vtu Piece file
    timer.phase("write piece file");
    std::ofstream vtu(pieceFilePath.string(),
                      std::ios::out | std::ios::binary | std::ios::trunc);
    if (!vtu) {
      std::cerr << "Error: cannot open piece file for writing: "
                << pieceFilePath << "\n";
      return iFail;
    }

    vtu << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\""
        << endianStr << "\" header_type=\"UInt64\""
        << (useCompression ? " compressor=\"vtkZLibDataCompressor\"" : "")
        << ">\n"
        << "  <UnstructuredGrid>\n"
        << "    <Piece NumberOfPoints=\"" << nCell << "\" NumberOfCells=\""
        << nBrick << "\">\n"
        << "      <PointData>\n";

    for (int i = 0; i < nVar; ++i) {
      vtu << "        <DataArray type=\"Float32\" Name=\""
          << xml_attribute(cleanNames[i]) << "\" format=\"appended\" offset=\""
          << varOffsets[i] << "\"/>\n";
    }

    vtu << "      </PointData>\n"
        << "      <CellData>\n"
        << "      </CellData>\n"
        << "      <Points>\n"
        << "        <DataArray type=\"Float32\" Name=\"Points\" "
           "NumberOfComponents=\"3\" format=\"appended\" offset=\""
        << pointsOffset << "\"/>\n"
        << "      </Points>\n"
        << "      <Cells>\n"
        << "        <DataArray type=\"Int64\" Name=\"connectivity\" "
           "format=\"appended\" offset=\""
        << connOffset << "\"/>\n"
        << "        <DataArray type=\"Int64\" Name=\"offsets\" "
           "format=\"appended\" offset=\""
        << offsetsOffset << "\"/>\n"
        << "        <DataArray type=\"UInt8\" Name=\"types\" "
           "format=\"appended\" offset=\""
        << typesOffset << "\"/>\n"
        << "      </Cells>\n"
        << "    </Piece>\n"
        << "  </UnstructuredGrid>\n"
        << "  <AppendedData encoding=\"raw\">\n"
        << "    _";

    if (useCompression) {
      for (int i = 0; i < nVar; ++i) {
        write_bytes(vtu, varBuffers[i].data(), varBuffers[i].size());
      }
      write_bytes(vtu, pointsBuf.data(), pointsBuf.size());
      write_bytes(vtu, connBuf.data(), connBuf.size());
      write_bytes(vtu, offsetsBuf.data(), offsetsBuf.size());
      write_bytes(vtu, typesBuf.data(), typesBuf.size());
    } else {
      auto write_block = [&](const void* data, uint64_t sizeInBytes) {
        write_bytes(vtu, &sizeInBytes, sizeof(uint64_t));
        write_bytes(vtu, data, sizeInBytes);
      };

      // Write PointData arrays
      amrex::Vector<float> varBuf(nCell);
      for (int i = 0; i < nVar; ++i) {
        for (size_t j = 0; j < nCell; ++j) {
          varBuf[j] = vars[j * nVar + i];
        }
        write_block(varBuf.data(),
                    static_cast<uint64_t>(nCell) * sizeof(float));
      }

      // Write Points, connectivity, offsets, types
      write_block(xyz.data(), pointsBytes);
      write_block(connectivity.data(), connBytes);
      write_block(offsets.data(), offsetsBytes);
      write_block(cellTypes.data(), typesBytes);
    }

    vtu << "\n  </AppendedData>\n</VTKFile>\n";
    vtu.close();

    if (!vtu) {
      std::cerr << "Error: writing piece file failed: " << pieceFilePath
                << "\n";
      return iFail;
    }

    // 2. Write the .vtm Manifest file
    timer.phase("write manifest");
    std::ofstream vtm(filename.c_str(),
                      std::ofstream::out | std::ofstream::trunc);
    if (!vtm) {
      std::cerr << "Error: cannot open VTM file for writing: " << filename
                << "\n";
      return iFail;
    }

    vtm << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" "
           "byte_order=\""
        << endianStr << "\" header_type=\"UInt64\">\n"
        << "  <vtkMultiBlockDataSet>\n"
        << "    <DataSet index=\"0\" name=\"Block0\" file=\""
        << xml_attribute(pieceFileName) << "\"/>\n"
        << "  </vtkMultiBlockDataSet>\n"
        << "</VTKFile>\n";
    vtm.close();

    if (!vtm) {
      std::cerr << "Error: writing VTM file failed: " << filename << "\n";
      return iFail;
    }

    return iSuccess;
  }

private:
  static std::string xml_attribute(const std::string& value) {
    std::string escaped;
    for (const char ch : value) {
      switch (ch) {
        case '&':
          escaped += "&amp;";
          break;
        case '<':
          escaped += "&lt;";
          break;
        case '>':
          escaped += "&gt;";
          break;
        case '"':
          escaped += "&quot;";
          break;
        case '\'':
          escaped += "&apos;";
          break;
        case '\t':
          escaped += "&#9;";
          break;
        case '\n':
          escaped += "&#10;";
          break;
        case '\r':
          escaped += "&#13;";
          break;
        default:
          escaped += ch;
          break;
      }
    }
    return escaped;
  }

  std::string pieceFileName;
  std::filesystem::path outputDir;
  bool useCompression = false;
};

#endif

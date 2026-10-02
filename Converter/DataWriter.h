#ifndef _DATAWRITER_H_
#define _DATAWRITER_H_

#include <cassert>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>

#include "DataContainer.h"
#include "ZLibCompressor.h"

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
  bool can_write() const {
    std::error_code error;
    // weakly_canonical resolves symlinks and normalises without requiring
    // that the output file already exists (unlike std::filesystem::equivalent).
    const auto src = std::filesystem::weakly_canonical(sourceFilename, error);
    const auto dst = std::filesystem::weakly_canonical(filename, error);
    if (src == dst) {
      std::cerr << "Error: output aliases the input file: " << filename << "\n";
      return false;
    }
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

    size_t nCell = dc->count_cell();
    size_t nBrick = dc->count_zone();

    amrex::Vector<float> vars;
    vars.resize(nCell * dc->n_var());
    dc->get_cell(vars);

    amrex::Vector<size_t> zones;

    zones.resize(nBrick * zoneType.n_vertex());

    dc->get_zones(zones);

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

    // Write cell data
    for (size_t i = 0; i < nCell; ++i) {
      for (int j = 0; j < dc->n_var(); ++j) {
        outFile << vars[i * dc->n_var() + j] << " ";
      }
      outFile << "\n";
    }

    // Write zone data
    for (size_t i = 0; i < nBrick; ++i) {
      for (int j = 0; j < zoneType.n_vertex(); ++j) {
        outFile << zones[i * zoneType.n_vertex() + j] << " ";
      }
      outFile << "\n";
    }

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

    size_t nCell = dc->count_cell();
    size_t nBrick = dc->count_zone();
    const int nVertex = zoneType.n_vertex();
    const size_t maxInt = std::numeric_limits<int>::max();
    if (nVertex <= 0 || nCell > maxInt / 3 ||
        nBrick > maxInt / static_cast<size_t>(nVertex + 1)) {
      std::cerr << "Error: mesh exceeds legacy VTK integer limits.\n";
      return iFail;
    }

    amrex::Vector<float> vars;
    vars.resize(nCell * dc->n_var());
    dc->get_cell(vars);

    amrex::Vector<float> xyz;
    // If nDim == 2, set the coordinates of the third dimension to 0.
    xyz.resize(3 * nCell);
    dc->get_loc(xyz);
    assert(xyz.size() == 3 * nCell);

    //=== Brick data ===
    amrex::Vector<size_t> zones;
    amrex::Vector<int> bricksInt;
    zones.resize(nBrick * zoneType.n_vertex());
    dc->get_zones(zones);
    bricksInt.reserve(zones.size());
    for (const auto index : zones) {
      if (index == 0 || index > nCell) {
        std::cerr << "Error: invalid mesh connectivity.\n";
        return iFail;
      }
      bricksInt.push_back(static_cast<int>(index - 1));
    }
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

    amrex::Vector<amrex::Vector<float> > valueStorage(
        dc->n_var(), amrex::Vector<float>(nCell));
    amrex::Vector<float*> v(dc->n_var());
    for (int i = 0; i < dc->n_var(); ++i) {
      v[i] = valueStorage[i].data();
      for (size_t j = 0; j < nCell; ++j) {
        v[i][j] = vars[j * dc->n_var() + i];
      }
    }

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
    const std::filesystem::path srcPath(sourceFilename);
    const std::string name = srcPath.filename().string();
    const std::string dirName = name + "_vtm";
    const std::filesystem::path parentDir = srcPath.parent_path();
    outputDir = parentDir.empty() ? std::filesystem::path(dirName)
                                  : parentDir / dirName;
    filename = (outputDir / (name + ".vtm")).string();
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

    const size_t nCell = dc->count_cell();
    const size_t nBrick = dc->count_zone();
    const int nVertex = zoneType.n_vertex();

    const std::filesystem::path srcPath(sourceFilename);
    const std::string name = srcPath.filename().string();
    const std::string pieceFileName = name + "_0.vtu";
    const std::filesystem::path pieceFilePath = outputDir / pieceFileName;

    // Check alias guard for output directory, vtm file, and piece file
    std::error_code ec;
    const auto src = std::filesystem::weakly_canonical(sourceFilename, ec);
    const auto outDirCanon = std::filesystem::weakly_canonical(outputDir, ec);
    const auto vtmCanon = std::filesystem::weakly_canonical(filename, ec);
    const auto pieceCanon =
        std::filesystem::weakly_canonical(pieceFilePath, ec);
    if (src == outDirCanon || src == vtmCanon || src == pieceCanon) {
      std::cerr << "Error: output aliases the input file: " << sourceFilename
                << "\n";
      return iFail;
    }

    std::filesystem::create_directories(outputDir, ec);
    if (ec) {
      std::cerr << "Error: cannot create directory " << outputDir << ": "
                << ec.message() << "\n";
      return iFail;
    }

    // Points
    amrex::Vector<float> xyz;
    xyz.resize(3 * nCell);
    dc->get_loc(xyz);
    assert(xyz.size() == 3 * nCell);

    // Connectivity
    amrex::Vector<size_t> zones;
    zones.resize(nBrick * nVertex);
    dc->get_zones(zones);
    amrex::Vector<int64_t> connectivity;
    connectivity.reserve(zones.size());
    for (const auto index : zones) {
      if (index == 0 || index > nCell) {
        std::cerr << "Error: invalid mesh connectivity.\n";
        return iFail;
      }
      connectivity.push_back(static_cast<int64_t>(index - 1));
    }

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

    amrex::Vector<float> vars;
    vars.resize(nCell * nVar);
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
      vtu << "        <DataArray type=\"Float32\" Name=\"" << cleanNames[i]
          << "\" format=\"appended\" offset=\"" << varOffsets[i] << "\"/>\n";
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
        vtu.write(reinterpret_cast<const char*>(varBuffers[i].data()),
                  varBuffers[i].size());
      }
      vtu.write(reinterpret_cast<const char*>(pointsBuf.data()),
                pointsBuf.size());
      vtu.write(reinterpret_cast<const char*>(connBuf.data()), connBuf.size());
      vtu.write(reinterpret_cast<const char*>(offsetsBuf.data()),
                offsetsBuf.size());
      vtu.write(reinterpret_cast<const char*>(typesBuf.data()),
                typesBuf.size());
    } else {
      auto write_block = [&](const void* data, uint64_t sizeInBytes) {
        vtu.write(reinterpret_cast<const char*>(&sizeInBytes),
                  sizeof(uint64_t));
        vtu.write(reinterpret_cast<const char*>(data), sizeInBytes);
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
        << "    <DataSet index=\"0\" name=\"Block0\" file=\"" << pieceFileName
        << "\"/>\n"
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
  std::filesystem::path outputDir;
  bool useCompression = false;
};

#endif

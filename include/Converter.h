#ifndef _CONVERTER_H_
#define _CONVERTER_H_

#include "DataWriter.h"

class Converter {
public:
  Converter(const std::string& in) {
    sourceFile = in;
    sourceType = find_source_type();

    switch (sourceType) {
      case FileType::AMREX: {
        amrex::Geometry gm;
        amrex::AmrInfo info;
        amrex::RealBox domain;
        amrex::Box cellBox;

        { // Read simulation domain, and cell number from header.
          std::string headerName = sourceFile + "/Header";

          int nVar, nDim, finest_lev;
          amrex::Real time;
          amrex::Vector<std::string> varNames;

          AMReXDataContainer::read_header(headerName, nVar, nDim, time,
                                          finest_lev, domain, cellBox,
                                          varNames);

          info.max_level = finest_lev;
          info.blocking_factor.clear();
          for (int iLev = 0; iLev <= info.max_level; iLev++) {
            info.blocking_factor.push_back(amrex::IntVect(1));
          }
        }

        gm.define(cellBox, &domain);

        dc = std::make_unique<AMReXDataContainer>(sourceFile, gm, info);

        dc->print();
        break;
      }
      case FileType::IDL: {
        // dc = std::make_unique<IDLDataContainer>(sourceFile);
        break;
      }
      case FileType::TECPLOT: {
        // dc = std::make_unique<TecplotDataContainer>(sourceFile);
        break;
      }
    }

    destType = FileType::TECPLOT;

    switch (destType) {
      case FileType::IDL: {
        break;
      }
      case FileType::TECPLOT: {
        writer = std::make_unique<TECWriter>(dc.get(), sourceFile);
        break;
      }
    }
  }

  FileType find_source_type() {
    if (sourceFile.size() > 6 &&
        sourceFile.substr(sourceFile.size() - 6, 6) == "_amrex") {
      // *_amrex
      return FileType::AMREX;
    } else if (sourceFile.size() > 4 &&
               sourceFile.substr(sourceFile.size() - 4, 4) == ".dat") {
      // *.dat
      return FileType::TECPLOT;
    } else if (sourceFile.size() > 4 &&
               sourceFile.substr(sourceFile.size() - 4, 4) == ".out") {
      // *.out
      return FileType::IDL;
    } else {
      amrex::Abort("Unknown sourceFile type");
    }
  }
  void read() { dc->read(); }
  void write() { writer->write(); }

private:
  std::string sourceFile;
  FileType sourceType;
  std::unique_ptr<DataContainer> dc;

  FileType destType;
  std::unique_ptr<DataWriter> writer;
};

#endif

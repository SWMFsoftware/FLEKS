#ifndef _DATAWRITER_H_
#define _DATAWRITER_H_

#include "DataContainer.h"

class DataWriter {
public:
  DataWriter(DataContainer* dcIn, const std::string& filenameIn) {
    dc = dcIn;
    filename = filenameIn;
  };

  virtual ~DataWriter();

  virtual void write();

  void print() {}

protected:
  DataContainer* dc;
  std::string filename;
  FileType fType;
};

#endif
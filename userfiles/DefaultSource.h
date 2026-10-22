#ifndef _USERSOURCE_H_
#define _USERSOURCE_H_

#include "SourceInterface.h"

class UserSource : public SourceInterface {
public:
  UserSource(Grid& gridIn, const FluidInterface& other, int id, std::string tag,
             FluidType typeIn)
      : SourceInterface(gridIn, other, id, tag, typeIn) {
    info = "Empty user source class";
  }
};

#endif

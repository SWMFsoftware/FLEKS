#include "DomainGrid.h"

using namespace amrex;

void DomainGrid::init() {
  // nCell and boxRange should have been set before calling this function.

  for (int i = 0; i < nDim; i++) {
    centerBoxLo[i] = 0;
    centerBoxHi[i] = nCell[i] - 1;
  }

  centerBox.setSmall(centerBoxLo);
  centerBox.setBig(centerBoxHi);

  geom.define(centerBox, &boxRange, coord, periodicity);

  centerBA.define(centerBox);
  centerBA.maxSize(maxBlockSize);

  dm.define(centerBA);

  nodeBA = convert(centerBA, amrex::IntVect{ AMREX_D_DECL(1, 1, 1) });

  costMF.define(centerBA, dm, 1, 0);
  costMF.setVal(0);

  amrex::Print() << "DomainGrid:: Domain range = " << boxRange << std::endl;
  amrex::Print() << "DomainGrid:: Total block #  = " << nodeBA.size()
                 << std::endl;
  // amrex::Print() << "DomainGrid:: centerBA = " << centerBA << std::endl;
}

void DomainGrid::resize_pic(){
    BoxArray baPic; 
}
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

  Print() << "DomainGrid:: Domain range = " << boxRange << std::endl;
  Print() << "DomaxinGrid:: Total block #  = " << nodeBA.size() << std::endl;
  // amrex::Print() << "DomainGrid:: centerBA = " << centerBA << std::endl;
}

BoxArray DomainGrid::resize_pic_ba(int iCycle) {
  BoxArray baPic;

  if (iCycle == 0) {
    baPic = centerBA;
  } else {

    IntVect quarterCell;
    for (int i = 0; i < nDim; ++i) {
      quarterCell[i] = (centerBoxHi[i] - centerBoxLo[i] + 1) / 4;
    }

    Box bxPic;
    if (iCycle % 2 == 0) {
      bxPic.setSmall(quarterCell / 2);
      bxPic.setBig(centerBoxHi - quarterCell / 2);
    } else {
      bxPic.setSmall(quarterCell);
      bxPic.setBig(centerBoxHi - quarterCell);
    }

    baPic.define(bxPic);
    baPic.maxSize(maxBlockSize);
  }

  Print() << "baPic = " << baPic << std::endl;

  return baPic;
}
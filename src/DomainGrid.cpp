#include "DomainGrid.h"
#include "GridUtility.h"

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
  std::string nameFunc = "Pic::resize_pic_ba";
  Print() << nameFunc << " is runing..." << std::endl;

  Vector<IntVect> tagVec;

  // Loop through the cells in the whole FLEKS domain on each processor.
  for (int i = centerBoxLo[ix_]; i <= centerBoxHi[ix_]; i++)
    for (int j = centerBoxLo[iy_]; j <= centerBoxHi[iy_]; j++)
      for (int k = centerBoxLo[iz_]; k <= centerBoxHi[iz_]; k++) {
        if (gridInfo.get_status(i, j, k) == GridInfo::iPicOn_) {
          tagVec.push_back({ i, j, k });
          Print()<<"tagVec: i = "<<i<<" j = "<<j<<" k = "<<k<<std::endl;
        }
      }

  BoxArray baNew; 
  Print()<<"baPicOld = "<<baPicOld<<std::endl;
  add_cells_to_BoxArray(baNew, tagVec); 
  baNew.maxSize(maxBlockSize); 
  Print()<<"baNew = "<<baNew<<std::endl;
  baPicOld = baNew; 

  return baNew; 

  // BoxArray baPic;

  // if (true || iCycle == 0) {
  //   baPic = centerBA;
  // } else {

  //   IntVect quarterCell;
  //   for (int i = 0; i < nDim; ++i) {
  //     quarterCell[i] = (centerBoxHi[i] - centerBoxLo[i] + 1) / 4;
  //   }

  //   Box bxPic;
  //   if (iCycle % 2 == 0) {
  //     bxPic.setSmall(quarterCell / 2);
  //     bxPic.setBig(centerBoxHi - quarterCell / 2);
  //   } else {
  //     bxPic.setSmall(quarterCell);
  //     bxPic.setBig(centerBoxHi - quarterCell);
  //   }

  //   baPic.define(bxPic);
  //   baPic.maxSize(maxBlockSize);
  // }

  // Print() << "baPic = " << baPic << std::endl;

  // return baPic;
}
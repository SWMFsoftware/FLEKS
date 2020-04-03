#include "DomainGrid.h"
#include "GridUtility.h"
#include "Timer.h"

using namespace amrex;

//=====================================================
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

  Print() << "Domain range = " << boxRange << std::endl;
}

BoxArray DomainGrid::resize_pic_ba(int iCycle) {  
  std::string nameFunc = "Pic::resize_pic_ba";
  Timer funcTimer(nameFunc); 
  Print() << nameFunc << " is runing..." << std::endl;

  // BoxList blInActive;

  // const int di = gridInfo.get_patch_size(ix_);
  // const int dj = gridInfo.get_patch_size(iy_);
  // const int dk = gridInfo.get_patch_size(iz_);

  // IntVect minBoxLo = centerBoxHi, minBoxHi = centerBoxLo;

  // // Loop through the cells in the whole FLEKS domain on each processor.
  // for (int i = centerBoxLo[ix_]; i <= centerBoxHi[ix_]; i += di)
  //   for (int j = centerBoxLo[iy_]; j <= centerBoxHi[iy_]; j += dj)
  //     for (int k = centerBoxLo[iz_]; k <= centerBoxHi[iz_]; k += dk) {
  //       if (gridInfo.get_status(i, j, k) == GridInfo::iPicOn_) {

  //         if (i < minBoxLo[ix_])
  //           minBoxLo[ix_] = i;
  //         if (j < minBoxLo[iy_])
  //           minBoxLo[iy_] = j;
  //         if (k < minBoxLo[iz_])
  //           minBoxLo[iz_] = k;

  //         if (i + di - 1 > minBoxHi[ix_])
  //           minBoxHi[ix_] = i + di - 1;
  //         if (j + dj - 1 > minBoxHi[iy_])
  //           minBoxHi[iy_] = j + dj - 1;
  //         if (k + dk - 1 > minBoxHi[iz_])
  //           minBoxHi[iz_] = k + dk - 1;

  //       }
  //     }

  BoxList bl;
  get_boxlist_from_region(bl, gridInfo, centerBoxLo, centerBoxHi);
  // get_boxlist_from_region(blInActive, gridInfo, centerBoxLo, centerBoxHi);

  // BoxArray batmp(bl);
  // Print() << "batmp = " << batmp << std::endl;

  // BoxList blActive;
  // BoxArray ba1(blInActive);
  // blActive = ba1.complementIn(centerBox);
  // for (int i = 0; i < 10; i++) {
  //   blActive.simplify();
  // }
  // Print() << "blActive" << blActive << std::endl;

  // for (int i = 0; i < 3; i++) {
  //   bl.simplify();
  // }
  // BoxArray batmp(bl); 
  // batmp.maxSize(maxBlockSize); 
  // bl.clear(); 
  // bl = BoxList(batmp); 
  // for (int i = 0; i < 3; i++) {
  //   bl.simplify();
  // }
  


  BoxArray baNew(bl);
  // add_boxes_to_BoxArray(baNew, vecBox);
  baNew.maxSize(maxBlockSize);
  Print() << "Total box # = " << baNew.size() << std::endl;
  baPicOld = baNew;

  return baNew;
}
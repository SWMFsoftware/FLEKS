#include <AMReX_MultiFabUtil.H>

#include "Domain.h"

void Domain::init() {
  define_domain();
  init_field();
  init_particles();
}
//---------------------------------------------------------

void Domain::define_domain() {
  iProc = ParallelDescriptor::MyProc();
  {
    //---- Geometry initialization -------
    nGst = 2;

    nCell[ix_] = 8;
    nCell[iy_] = 8;
    if (nDim > 2)
      nCell[iz_] = 8;

    nCellBlockMax = 4;

    for (auto& x : periodicity)
      x = 1;

    for (int i = 0; i < nDim; i++) {
      centerBoxLo[i] = 0;
      centerBoxHi[i] = nCell[i] - 1;

      boxRange.setLo(i, 0);
      boxRange.setHi(i, 8);
    }

    centerBox.setSmall(centerBoxLo);
    centerBox.setBig(centerBoxHi);

    coord = 0; // Cartesian

    geom.define(centerBox, &boxRange, coord, periodicity);

    centerBA.define(centerBox);
    centerBA.maxSize(nCellBlockMax);

    dm.define(centerBA);

    nodeBA = convert(centerBA, IntVect{ AMREX_D_DECL(1, 1, 1) });
  }

  {
    // EM field
    nodeE.define(nodeBA, dm, 3, nGst);
    nodeE.setVal(0.0);
    nodeB.define(nodeBA, dm, 3, nGst);
    nodeB.setVal(0.0);

    centerB.define(centerBA, dm, 3, nGst);
    centerB.setVal(0.0);
  }

  {
    // Plasma
    nSpecies = 2;
    nodePlasma.resize(nSpecies);
    for (auto& pl : nodePlasma) {
      pl.define(nodeBA, dm, nMoments, nGst);
      pl.setVal(0.0);
    }

    // Only 1 ghost cell layer is needed!
    nodeMMatrix.define(nodeBA, dm, 27 * 9, 1);
    nodeMMatrix.setVal(0.0);

    partVect.resize(nSpecies);

  }

  //   AllPrint() << "iproc =  " << ParallelDescriptor::MyProc()
  //              << " local fab size = " << nodeE.local_size() << " dm = " <<
  //              dm
  //              << "\n";
}
//---------------------------------------------------------

void Domain::init_field() {
  const Real* dx = geom.CellSize();

  auto init_E = [&dx](const Box& box, Array4<Real> const& E) {
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          const Real x = dx[ix_] * i;
          const Real y = dx[iy_] * j;
          const Real z = dx[iz_] * k;

          E(i, j, k, ix_) = x;
          E(i, j, k, iy_) = y;
          E(i, j, k, iz_) = z;
          Print() << " dx = " << dx[ix_] << " y = " << y
                  << " Ex = " << E(i, j, k, ix_) << std::endl;
        }
  };
  //--------------------------

  auto init_B = [&dx](const Box& box, Array4<Real> const& B) {
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          const Real x = dx[ix_] * i;
          const Real y = dx[iy_] * j;
          const Real z = dx[iz_] * k;

          B(i, j, k, ix_) = 10 * x;
          B(i, j, k, iy_) = 10 * y;
          B(i, j, k, iz_) = 10 * z;
        }
  };
  //--------------------------

  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) // Loop over grids
  {
    FArrayBox& fab = nodeE[mfi];
    init_E(mfi.validbox(), fab.array());
  }

  for (MFIter mfi(nodeB); mfi.isValid(); ++mfi) // Loop over grids
  {
    FArrayBox& fab = nodeB[mfi];
    init_B(mfi.validbox(), fab.array());
  }

  // Interpolate from node to cell center.
  average_node_to_cellcenter(centerB, 0, nodeB, 0, 3);

  nodeE.FillBoundary(geom.periodicity());
  nodeB.FillBoundary(geom.periodicity());
  centerB.FillBoundary(geom.periodicity());
  for (auto& pl : nodePlasma) {
    pl.FillBoundary(geom.periodicity());
  }

  auto print_B = [&dx](const Box& box, Array4<Real> const& B) {
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          AllPrint() << " i = " << i << " j = " << j << " k = " << k
                     << " Bx = " << B(i, j, k, ix_)
                     << " By = " << B(i, j, k, iy_)
                     << " Bz = " << B(i, j, k, iz_) << std::endl;
        }
  };

  for (MFIter mfi(centerB); mfi.isValid(); ++mfi) // Loop over grids
  {
    FArrayBox& fab = centerB[mfi];
    print_B(mfi.fabbox(), fab.array());
  }
}
//---------------------------------------------------------


void Domain::init_particles() {

for(auto& parts: partVect){
  //parts.ini
}
}

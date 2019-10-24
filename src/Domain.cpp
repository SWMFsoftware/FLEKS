#include <AMReX_MultiFabUtil.H>

#include "Domain.h"

void Domain::init(Real timeIn, const std::string& paramString, int* paramInt,
                  double* gridDim, double* paramReal, int iDomain) {

  fluidInterface.set_myrank(ParallelDescriptor::MyProc());
  fluidInterface.set_nProcs(ParallelDescriptor::NProcs());

  std::stringstream* ss = nullptr;
  fluidInterface.ReadFromGMinit(paramInt, gridDim, paramReal, ss);

  fluidInterface.readParam = paramString;

  timeNowSI = timeIn;
  timeNow = timeIn * fluidInterface.getSi2NoT();
  dt = 1;
  dtSI = dt * fluidInterface.getNo2SiT();

  read_param();

  fluidInterface.PrintFluidPicInterface();

  define_domain();

  Print() << "Domain::init() end" << std::endl;
}
//---------------------------------------------------------

void Domain::read_param() {
  qom = new double[1];
  qom[0] = -777;

  std::string command;
  ReadParam& readParam = fluidInterface.readParam;
  readParam.set_verbose(iProc == 0);
  while (readParam.get_next_command(command)) {
    if (command == "#NSYNC") {

    } else if (command == "#PARTICLES") {
      npcelx = new int[1];
      npcely = new int[1];
      npcelz = new int[1];
      readParam.read_var("npcelx", npcelx[0]);
      readParam.read_var("npcely", npcely[0]);
      readParam.read_var("npcelz", npcelz[0]);
    } else if (command == "#ELECTRON") {
      // iones info comes from BATSRUS
      readParam.read_var("qom", qom[0]);
    }
  }

  // Passing qom npcelx... into this function and re-allocating memory
  // for them is extremely ugly!! --Yuxi
  fluidInterface.fixPARAM(qom, npcelx, npcely, npcelz, &nSpecies);

  // Print()<<" nSpecies = "<<nSpecies<<std::endl;
}

void Domain::set_ic() {
  init_field();
  init_particles();
}

void Domain::define_domain() {
  iProc = ParallelDescriptor::MyProc();
  {
    //---- Geometry initialization -------
    nGst = 1;

    nCell[ix_] = fluidInterface.getFluidNxc();
    nCell[iy_] = fluidInterface.getFluidNyc();
    nCell[iz_] = fluidInterface.getFluidNzc();

    nCellBlockMax = 8;

    for (auto& x : periodicity)
      x = 1;

    for (int i = 0; i < nDim; i++) {
      centerBoxLo[i] = 0;
      centerBoxHi[i] = nCell[i] - 1;

      boxRange.setLo(i, fluidInterface.getphyMin(i));
      boxRange.setHi(i, fluidInterface.getphyMax(i));
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
    // nSpecies = 2;
    iTot = nSpecies;
    // The last one is the sum of all species.
    nodePlasma.resize(nSpecies + 1);
    for (auto& pl : nodePlasma) {
      pl.define(nodeBA, dm, nMoments, nGst);
      pl.setVal(0.0);
    }

    // Only 1 ghost cell layer is needed!
    nodeMMatrix.define(nodeBA, dm, 27 * 9, 1);
    nodeMMatrix.setVal(0.0);

    // partVect.resize(nSpecies);
  }

  Print() << " centerBox = " << centerBox << " boxRange = " << boxRange
          << std::endl;

  //   AllPrint() << "iproc =  " << ParallelDescriptor::MyProc()
  //              << " local fab size = " << nodeE.local_size() << " dm = " <<
  //              dm
  //              << "\n";
}
//---------------------------------------------------------

void Domain::init_field() {
  const Real* dx = geom.CellSize();

  int iBlock = 0;
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) // Loop over grids
  {
    FArrayBox& fab = nodeE[mfi];
    const Box& box = mfi.validbox();
    const Array4<Real>& E = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          E(i, j, k, ix_) =
              fluidInterface.getEx(iBlock, i - lo.x, j - lo.y, k - lo.z);
          E(i, j, k, iy_) =
              fluidInterface.getEy(iBlock, i - lo.x, j - lo.y, k - lo.z);
          E(i, j, k, iz_) =
              fluidInterface.getEz(iBlock, i - lo.x, j - lo.y, k - lo.z);
        }

    iBlock++;
  }

  iBlock = 0;
  for (MFIter mfi(nodeB); mfi.isValid(); ++mfi) // Loop over grids
  {
    FArrayBox& fab = nodeB[mfi];

    const Box& box = mfi.validbox();
    const Array4<Real>& B = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          B(i, j, k, ix_) =
              fluidInterface.getBx(iBlock, i - lo.x, j - lo.y, k - lo.z);
          B(i, j, k, iy_) =
              fluidInterface.getBy(iBlock, i - lo.x, j - lo.y, k - lo.z);
          B(i, j, k, iz_) =
              fluidInterface.getBz(iBlock, i - lo.x, j - lo.y, k - lo.z);
        }

    iBlock++;
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
  for (int i = 0; i < nSpecies; i++) {
    IntVect nPartPerCell = { npcelx[i], npcely[i], npcelz[i] };
    auto ptr = std::make_unique<Particles>(
        geom, dm, centerBA, i, fluidInterface.getQiSpecies(i),
        fluidInterface.getMiSpecies(i), nPartPerCell);
    parts.push_back(std::move(ptr));
  }

  for (auto& pts : parts) {
    pts->add_particles_domain(fluidInterface);
  }

  sum_moments();
  Print() << " parts size = " << parts.size() << std::endl;

  // for(auto& parts: partVect){
  // parts.ini
  //}
}

void Domain::sum_moments() {
  nodePlasma[nSpecies].setVal(0.0);
  for (int i = 0; i < nSpecies; i++) {
    parts[i]->sum_moments(nodePlasma[i], dt);
    MultiFab::Add(nodePlasma[nSpecies], nodePlasma[i], 0, 0, nMoments, 0);
  }

    for (MFIter mfi(nodePlasma[0]); mfi.isValid(); ++mfi) // Loop over grids
  {
    FArrayBox& fab = nodePlasma[0][mfi];
  

    const Box& box = mfi.validbox();
    const Array4<Real>& data = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) 
        for(int iVar = 0; iVar<fab.nComp(); iVar++)
        {
          Print()<<" i = "<<i<<" j = "<<j<<" k = "<<k
          <<" iVar = "<<iVar<<" val = "<<data(i,j,k,iVar)<<std::endl;

        }
  }

  // sum to total here
}

void Domain::update() {
  timeNow += dt;
  timeNowSI += dtSI;
}

void Domain::set_state_var(double* data, int* index) {
  std::string nameFunc = "Domain::set_state_var";
  Print() << nameFunc << " begin" << std::endl;

  int nBlockLocal = nodeE.local_size();
  fluidInterface.BlockMin_BD.clear();
  fluidInterface.BlockMax_BD.clear();
  fluidInterface.CellSize_BD.clear();

  if (nBlockLocal > 0) {
    const int nDimMax = 3;
    fluidInterface.BlockMin_BD.init(nBlockLocal, nDimMax);
    fluidInterface.BlockMax_BD.init(nBlockLocal, nDimMax);
    fluidInterface.CellSize_BD.init(nBlockLocal, nDimMax);

    const Real* dx = geom.CellSize();

    Print() << " 1dxz = " << dx[iz_] << std::endl;
    int iBlock = 0;
    for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto lo = lbound(box);
      const auto hi = ubound(box);

      double xMin[3] = { lo.x * dx[ix_], lo.y * dx[iy_], lo.z * dx[iz_] };
      double xMax[3] = { hi.x * dx[ix_], hi.y * dx[iy_], hi.z * dx[iz_] };

      for (int iDim = 0; iDim < nDimMax; iDim++) {
        fluidInterface.BlockMin_BD(iBlock, iDim) = xMin[iDim];
        fluidInterface.BlockMax_BD(iBlock, iDim) = xMax[iDim];
        fluidInterface.CellSize_BD(iBlock, iDim) = dx[iDim];

        fluidInterface.nG_D[iDim] = nGst;
      }
      iBlock++;
    }

    MFIter mfi(nodeE);
    const Box& box = mfi.fabbox();
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    fluidInterface.set_State_BGV(nBlockLocal, hi.x - lo.x + 1, hi.y - lo.y + 1,
                                 hi.z - lo.z + 1, data, index);
  }
  Print() << nameFunc << " end" << std::endl;
}

int Domain::get_grid_nodes_number() {

  int nPoint = 0;
  int nBlock = nodeE.local_size();
  if (nBlock > 0) {
    MFIter mfi(nodeE);
    BL_ASSERT(mfi.isValid());
    const auto lo = lbound(mfi.fabbox());
    const auto hi = ubound(mfi.fabbox());

    nPoint = (hi.x - lo.x + 1) * (hi.y - lo.y + 1);

    if (fluidInterface.getnDim() > 2)
      nPoint *= (hi.z - lo.z + 1);
  }
  nPoint *= nBlock;
  return nPoint;
}

void Domain::get_grid(double* pos_DI) {
  std::string nameFunc = "Domain::get_grid";
  Print() << nameFunc << " begin" << std::endl;
  const Real* dx = geom.CellSize();
  int iBlock = 0;
  int nCount = 0;
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
    const Box& box = mfi.fabbox();
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    int kMin = 0, kMax = 0;
    if (fluidInterface.getnDim() > 2) {
      kMin = lo.z;
      kMax = hi.z;
    }
    for (int i = lo.x; i <= hi.x; ++i)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int k = kMin; k <= kMax; ++k) {
          pos_DI[nCount++] = i * dx[ix_] + boxRange.lo(ix_);
          pos_DI[nCount++] = j * dx[iy_] + boxRange.lo(iy_);
          if (fluidInterface.getnDim() > 2)
            pos_DI[nCount++] = k * dx[iz_] + boxRange.lo(iz_);
        }
  }
  Print() << nameFunc << " end" << std::endl;
}
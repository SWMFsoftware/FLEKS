#include <math.h>

#include <AMReX_MultiFabUtil.H>

#include "../../srcInterface/multi_ipic3d_domain.h"
#include "Domain.h"
#include "LinearSolver.h"
#include "Timing_c.h"
#include "Utility.h"

using namespace amrex;

void Domain::init(Real timeIn, const std::string& paramString, int* paramInt,
                  double* gridDim, double* paramReal, int iDomain) {

  fluidInterface.set_myrank(ParallelDescriptor::MyProc());
  fluidInterface.set_nProcs(ParallelDescriptor::NProcs());

  std::stringstream* ss = nullptr;

  fluidInterface.ReadFromGMinit(paramInt, gridDim, paramReal, ss);
  fluidInterface.readParam = paramString;

  int iCycle = 0;
  fluidInterface.setCycle(iCycle);

  iProc = ParallelDescriptor::MyProc();
  read_param();

  fluidInterface.PrintFluidPicInterface();

  define_domain();

  init_time_ctr();
}
//---------------------------------------------------------

void Domain::init_time_ctr() {
  tc.set_si2no(fluidInterface.getSi2NoT());

  { //----------Init plot data------------------------

    //------ Scalar parameters.----------
    std::vector<std::string> scalarName_I;
    std::vector<double> scalarVar_I;
    std::string ms = "mS", qs = "qS";
    const int nS = fluidInterface.get_nS();
    for (int i = 0; i < nS; ++i) {
      scalarName_I.push_back(ms + std::to_string(i));
      scalarName_I.push_back(qs + std::to_string(i));
      scalarVar_I.push_back(fluidInterface.getMiSpecies(i));
      scalarVar_I.push_back(fluidInterface.getQiSpecies(i));
    }
    scalarName_I.push_back("cLight");
    scalarVar_I.push_back(fluidInterface.getcLightSI());
    scalarName_I.push_back("rPlanet");
    scalarVar_I.push_back(fluidInterface.getrPlanet());
    //-------------------------------------

    for (auto& plot : tc.plots) {
      auto& writer = plot.writer;

      // Pass information to writers.
      writer.set_rank(ParallelDescriptor::MyProc());
      writer.set_nProcs(ParallelDescriptor::NProcs());
      writer.set_nDim(fluidInterface.getnDim());
      writer.set_iRegion(0);
      writer.set_domainMin_D(
          { { boxRange.lo(ix_), boxRange.lo(iy_), boxRange.lo(iz_) } });

      writer.set_domainMax_D(
          { { boxRange.hi(ix_), boxRange.hi(iy_), boxRange.hi(iz_) } });

      const Real* dx = geom.CellSize();
      writer.set_dx_D({ { dx[ix_], dx[iy_], dx[iz_] } });

      writer.set_axisOrigin_D({ { 0, 0, 0 } });
      writer.set_nSpecies(nS);
      writer.set_units(fluidInterface.getNo2SiL(), fluidInterface.getNo2SiV(),
                       fluidInterface.getNo2SiB(), fluidInterface.getNo2SiRho(),
                       fluidInterface.getNo2SiP(), fluidInterface.getNo2SiJ(),
                       fluidInterface.getrPlanet());
      writer.set_scalarValue_I(scalarVar_I);
      writer.set_scalarName_I(scalarName_I);
      //--------------------------------------------------
      writer.init();
      writer.print();
    }
  }
}

//---------------------------------------------------------

void Domain::read_param() {

  nCellBlockMax = 8;

  qom = new double[1];
  qom[0] = -777;

  tc.set_dt(1);

  std::string command;
  ReadParam& readParam = fluidInterface.readParam;
  readParam.set_verbose(iProc == 0);
  while (readParam.get_next_command(command)) {
    if (command == "#MAXBLOCKSIZE") {
      // The block size in each direction can not larger than nCellBlockMax.
      readParam.read_var("nCellBlockMax", nCellBlockMax);
    } else if (command == "#TIMECONTROL") {
      Real dtSI;
      readParam.read_var("dtSI", dtSI);
      tc.set_dt_si(dtSI);
    } else if (command == "#PARTICLES") {
      npcelx = new int[1];
      npcely = new int[1];
      npcelz = new int[1];
      readParam.read_var("npcelx", npcelx[0]);
      readParam.read_var("npcely", npcely[0]);
      readParam.read_var("npcelz", npcelz[0]);
    } else if (command == "#ELECTRON") {
      readParam.read_var("qom", qom[0]);
    } else if (command == "#DISCRETIZE") {
      readParam.read_var("theta", fsolver.theta);
      readParam.read_var("coefDiff", fsolver.coefDiff);
    } else if (command == "#SAVEPLOT") {
      int nPlot;
      readParam.read_var("nPlotFile", nPlot);

      for (int iPlot = 0; iPlot < nPlot; iPlot++) {

        std::string plotString;
        readParam.read_var("plotString", plotString);
        {
          std::string::size_type pos = plotString.find_first_not_of(' ');
          if (pos != std::string::npos)
            plotString.erase(0, pos);
        }

        int dnSave;
        readParam.read_var("dnSavePlot", dnSave);

        Real dtSave;
        readParam.read_var("dtSavePlot", dtSave);

        int dxSave;
        readParam.read_var("dxSavePlot", dxSave);

        std::string plotVar;
        if (plotString.find("var") != std::string::npos) {
          readParam.read_var("plotVar", plotVar);
        }

        PlotCtr pcTmp(&tc, iPlot, dtSave, dnSave, plotString, dxSave, plotVar);
        tc.plots.push_back(pcTmp);
      }
    }
  }

  // Passing qom npcelx... into this function and re-allocating memory
  // for them is extremely ugly!! --Yuxi
  fluidInterface.fixPARAM(qom, npcelx, npcely, npcelz, &nSpecies);
}

void Domain::set_ic() {
  init_field();
  init_particles();

  tc.write_plots(true);
}

void Domain::define_domain() {
  {
    //---- Geometry initialization -------
    nGst = 1;

    nCell[ix_] = fluidInterface.getFluidNxc();
    nCell[iy_] = fluidInterface.getFluidNyc();
    nCell[iz_] = fluidInterface.getFluidNzc();

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
    nodeEth.define(nodeBA, dm, 3, nGst);
    nodeEth.setVal(0.0);

    nodeB.define(nodeBA, dm, 3, nGst);
    nodeB.setVal(0.0);

    centerB.define(centerBA, dm, 3, nGst);
    centerB.setVal(0.0);

    centerDivE.define(centerBA, dm, 1, nGst);
    centerDivE.setVal(0.0);

    centerPhi.define(centerBA, dm, 1, nGst);
    centerPhi.setVal(0.0);

    tempNode3.define(nodeBA, dm, 3, nGst);
    tempNode3.setVal(0.0);

    tempCenter3.define(centerBA, dm, 3, nGst);
    tempCenter3.setVal(0.0);

    tempCenter1.define(centerBA, dm, 1, nGst);
    tempCenter1.setVal(0.0);

    tempCenter1_1.define(centerBA, dm, 1, nGst);
    tempCenter1_1.setVal(0.0);
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
    nodeMM.define(nodeBA, dm, 1, 1);
    const RealMM mm0(0.0);
    nodeMM.setVal(mm0);

    //-------divE correction----------------
    centerNetChargeOld.define(centerBA, dm, 1, 1);
    centerNetChargeOld.setVal(0.0);
    centerNetChargeN.define(centerBA, dm, 1, 1);
    centerNetChargeN.setVal(0.0);
    centerNetChargeNew.define(centerBA, dm, 1, 1);
    centerNetChargeNew.setVal(0.0);

    centerMM.define(centerBA, dm, 1, 1);
    centerMM.setVal(0.0);
    //--------------------------------------
  }

  Print() << "Domain:: Domain range = " << boxRange << std::endl;
  Print() << "Domain:: Total block #  = " << nodeBA.size() << std::endl;
}
//---------------------------------------------------------

void Domain::init_field() {
  const Real* dx = geom.CellSize();

  int iBlock = 0;
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) // Loop over grids
  {
    FArrayBox& fab = nodeE[mfi];
    const Box& box = mfi.fabbox();
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

    const Box& box = mfi.fabbox();
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
  sum_to_center(false);
}

void Domain::particle_mover() {
  std::string nameFunc = "Domain::mover";
  timing_start(nameFunc);

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->mover(nodeEth, nodeB, tc.get_dt());
  }

  timing_stop(nameFunc);
}

void Domain::sum_moments() {
  std::string nameFunc = "Domain::sum_moments";
  timing_start(nameFunc);

  nodePlasma[nSpecies].setVal(0.0);
  const RealMM mm0(0.0);
  nodeMM.setVal(mm0);

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->sum_moments(nodePlasma[i], nodeMM, nodeB, tc.get_dt());
    MultiFab::Add(nodePlasma[nSpecies], nodePlasma[i], 0, 0, nMoments, 0);
  }

  nodeMM.SumBoundary(geom.periodicity());
  nodeMM.FillBoundary(geom.periodicity());

  timing_stop(nameFunc);
}

void Domain::divE_correction() {
  std::string nameFunc = "Domain::divE_correction";
  timing_start(nameFunc);

  nSolveCenter =
      get_fab_grid_points_number(centerDivE) * centerDivE.local_size();

  for (int iIter = 0; iIter < 3; iIter++) {
    sum_to_center(true);

    Print() << "----------- div(E) correction ----------" << std::endl;
    calculate_phi(matvec_divE_accurate);

    divE_correct_particle_position();
  }
  // DO CORRECTION
  sum_to_center(false);
  timing_stop(nameFunc);
}

void Domain::divE_correct_particle_position() {

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->divE_correct_position(centerPhi);
  }
}

void Domain::calculate_phi(MATVEC fMatvec, Real tol, int nIter) {
  double xLeft[nSolveCenter];
  double b[nSolveCenter];

  for (int i = 0; i < nSolveCenter; i++) {
    xLeft[i] = 0;
  }

  div_node_to_center(nodeE, tempCenter1, geom.InvCellSize());

  // lap(phi) = (divE - rhoc*(4pi))/rhoTheta, where,
  // rhoc = (1-rhoTheta)*rhoc(n+1/2) + rhoTheta*rhoc(n+3/2)
  // The correction density Delt_rhoc(n+3/2) = 4*pi*lap(phi)

  MultiFab::LinComb(tempCenter1, 1.0 / rhoTheta, tempCenter1, 0,
                    -fourPI / rhoTheta, centerNetChargeN, 0, 0,
                    tempCenter1.nComp(), 0);

  tempCenter1.FillBoundary(geom.periodicity());

  convert_3d_to_1d(tempCenter1, b, geom);

  {
    linear_solver_matvec_c = fMatvec;

    int nVarSolve = 1;
    int nDimIn = 3;

    MFIter mfi(tempCenter1);
    auto lo = lbound(mfi.validbox());
    auto hi = ubound(mfi.validbox());

    int nI = hi.x - lo.x + 1;
    int nJ = hi.y - lo.y + 1;
    int nK = hi.z - lo.z + 1;
    int nBlock = tempCenter1.local_size();
    MPI_Fint iComm = MPI_Comm_c2f(ParallelDescriptor::Communicator());
    double precond_matrix_II[1][1];
    precond_matrix_II[0][0] = 0;
    // parameter to choose preconditioner types
    // 0:No precondition; 1: BILU; 2:DILU;
    //[-1,0): MBILU;
    double PrecondParam = 0;
    int lTest = ParallelDescriptor::MyProc() == 0;

    linear_solver_wrapper("GMRES", &tol, &nIter, &nVarSolve, &nDimIn, &nI, &nJ,
                          &nK, &nBlock, &iComm, b, xLeft, &PrecondParam,
                          precond_matrix_II[0], &lTest);
  }

  convert_1d_to_3d(xLeft, centerPhi, geom);
  centerPhi.FillBoundary(geom.periodicity());
}

void Domain::divE_accurate_matvec(double* vecIn, double* vecOut) {

  convert_1d_to_3d(vecIn, tempCenter1, geom);
  tempCenter1.FillBoundary(geom.periodicity());

  tempCenter1_1.setVal(0.0);
  for (amrex::MFIter mfi(tempCenter1); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const amrex::Array4<amrex::Real>& lArr = tempCenter1_1[mfi].array();
    const amrex::Array4<amrex::Real const>& rArr = tempCenter1[mfi].array();
    const amrex::Array4<RealCMM>& mmArr = centerMM[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i)

          for (int i2 = i - 1; i2 <= i + 1; i2++)
            for (int j2 = j - 1; j2 <= j + 1; j2++)
              for (int k2 = k - 1; k2 <= k + 1; k2++) {
                const int gp = (i2 - i + 1) * 9 + (j2 - j + 1) * 3 + k2 - k + 1;
                lArr(i, j, k) += rArr(i2, j2, k2) * mmArr(i, j, k).data[gp];
              }
  }
  tempCenter1_1.mult(fourPI * fourPI);
  convert_3d_to_1d(tempCenter1_1, vecOut, geom);
}

void Domain::sum_to_center(bool isBeforeCorrection) {
  centerNetChargeNew.setVal(0.0);

  const RealCMM mm0(0.0);
  centerMM.setVal(mm0);

  bool doNetChargeOnly = !isBeforeCorrection;

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->sum_to_center(centerNetChargeNew, centerMM, doNetChargeOnly);
  }

  if (!doNetChargeOnly) {
    centerMM.SumBoundary(geom.periodicity());
    centerMM.FillBoundary(geom.periodicity());
  }

  centerNetChargeNew.SumBoundary(geom.periodicity());
  centerNetChargeNew.FillBoundary(geom.periodicity());

  MultiFab::LinComb(centerNetChargeN, 1 - rhoTheta, centerNetChargeOld, 0,
                    rhoTheta, centerNetChargeNew, 0, 0,
                    centerNetChargeN.nComp(), centerNetChargeN.nGrow());

  if (!isBeforeCorrection) {
    MultiFab::Copy(centerNetChargeOld, centerNetChargeNew, 0, 0,
                   centerNetChargeOld.nComp(), centerNetChargeOld.nGrow());
  }

  // centerNetCharge.SumBoundary(geom.periodicity());
  // centerNetCharge.FillBoundary(geom.periodicity());

  // for (MFIter mfi(centerMM); mfi.isValid(); ++mfi) {
  //   const Box& box = mfi.validbox();
  //   auto const& data = centerMM[mfi].array();

  //   const auto lo = lbound(box);
  //   const auto hi = ubound(box);

  //   Real sum=0;
  //   for (int i = lo.x; i <= hi.x; ++i)
  //     for (int j = lo.y; j <= hi.y; ++j)
  //       for (int k = lo.z; k <= hi.z; ++k) {
  //         AllPrint() << " i = " << i << " j = " << j << " k = " << k
  //                    << " centerMM = " << data(i, j, k) << std::endl;
  //                    sum += data(i,j,k).data[0];
  //       }
  //       Print()<<" sum = "<<sum<<std::endl;
  // }

  // centerNetCharge.FillBoundary
}

void Domain::update() {
  std::string nameFunc = "Domain::update";
  timing_start(nameFunc);

  Print() << "================ Begin cycle = " << tc.get_cycle()
          << " at time = " << tc.get_time_si()
          << " (s) ======================" << std::endl;
  update_E();

  particle_mover();
  update_B();

  divE_correction();
  sum_moments();

  tc.update();
  tc.write_plots();

  timing_stop(nameFunc);
}

void Domain::set_state_var(double* data, int* index) {
  std::string nameFunc = "Domain::set_state_var";

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
}

void Domain::update_E() {
  std::string nameFunc = "Domain::update_E";
  timing_start(nameFunc);

  const int nNode = get_fab_grid_points_number(nodeE) * nodeE.local_size();
  nSolveNode = 3 * nNode;

  double rhs[nSolveNode];
  double xLeft[nSolveNode];
  double matvec[nSolveNode];

  for (int i = 0; i < nSolveNode; i++) {
    rhs[i] = 0;
    xLeft[i] = 0;
    matvec[i] = 0;
  }

  update_E_rhs(rhs);

  convert_3d_to_1d(nodeE, xLeft, geom);

  update_E_matvec(xLeft, matvec);

  for (int i = 0; i < nSolveNode; i++) {
    rhs[i] -= matvec[i];
    xLeft[i] = 0;
  }

  {
    linear_solver_matvec_c = matvec_E_solver;

    int nVarSolve = 3;
    int nIter = 200;
    double EFieldTol = 1e-6;
    int nDimIn = 3;

    MFIter mfi(nodeE);
    auto lo = lbound(mfi.validbox());
    auto hi = ubound(mfi.validbox());

    int nI = hi.x - lo.x + 1;
    int nJ = hi.y - lo.y + 1;
    int nK = hi.z - lo.z + 1;
    int nBlock = nodeE.local_size();
    MPI_Fint iComm = MPI_Comm_c2f(ParallelDescriptor::Communicator());
    double precond_matrix_II[1][1];
    precond_matrix_II[0][0] = 0;
    // parameter to choose preconditioner types
    // 0:No precondition; 1: BILU; 2:DILU;
    //[-1,0): MBILU;
    double PrecondParam = 0;
    int lTest = ParallelDescriptor::MyProc() == 0;

    Print() << "-------------- Electric field solver ---------------" << std::endl;
    linear_solver_wrapper("GMRES", &EFieldTol, &nIter, &nVarSolve, &nDimIn, &nI,
                          &nJ, &nK, &nBlock, &iComm, rhs, xLeft, &PrecondParam,
                          precond_matrix_II[0], &lTest);    
  }

  nodeEth.setVal(0.0);
  convert_1d_to_3d(xLeft, nodeEth, geom);
  nodeEth.SumBoundary(geom.periodicity());
  nodeEth.FillBoundary(geom.periodicity());

  MultiFab::Add(nodeEth, nodeE, 0, 0, nodeEth.nComp(), nGst);

  MultiFab::LinComb(nodeE, -(1.0 - fsolver.theta) / fsolver.theta, nodeE, 0,
                    1. / fsolver.theta, nodeEth, 0, 0, nodeE.nComp(), nGst);

  timing_stop(nameFunc);
}

void Domain::update_E_matvec(const double* vecIn, double* vecOut) {
  zero_array(vecOut, nSolveNode);

  MultiFab vecMF(nodeBA, dm, 3, nGst);
  vecMF.setVal(0.0);

  MultiFab imageMF(nodeBA, dm, 3, nGst);
  imageMF.setVal(0.0);

  convert_1d_to_3d(vecIn, vecMF, geom);

  // The right side edges should be filled in.
  vecMF.SumBoundary(geom.periodicity());

  // M*E needs ghost cell information.
  vecMF.FillBoundary(geom.periodicity());

  lap_node_to_node(vecMF, imageMF, dm, geom);

  Real delt2 = pow(fsolver.theta * tc.get_dt(), 2);
  imageMF.mult(-delt2);

  { // grad(divE)
    div_node_to_center(vecMF, centerDivE, geom.InvCellSize());
    if (fsolver.coefDiff > 0) {
      // Calculate cell center E for center-to-center divE.
      average_node_to_cellcenter(tempCenter3, 0, vecMF, 0, 3);
      tempCenter3.FillBoundary(geom.periodicity());
      div_center_to_center(tempCenter3, tempCenter1, geom.InvCellSize());

      MultiFab::LinComb(centerDivE, 1 - fsolver.coefDiff, centerDivE, 0,
                        fsolver.coefDiff, tempCenter1, 0, 0, 1, 0);
    }
    centerDivE.FillBoundary(geom.periodicity());

    grad_center_to_node(centerDivE, tempNode3, geom.InvCellSize());

    tempNode3.mult(delt2);
    MultiFab::Add(imageMF, tempNode3, 0, 0, imageMF.nComp(), 0);
  }

  update_E_M_dot_E(vecMF, tempNode3);
  MultiFab::Add(imageMF, tempNode3, 0, 0, imageMF.nComp(), 0);

  MultiFab::Add(imageMF, vecMF, 0, 0, imageMF.nComp(), 0);

  imageMF.FillBoundary(geom.periodicity());

  convert_3d_to_1d(imageMF, vecOut, geom);
}

void Domain::update_E_M_dot_E(const MultiFab& inMF, MultiFab& outMF) {

  outMF.setVal(0.0);
  Real c0 = fourPI * fsolver.theta * tc.get_dt();
  for (amrex::MFIter mfi(outMF); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const amrex::Array4<amrex::Real const>& inArr = inMF[mfi].array();
    const amrex::Array4<amrex::Real>& ourArr = outMF[mfi].array();
    const amrex::Array4<RealMM>& mmArr = nodeMM[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i)
          for (int k2 = k - 1; k2 <= k + 1; k2++)
            for (int j2 = j - 1; j2 <= j + 1; j2++)
              for (int i2 = i - 1; i2 <= i + 1; i2++) {

                const int gp = (k2 - k + 1) * 9 + (j2 - j + 1) * 3 + i2 - i + 1;
                const int idx0 = gp * 9;

                Real M_I[9];
                for (int ii = 0; ii < 9; ii++) {
                  M_I[ii] = mmArr(i, j, k).data[idx0 + ii];
                }

                const double& vctX =
                    inArr(i2, j2, k2, ix_); // vectX[i2][j2][k2];
                const double& vctY = inArr(i2, j2, k2, iy_);
                const double& vctZ = inArr(i2, j2, k2, iz_);
                ourArr(i, j, k, ix_) +=
                    (vctX * M_I[0] + vctY * M_I[1] + vctZ * M_I[2]) * c0;
                ourArr(i, j, k, iy_) +=
                    (vctX * M_I[3] + vctY * M_I[4] + vctZ * M_I[5]) * c0;
                ourArr(i, j, k, iz_) +=
                    (vctX * M_I[6] + vctY * M_I[7] + vctZ * M_I[8]) * c0;
              }
  }
}

void Domain::update_E_rhs(double* rhs) {
  MultiFab tempNode(nodeBA, dm, 3, nGst);
  tempNode.setVal(0.0);
  MultiFab temp2Node(nodeBA, dm, 3, nGst);
  temp2Node.setVal(0.0);

  const Real* invDx = geom.InvCellSize();
  curl_center_to_node(centerB, tempNode, invDx);

  MultiFab::Saxpy(temp2Node, -fourPI, nodePlasma[iTot], iJhx_, 0,
                  temp2Node.nComp(), 0);

  MultiFab::Add(temp2Node, tempNode, 0, 0, tempNode.nComp(), 0);
  temp2Node.mult(fsolver.theta * tc.get_dt());

  MultiFab::Add(temp2Node, nodeE, 0, 0, nodeE.nComp(), 0);

  convert_3d_to_1d(temp2Node, rhs, geom);
}

void Domain::update_B() {
  std::string nameFunc = "Domain::update_B";
  timing_start(nameFunc);

  MultiFab dB(centerBA, dm, 3, nGst);

  curl_node_to_center(nodeEth, dB, geom.InvCellSize());

  MultiFab::Saxpy(centerB, -tc.get_dt(), dB, 0, 0, centerB.nComp(), 0);
  centerB.FillBoundary(geom.periodicity());

  average_center_to_node(centerB, nodeB);
  nodeB.FillBoundary(geom.periodicity());

  timing_stop(nameFunc);
}



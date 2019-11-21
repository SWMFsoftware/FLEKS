#include <math.h>

#include <AMReX_MultiFabUtil.H>

#include "Domain.h"
#include "LinearSolver.h"
#include "SWMFDomains.h"
#include "Timing_c.h"
#include "Utility.h"

using namespace amrex;

void Domain::init(Real timeIn, const std::string& paramString, int* paramInt,
                  double* gridDim, double* paramReal, int iDomain) {

  fluidInterface.init();
  fluidInterface.receive_info_from_gm(paramInt, gridDim, paramReal,
                                      paramString);

  // Read from PARAM.in
  read_param();

  fluidInterface.PrintFluidPicInterface();

  make_grid();
  make_data();

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
  // The default values shoudl be set in the constructor.

  qom = new double[1];
  qom[0] = -777;

  std::string command;
  ReadParam& readParam = fluidInterface.readParam;
  readParam.set_verbose(ParallelDescriptor::MyProc() == 0);
  while (readParam.get_next_command(command)) {
    if (command == "#MAXBLOCKSIZE") {
      // The block size in each direction can not larger than maxBlockSize.
      int tmp;
      readParam.read_var("maxBlockSize", tmp);
      set_maxBlockSize(tmp);
    } else if (command == "#TIMECONTROL") {
      Real dtSI;
      readParam.read_var("dtSI", dtSI);
      tc.set_dt_si(dtSI);
    } else if (command == "#DIVE") {
      readParam.read_var("doCorrectDivE", doCorrectDivE);
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
    } else if (command == "#PERIODICITY") {
      for (int i = 0; i < nDimMax; i++) {
        bool isPeriodic;
        readParam.read_var("isPeriodic", isPeriodic);
        set_periodicity(i, isPeriodic);
      }
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

      //--------- The commands below exist in restart.H only --------
    } else if (command == "#RESTART") {
      readParam.read_var("doRestart", doRestart);
    } else if (command == "#NSTEP") {
      int nStep;
      readParam.read_var("nStep", nStep);
      tc.set_cycle(nStep);
    } else if (command == "#TIMESIMULATION") {
      Real time;
      readParam.read_var("time", time);
      tc.set_time_si(time);
    } else if (command == "#GEOMETRY") {
      RealBox bxTmp;
      for (int i = 0; i < nDimMax; ++i) {
        Real lo, hi;
        readParam.read_var("min", lo);
        readParam.read_var("max", hi);
        bxTmp.setLo(i, lo);
        bxTmp.setHi(i, hi);
      }
      set_boxRange(bxTmp);
    } else if (command == "#NCELL") {
      IntVect tmp;
      for (int i = 0; i < nDimMax; ++i) {
        readParam.read_var("nCell", tmp[i]);
      }
      set_nCell(tmp);
    }
    //--------- The commands above exist in restart.H only --------
  }

  // Passing qom npcelx... into this function and re-allocating memory
  // for them is extremely ugly!! --Yuxi
  fluidInterface.fixPARAM(qom, npcelx, npcely, npcelz, &nSpecies);
}

void Domain::set_ic() {

  if (doRestart) {
    read_restart();
  } else {
    set_ic_field();
    set_ic_particles();
  }

  tc.write_plots(true);
}

void Domain::make_grid() {
  set_nGst(1);

  if (!doRestart) {
    IntVect nCellTmp;
    RealBox boxRangeTmp;

    // If restart, these variables read from restart.H
    nCellTmp[ix_] = fluidInterface.getFluidNxc();
    nCellTmp[iy_] = fluidInterface.getFluidNyc();
    nCellTmp[iz_] = fluidInterface.getFluidNzc();

    for (int i = 0; i < nDim; i++) {
      boxRangeTmp.setLo(i, fluidInterface.getphyMin(i));
      boxRangeTmp.setHi(i, fluidInterface.getphyMax(i));
    }
    set_nCell(nCellTmp);
    set_boxRange(boxRangeTmp);
  }

  DomainGrid::init();

  // Also make grid for the fluid interface.
  fluidInterface.make_grid(dm, geom, centerBA, nodeBA, nGst);
}

void Domain::make_data() {

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

    // Create particle containers.
    for (int i = 0; i < nSpecies; i++) {
      IntVect nPartPerCell = { npcelx[i], npcely[i], npcelz[i] };
      auto ptr = std::make_unique<Particles>(
          geom, dm, centerBA, &tc, i, fluidInterface.getQiSpecies(i),
          fluidInterface.getMiSpecies(i), nPartPerCell);
      parts.push_back(std::move(ptr));
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

  {

    int nGrid = get_local_node_or_cell_number(nodeE);
    double tol = 1e-6;
    int nIter = 200;
    eSolver.init(nGrid, nDimMax, nDimMax, tol, nIter, matvec_E_solver);
  }

  {
    int nGrid = get_local_node_or_cell_number(centerDivE);
    double tol = 1e-2;
    int nIter = 20;
    divESolver.init(nGrid, 1, nDimMax, tol, nIter, matvec_divE_accurate);
  }
}
//---------------------------------------------------------

void Domain::set_ic_field() {
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) // Loop over grids
  {
    FArrayBox& fab = nodeE[mfi];
    const Box& box = mfi.fabbox();
    const Array4<Real>& arrE = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          arrE(i, j, k, ix_) = fluidInterface.get_ex(mfi, i, j, k);
          arrE(i, j, k, iy_) = fluidInterface.get_ey(mfi, i, j, k);
          arrE(i, j, k, iz_) = fluidInterface.get_ez(mfi, i, j, k);
        }
  }

  MultiFab::Copy(nodeB, fluidInterface.get_nodeFluid(), fluidInterface.iBx, 0,
                 nodeB.nComp(), nodeB.nGrow());

  // Interpolate from node to cell center.
  average_node_to_cellcenter(centerB, 0, nodeB, 0, centerB.nComp(),
                             centerB.nGrow());
}
//---------------------------------------------------------

void Domain::set_ic_particles() {
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
    parts[i]->inject_particles_at_boundary(fluidInterface);
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

  for (int iIter = 0; iIter < 3; iIter++) {
    sum_to_center(true);

    Print() << "\n----------- div(E) correction at iter " << iIter
            << "----------" << std::endl;
    calculate_phi(divESolver);

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

void Domain::calculate_phi(LinearSolver& solver) {
  solver.reset(get_local_node_or_cell_number(centerDivE));

  div_node_to_center(nodeE, tempCenter1, geom.InvCellSize());

  MultiFab::LinComb(tempCenter1, 1.0 / rhoTheta, tempCenter1, 0,
                    -fourPI / rhoTheta, centerNetChargeN, 0, 0,
                    tempCenter1.nComp(), 0);

  tempCenter1.FillBoundary(geom.periodicity());

  convert_3d_to_1d(tempCenter1, solver.rhs, geom);

  solver.solve();

  convert_1d_to_3d(solver.xLeft, centerPhi, geom);
  centerPhi.FillBoundary(geom.periodicity());
}

void Domain::divE_accurate_matvec(double* vecIn, double* vecOut) {
  zero_array(vecOut, divESolver.get_nSolve());

  convert_1d_to_3d(vecIn, tempCenter1, geom);
  tempCenter1.FillBoundary(geom.periodicity());

  apply_float_boundary(tempCenter1, geom, 0, tempCenter1.nComp());

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

  apply_external_BC(centerNetChargeNew, 0, centerNetChargeNew.nComp(),
                    &Domain::get_zero);

  MultiFab::LinComb(centerNetChargeN, 1 - rhoTheta, centerNetChargeOld, 0,
                    rhoTheta, centerNetChargeNew, 0, 0,
                    centerNetChargeN.nComp(), centerNetChargeN.nGrow());

  if (!isBeforeCorrection) {
    MultiFab::Copy(centerNetChargeOld, centerNetChargeNew, 0, 0,
                   centerNetChargeOld.nComp(), centerNetChargeOld.nGrow());
  }
}

void Domain::update() {
  std::string nameFunc = "Domain::update";
  timing_start(nameFunc);

  Print() << "\n================ Begin cycle = " << tc.get_cycle()
          << " at time = " << tc.get_time_si()
          << " (s) ======================" << std::endl;
  update_E();

  particle_mover();
  update_B();

  if (doCorrectDivE) {
    divE_correction();
  }

  sum_moments();

  tc.update();
  tc.write_plots();

  timing_stop(nameFunc);
}

void Domain::update_E() {
  std::string nameFunc = "Domain::update_E";
  timing_start(nameFunc);

  eSolver.reset(get_local_node_or_cell_number(nodeE));

  update_E_rhs(eSolver.rhs);

  convert_3d_to_1d(nodeE, eSolver.xLeft, geom);

  update_E_matvec(eSolver.xLeft, eSolver.matvec, false);

  for (int i = 0; i < eSolver.get_nSolve(); i++) {
    eSolver.rhs[i] -= eSolver.matvec[i];
    eSolver.xLeft[i] = 0;
  }

  Print() << "\n----------- E solver ------------------" << std::endl;
  eSolver.solve();

  nodeEth.setVal(0.0);
  convert_1d_to_3d(eSolver.xLeft, nodeEth, geom);
  nodeEth.SumBoundary(geom.periodicity());
  nodeEth.FillBoundary(geom.periodicity());

  MultiFab::Add(nodeEth, nodeE, 0, 0, nodeEth.nComp(), nGst);

  MultiFab::LinComb(nodeE, -(1.0 - fsolver.theta) / fsolver.theta, nodeE, 0,
                    1. / fsolver.theta, nodeEth, 0, 0, nodeE.nComp(), nGst);

  apply_external_BC(nodeE, 0, nDimMax, &Domain::get_node_E);
  apply_external_BC(nodeEth, 0, nDimMax, &Domain::get_node_E);

  // Apply float BC in order to compare with iPIC3D. It is not right to apply
  // float BC here!!!!!!!--Yuxi
  apply_float_boundary(nodeE, geom, 0, nodeE.nComp());
  apply_float_boundary(nodeEth, geom, 0, nodeEth.nComp());

  timing_stop(nameFunc);
}

void Domain::update_E_matvec(const double* vecIn, double* vecOut,
                             const bool useZeroBC) {
  zero_array(vecOut, eSolver.get_nSolve());

  MultiFab vecMF(nodeBA, dm, 3, nGst);
  vecMF.setVal(0.0);

  MultiFab matvecMF(nodeBA, dm, 3, nGst);
  matvecMF.setVal(0.0);

  convert_1d_to_3d(vecIn, vecMF, geom);

  // The right side edges should be filled in.
  vecMF.SumBoundary(geom.periodicity());

  // M*E needs ghost cell information.
  vecMF.FillBoundary(geom.periodicity());

  if (useZeroBC) {
    // The boundary nodes would not be filled in by convert_1d_3d. So, there is
    // not need to apply zero boundary conditions again here.
  } else {
    apply_external_BC(vecMF, 0, nDimMax, &Domain::get_node_E);
  }

  lap_node_to_node(vecMF, matvecMF, dm, geom);

  Real delt2 = pow(fsolver.theta * tc.get_dt(), 2);
  matvecMF.mult(-delt2);

  { // grad(divE)
    div_node_to_center(vecMF, centerDivE, geom.InvCellSize());
    if (fsolver.coefDiff > 0) {
      // Calculate cell center E for center-to-center divE.
      average_node_to_cellcenter(tempCenter3, 0, vecMF, 0, 3,
                                 tempCenter3.nGrow());

      tempCenter3.FillBoundary(geom.periodicity());

      // It seems not necessary, nor a good idea to apply float BC here. --Yuxi
      apply_float_boundary(tempCenter3, geom, 0, tempCenter3.nComp());

      div_center_to_center(tempCenter3, tempCenter1, geom.InvCellSize());

      MultiFab::LinComb(centerDivE, 1 - fsolver.coefDiff, centerDivE, 0,
                        fsolver.coefDiff, tempCenter1, 0, 0, 1, 0);
    }

    centerDivE.FillBoundary(geom.periodicity());

    grad_center_to_node(centerDivE, tempNode3, geom.InvCellSize());

    tempNode3.mult(delt2);
    MultiFab::Add(matvecMF, tempNode3, 0, 0, matvecMF.nComp(), 0);
  }

  tempNode3.setVal(0);
  update_E_M_dot_E(vecMF, tempNode3);
  MultiFab::Add(matvecMF, tempNode3, 0, 0, matvecMF.nComp(), 0);

  MultiFab::Add(matvecMF, vecMF, 0, 0, matvecMF.nComp(), 0);

  convert_3d_to_1d(matvecMF, vecOut, geom);
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

  apply_external_BC(centerB, 0, centerB.nComp(), &Domain::get_center_B);
  apply_external_BC(nodeB, 0, nodeB.nComp(), &Domain::get_node_B);

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

  apply_external_BC(centerB, 0, centerB.nComp(), &Domain::get_center_B);

  average_center_to_node(centerB, nodeB);
  nodeB.FillBoundary(geom.periodicity());

  apply_external_BC(nodeB, 0, nodeB.nComp(), &Domain::get_node_B);

  timing_stop(nameFunc);
}

void Domain::apply_external_BC(amrex::MultiFab& mf, const int iStart,
                               const int nComp, GETVALUE func) {
  if (geom.isAllPeriodic())
    return;
  if (mf.nGrow() == 0)
    return;

  //! create a grown domain box containing valid + periodic cells
  const Box& domain = geom.Domain();
  Box gdomain = amrex::convert(domain, mf.boxArray().ixType());
  const IntVect& ngrow = mf.nGrowVect();
  for (int i = 0; i < nDimMax; ++i) {
    if (geom.isPeriodic(i)) {
      gdomain.grow(i, ngrow[i]);
    }
  }

  {
    Vector<BCRec> bcDomain(1);
    for (int idim = 0; idim < nDimMax; ++idim) {
      if (geom.isPeriodic(idim)) {
        bcDomain[0].setLo(idim, BCType::int_dir);
        bcDomain[0].setHi(idim, BCType::int_dir);
      } else {
        bcDomain[0].setLo(idim, BCType::ext_dir);
        bcDomain[0].setHi(idim, BCType::ext_dir);
      }
    }

    //"g" means "global".
    const Box& gbx = convert(geom.Domain(), mf.boxArray().ixType());
    const auto glo = gbx.loVect(); // Do not include ghost cells.
    const auto ghi = gbx.hiVect();

    int igMin = glo[ix_], igMax = ghi[ix_];
    int jgMin = glo[iy_], jgMax = ghi[iy_];
    int kgMin = glo[iz_], kgMax = ghi[iz_];

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.fabbox();

      //! if there are cells not in the valid + periodic grown box
      //! we need to fill them here
      //!
      if (!gdomain.contains(bx)) {
        //! Based on bcDomain for the domain, we need to make bcr for this Box
        Vector<BCRec> bcr(1);
        amrex::setBC(bx, domain, 0, 0, 1, bcDomain, bcr);

        amrex::Array4<amrex::Real> const& arr = mf[mfi].array();

        // Include ghost cells.
        const auto lo = bx.loVect();
        const auto hi = bx.hiVect();

        int iMin = lo[ix_], iMax = hi[ix_];
        int jMin = lo[iy_], jMax = hi[iy_];
        int kMin = lo[iz_], kMax = hi[iz_];

        // x left
        if (bcr[0].lo(ix_) == BCType::ext_dir)
          for (int iVar = iStart; iVar < nComp; iVar++)
            for (int k = kMin; k <= kMax; k++)
              for (int j = jMin; j <= jMax; j++)
                for (int i = iMin; i <= igMin - 1 + nVirGst; i++) {
                  arr(i, j, k, iVar) =
                      (this->*func)(mfi, i, j, k, iVar - iStart);
                }

        // x right
        if (bcr[0].hi(ix_) == BCType::ext_dir)
          for (int iVar = iStart; iVar < nComp; iVar++)
            for (int k = kMin; k <= kMax; k++)
              for (int j = jMin; j <= jMax; j++)
                for (int i = igMax + 1 - nVirGst; i <= iMax; i++) {
                  arr(i, j, k, iVar) =
                      (this->*func)(mfi, i, j, k, iVar - iStart);
                }

        // y left
        if (bcr[0].lo(iy_) == BCType::ext_dir)
          for (int iVar = iStart; iVar < nComp; iVar++)
            for (int k = kMin; k <= kMax; k++)
              for (int j = jMin; j <= jgMin - 1 + nVirGst; j++)
                for (int i = iMin; i <= iMax; i++) {
                  arr(i, j, k, iVar) =
                      (this->*func)(mfi, i, j, k, iVar - iStart);
                }

        // y right
        if (bcr[0].hi(iy_) == BCType::ext_dir)
          for (int iVar = iStart; iVar < nComp; iVar++)
            for (int k = kMin; k <= kMax; k++)
              for (int j = jgMax + 1 - nVirGst; j <= jMax; j++)
                for (int i = iMin; i <= iMax; i++) {
                  arr(i, j, k, iVar) =
                      (this->*func)(mfi, i, j, k, iVar - iStart);
                }

        // z left
        if (bcr[0].lo(iz_) == BCType::ext_dir)
          for (int iVar = iStart; iVar < nComp; iVar++)
            for (int k = kMin; k <= kgMin - 1 + nVirGst; k++)
              for (int j = jMin; j <= jMax; j++)
                for (int i = iMin; i <= iMax; i++) {
                  arr(i, j, k, iVar) =
                      (this->*func)(mfi, i, j, k, iVar - iStart);
                }

        // z right
        if (bcr[0].hi(iz_) == BCType::ext_dir)
          for (int iVar = iStart; iVar < nComp; iVar++)
            for (int k = kgMax + 1 - nVirGst; k <= kMax; k++)
              for (int j = jMin; j <= jMax; j++)
                for (int i = iMin; i <= iMax; i++) {
                  arr(i, j, k, iVar) =
                      (this->*func)(mfi, i, j, k, iVar - iStart);
                }
      }
    }
  }
}

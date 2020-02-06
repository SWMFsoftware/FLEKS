#include <math.h>

#include <AMReX_MultiFabUtil.H>

#include "Pic.h"
#include "GridUtility.h"
#include "LinearSolver.h"
#include "SWMFDomains.h"
#include "Timing_c.h"
#include "Utility.h"

using namespace amrex;

void Pic::init(Real timeIn, const std::string& paramString, int* paramInt,
               double* gridDim, double* paramReal,
               std::shared_ptr<FluidInterface>& fluidIn,
               std::shared_ptr<TimeCtr>& tcIn) {
  tc = tcIn;
  fluidInterface = fluidIn;
  // read_param();
}
//---------------------------------------------------------

void Pic::read_param(const std::string& command, ReadParam& readParam) {

  if (command == "#DIVE") {
    readParam.read_var("doCorrectDivE", doCorrectDivE);
  } else if (command == "#EFIELDSOLVER") {
    Real tol;
    int nIter;
    readParam.read_var("tol", tol);
    readParam.read_var("nIter", nIter);
    eSolver.set_tol(tol);
    eSolver.set_nIter(nIter);
  } else if (command == "#PARTICLES") {
    readParam.read_var("npcelx", nPartPerCell[ix_]);
    readParam.read_var("npcely", nPartPerCell[iy_]);
    readParam.read_var("npcelz", nPartPerCell[iz_]);
  } else if (command == "#ELECTRON") {
    readParam.read_var("qom", qomEl);
  } else if (command == "#DISCRETIZE" || command == "#DISCRETIZATION") {
    readParam.read_var("theta", fsolver.theta);
    readParam.read_var("coefDiff", fsolver.coefDiff);
  } else if (command == "#RESAMPLING") {
    readParam.read_var("doReSampling", doReSampling);
    readParam.read_var("reSamplingLowLimit", reSamplingLowLimit);
    readParam.read_var("reSamplingHighLimit", reSamplingHighLimit);
  }

  fluidInterface->set_plasma_charge_and_mass(qomEl);
  nSpecies = fluidInterface->get_nS();
}

void Pic::set_ic() {
  std::string nameFunc = "Pic::set_ic";

  Print() << nameFunc << " begin" << std::endl;

  set_ic_field();
  set_ic_particles();

  Print() << nameFunc << " end" << std::endl;
}

void Pic::make_grid(int nGstIn, const BoxArray& centerBAIn,
                    const Geometry& geomIn) {
  set_nGst(nGstIn);

  centerBA = centerBAIn;

  geom = geomIn;

  nodeBA = convert(centerBA, amrex::IntVect{ AMREX_D_DECL(1, 1, 1) });

  dm.define(centerBA);
  costMF.define(centerBA, dm, 1, 0);
  costMF.setVal(0);
}

void Pic::make_data() {

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

    plasmaEnergy.resize(nSpecies + 1);

    // The last one is the sum of all species.
    nodePlasma.resize(nSpecies + 1);
    for (auto& pl : nodePlasma) {
      pl.define(nodeBA, dm, nMoments, nGst);
      pl.setVal(0.0);
    }

    // Create particle containers.
    for (int i = 0; i < nSpecies; i++) {
      auto ptr = std::make_unique<Particles>(
          geom, dm, centerBA, tc.get(), i, fluidInterface->getQiSpecies(i),
          fluidInterface->getMiSpecies(i), nPartPerCell);
      parts.push_back(std::move(ptr));
    }

    nodeMM.define(nodeBA, dm, 1, nGst);
    const RealMM mm0(0.0);
    nodeMM.setVal(mm0);

    //-------divE correction----------------
    centerNetChargeOld.define(centerBA, dm, 1, nGst);
    centerNetChargeOld.setVal(0.0);
    centerNetChargeN.define(centerBA, dm, 1, nGst);
    centerNetChargeN.setVal(0.0);
    centerNetChargeNew.define(centerBA, dm, 1, nGst);
    centerNetChargeNew.setVal(0.0);

    centerMM.define(centerBA, dm, 1, nGst);
    centerMM.setVal(0.0);
    //--------------------------------------
  }

  {
    int nGrid = get_local_node_or_cell_number(nodeE);
    eSolver.init(nGrid, nDimMax, nDimMax, matvec_E_solver);
  }

  {
    int nGrid = get_local_node_or_cell_number(centerDivE);
    divESolver.init(nGrid, 1, nDimMax, matvec_divE_accurate);
  }
}
//---------------------------------------------------------

void Pic::set_ic_field() {
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
          arrE(i, j, k, ix_) = fluidInterface->get_ex(mfi, i, j, k);
          arrE(i, j, k, iy_) = fluidInterface->get_ey(mfi, i, j, k);
          arrE(i, j, k, iz_) = fluidInterface->get_ez(mfi, i, j, k);
        }
  }

  MultiFab::Copy(nodeB, fluidInterface->get_nodeFluid(), fluidInterface->iBx, 0,
                 nodeB.nComp(), nodeB.nGrow());

  // Interpolate from node to cell center.
  average_node_to_cellcenter(centerB, 0, nodeB, 0, centerB.nComp(),
                             centerB.nGrow());

  nodeE.FillBoundary(geom.periodicity());
  nodeB.FillBoundary(geom.periodicity());
  centerB.FillBoundary(geom.periodicity());
}
//---------------------------------------------------------

void Pic::set_ic_particles() {
  for (auto& pts : parts) {
    pts->add_particles_domain(*fluidInterface);
  }

  sum_moments();
  sum_to_center(false);
}

void Pic::particle_mover() {
  std::string nameFunc = "Pic::mover";
  BL_PROFILE(nameFunc);

  timing_start(nameFunc);

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->mover(nodeEth, nodeB, tc->get_dt());

    if (doReSampling) {
      parts[i]->split_particles(reSamplingLowLimit);
      parts[i]->combine_particles(reSamplingHighLimit);
    }

    parts[i]->inject_particles_at_boundary(*fluidInterface);
  }

  timing_stop(nameFunc);
}

void Pic::sum_moments() {
  std::string nameFunc = "Pic::sum_moments";
  BL_PROFILE(nameFunc);

  timing_start(nameFunc);

  nodePlasma[nSpecies].setVal(0.0);
  const RealMM mm0(0.0);
  nodeMM.setVal(mm0);

  const Real dt = tc->get_dt();
  const auto& invDx = geom.InvCellSize();
  plasmaEnergy[iTot] = 0;
  for (int i = 0; i < nSpecies; i++) {
    PartInfo pinfo = parts[i]->sum_moments(nodePlasma[i], nodeMM, nodeB, dt);
    plasmaEnergy[i] = pinfo.energy;
    Print() << std::setprecision(5) << "Species " << i
            << ": max(uth) = " << pinfo.uMax
            << ", CFL_x = " << pinfo.uMax * dt * invDx[ix_]
            << ", CFL_y = " << pinfo.uMax * dt * invDx[iy_]
            << ", CFL_z = " << pinfo.uMax * dt * invDx[iz_] << std::endl;
    plasmaEnergy[iTot] += plasmaEnergy[i];
    MultiFab::Add(nodePlasma[nSpecies], nodePlasma[i], 0, 0, nMoments, 0);

    // Applying float boundary so that the plasma variables look right in the
    // output. It should have no influenece on the simulation results.
    apply_float_boundary(nodePlasma[i], geom, 0, nodePlasma[i].nComp(),
                         nVirGst);
  }
  Print() << "nodeMM 1 = " << nodeMM[0].array()(0, 0, 0).data[0] << std::endl;

  for (int i = 0 - 1; i <= 2 - 1; i++)
    for (int j = 0 - 1; j <= 2 - 1; j++)
      for (int k = -1; k <= 1; k++) {

        Print() << "i = " << i << " j = " << j << " k = " << k
                << "nodeMM = " << nodeMM[0].array()(i, j, k).data[0]
                << std::endl;
      }

  nodeMM.SumBoundary(geom.periodicity());

  Print() << "nodeMM 2 = " << nodeMM[0].array()(0, 0, 0).data[0] << std::endl;
  nodeMM.FillBoundary(geom.periodicity());

  Print() << "nodeMM 3 = " << nodeMM[0].array()(0, 0, 0).data[0] << std::endl;

  timing_stop(nameFunc);
}

void Pic::divE_correction() {
  std::string nameFunc = "Pic::divE_correction";
  BL_PROFILE(nameFunc);
  timing_start(nameFunc);

  for (int iIter = 0; iIter < 3; iIter++) {
    sum_to_center(true);

    Print() << "\n----------- div(E) correction at iter " << iIter
            << "----------" << std::endl;
    calculate_phi(divESolver);

    divE_correct_particle_position();
  }

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->Redistribute();
  }

  // DO CORRECTION
  sum_to_center(false);
  timing_stop(nameFunc);
}

void Pic::divE_correct_particle_position() {
  std::string nameFunc = "Pic::correct_position";
  BL_PROFILE(nameFunc);
  timing_start(nameFunc);

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->divE_correct_position(centerPhi);
  }

  timing_stop(nameFunc);
}

void Pic::calculate_phi(LinearSolver& solver) {
  std::string nameFunc = "Pic::calculate_phi";
  BL_PROFILE(nameFunc);
  timing_start(nameFunc);

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

  timing_stop(nameFunc);
}

void Pic::divE_accurate_matvec(double* vecIn, double* vecOut) {
  std::string nameFunc = "Pic::divE_matvec";
  BL_PROFILE(nameFunc);

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

void Pic::sum_to_center(bool isBeforeCorrection) {
  std::string nameFunc = "Pic::sum_to_center";
  BL_PROFILE(nameFunc);
  timing_start(nameFunc);

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
                    &Pic::get_zero);

  MultiFab::LinComb(centerNetChargeN, 1 - rhoTheta, centerNetChargeOld, 0,
                    rhoTheta, centerNetChargeNew, 0, 0,
                    centerNetChargeN.nComp(), centerNetChargeN.nGrow());

  if (!isBeforeCorrection) {
    MultiFab::Copy(centerNetChargeOld, centerNetChargeNew, 0, 0,
                   centerNetChargeOld.nComp(), centerNetChargeOld.nGrow());
  }

  timing_stop(nameFunc);
}

void Pic::update() {
  std::string nameFunc = "Pic::update";
  BL_PROFILE(nameFunc);
  timing_start(nameFunc);

  Print() << "\n================ Begin cycle = " << tc->get_cycle()
          << " at time = " << tc->get_time_si()
          << " (s) ======================" << std::endl;
  update_E();

  particle_mover();
  update_B();

  if (doCorrectDivE) {
    divE_correction();
  }

  // update time, step number.
  tc->update();

  // Apply load balance before sum_moments so that the moments and mass matrix
  // on the gird do NOT require MPI communication.
  load_balance();

  sum_moments();

  timing_stop(nameFunc);
}

void Pic::update_E() {
  std::string nameFunc = "Pic::update_E";
  BL_PROFILE(nameFunc);
  timing_start(nameFunc);

  eSolver.reset(get_local_node_or_cell_number(nodeE));

  update_E_rhs(eSolver.rhs);

  convert_3d_to_1d(nodeE, eSolver.xLeft, geom);

  update_E_matvec(eSolver.xLeft, eSolver.matvec, false);

  Real sum = 0;
  for (int i = 0; i < eSolver.get_nSolve(); i++) {
    eSolver.rhs[i] -= eSolver.matvec[i];
    eSolver.xLeft[i] = 0;
    sum += pow(eSolver.rhs[i], 2);
    Print() << "i = " << i << " rhs = " << eSolver.rhs[i]
            << " matvec = " << eSolver.matvec[i] << " sum = " << sqrt(sum)
            << std::endl;
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

  apply_external_BC(nodeE, 0, nDimMax, &Pic::get_node_E);
  apply_external_BC(nodeEth, 0, nDimMax, &Pic::get_node_E);

  // Apply float BC in order to compare with iPIC3D. It is not right to apply
  // float BC here!!!!!!!--Yuxi
  apply_float_boundary(nodeE, geom, 0, nodeE.nComp());
  apply_float_boundary(nodeEth, geom, 0, nodeEth.nComp());

  timing_stop(nameFunc);
}

void Pic::update_E_matvec(const double* vecIn, double* vecOut,
                          const bool useZeroBC) {
  std::string nameFunc = "Pic::E_matvec";
  BL_PROFILE(nameFunc);
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
    apply_external_BC(vecMF, 0, nDimMax, &Pic::get_node_E);
  }

  print_MultiFab(vecMF, "vecMF");

  lap_node_to_node(vecMF, matvecMF, dm, geom);

  Real delt2 = pow(fsolver.theta * tc->get_dt(), 2);
  matvecMF.mult(-delt2);

  print_MultiFab(matvecMF, "matvecMF1");

  { // grad(divE)
    div_node_to_center(vecMF, centerDivE, geom.InvCellSize());
    if (fsolver.coefDiff > 0) {
      // Calculate cell center E for center-to-center divE.
      average_node_to_cellcenter(tempCenter3, 0, vecMF, 0, 3,
                                 tempCenter3.nGrow());

      print_MultiFab(tempCenter3, "tempcenter3_1");

      tempCenter3.FillBoundary(geom.periodicity());

      // It seems not necessary, nor a good idea to apply float BC here. --Yuxi
      apply_float_boundary(tempCenter3, geom, 0, tempCenter3.nComp(), -1);

      print_MultiFab(tempCenter3, "tempcenter3_2");

      div_center_to_center(tempCenter3, tempCenter1, geom.InvCellSize());

      print_MultiFab(tempCenter1, "tempcenter1");

      MultiFab::LinComb(centerDivE, 1 - fsolver.coefDiff, centerDivE, 0,
                        fsolver.coefDiff, tempCenter1, 0, 0, 1,
                        centerDivE.nGrow());
    }

    centerDivE.FillBoundary(geom.periodicity());

    print_MultiFab(centerDivE, "centerDivE");

    grad_center_to_node(centerDivE, tempNode3, geom.InvCellSize());

    print_MultiFab(tempNode3, "tempnode3_1");

    tempNode3.mult(delt2);
    MultiFab::Add(matvecMF, tempNode3, 0, 0, matvecMF.nComp(), 0);
    print_MultiFab(matvecMF, "matvecMF2");
  }

  tempNode3.setVal(0);
  update_E_M_dot_E(vecMF, tempNode3);

  print_MultiFab(tempNode3, "tempnode3");

  MultiFab::Add(matvecMF, tempNode3, 0, 0, matvecMF.nComp(), 0);

  MultiFab::Add(matvecMF, vecMF, 0, 0, matvecMF.nComp(), 0);

  print_MultiFab(matvecMF, "matvecMF4");

  convert_3d_to_1d(matvecMF, vecOut, geom);
}

void Pic::update_E_M_dot_E(const MultiFab& inMF, MultiFab& outMF) {

  Print() << "nodeMM 4 = " << nodeMM[0].array()(0, 0, 0).data[0] << std::endl;
  outMF.setVal(0.0);
  Real c0 = fourPI * fsolver.theta * tc->get_dt();
  for (amrex::MFIter mfi(outMF); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const amrex::Array4<amrex::Real const>& inArr = inMF[mfi].array();
    const amrex::Array4<amrex::Real>& ourArr = outMF[mfi].array();
    const amrex::Array4<RealMM>& mmArr = nodeMM[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          Real* const data0 = mmArr(i, j, k).data;
          for (int k2 = k - 1; k2 <= k + 1; k2++)
            for (int j2 = j - 1; j2 <= j + 1; j2++)
              for (int i2 = i - 1; i2 <= i + 1; i2++) {
                const int gp = (k2 - k + 1) * 9 + (j2 - j + 1) * 3 + i2 - i + 1;
                const int idx0 = gp * 9;

                // Real* const data =
                Real* const M_I = &(data0[idx0]);
                // memccpy(M_I, &(data0[idx0]), sizeof(Real)*9);
                // for (int ii = 0; ii < 9; ii++) {
                //   M_I[ii] = data[ii];
                // }

                const double& vctX =
                    inArr(i2, j2, k2, ix_); // vectX[i2][j2][k2];
                const double& vctY = inArr(i2, j2, k2, iy_);
                const double& vctZ = inArr(i2, j2, k2, iz_);
                ourArr(i, j, k, ix_) +=
                    (vctX * M_I[0] + vctY * M_I[1] + vctZ * M_I[2]) * c0;

                if (i == 0 && j == 0 && k == 0) {
                  Print() << "----------------------------------" << std::endl;
                  Print() << "i2 = " << i2 << " j2 = " << j2 << " k2 = " << k2
                          << " vctX = " << vctX << " vctY = " << vctY
                          << " vctZ = " << vctZ << " M0 = " << M_I[0]
                          << " M1 = " << M_I[1] << " M2 = " << M_I[2]
                          << " idx0 " << idx0 << " data = "
                          << (vctX * M_I[0] + vctY * M_I[1] + vctZ * M_I[2]) *
                                 c0 << " outArr = " << ourArr(i, j, k, ix_)
                          << std::endl;

                  Print() << "----------------------------------" << std::endl;
                }

                ourArr(i, j, k, iy_) +=
                    (vctX * M_I[3] + vctY * M_I[4] + vctZ * M_I[5]) * c0;
                ourArr(i, j, k, iz_) +=
                    (vctX * M_I[6] + vctY * M_I[7] + vctZ * M_I[8]) * c0;
              }
        }
  }
}

void Pic::update_E_rhs(double* rhs) {
  MultiFab tempNode(nodeBA, dm, 3, nGst);
  tempNode.setVal(0.0);
  MultiFab temp2Node(nodeBA, dm, 3, nGst);
  temp2Node.setVal(0.0);

  apply_external_BC(centerB, 0, centerB.nComp(), &Pic::get_center_B);
  apply_external_BC(nodeB, 0, nodeB.nComp(), &Pic::get_node_B);

  const Real* invDx = geom.InvCellSize();
  curl_center_to_node(centerB, tempNode, invDx);

  MultiFab::Saxpy(temp2Node, -fourPI, nodePlasma[iTot], iJhx_, 0,
                  temp2Node.nComp(), 0);

  MultiFab::Add(temp2Node, tempNode, 0, 0, tempNode.nComp(), 0);

  temp2Node.mult(fsolver.theta * tc->get_dt());

  MultiFab::Add(temp2Node, nodeE, 0, 0, nodeE.nComp(), 0);

  convert_3d_to_1d(temp2Node, rhs, geom);
}

void Pic::update_B() {
  std::string nameFunc = "Pic::update_B";
  timing_start(nameFunc);

  MultiFab dB(centerBA, dm, 3, nGst);

  curl_node_to_center(nodeEth, dB, geom.InvCellSize());

  MultiFab::Saxpy(centerB, -tc->get_dt(), dB, 0, 0, centerB.nComp(), 0);
  centerB.FillBoundary(geom.periodicity());

  apply_external_BC(centerB, 0, centerB.nComp(), &Pic::get_center_B);

  average_center_to_node(centerB, nodeB);
  nodeB.FillBoundary(geom.periodicity());

  apply_external_BC(nodeB, 0, nodeB.nComp(), &Pic::get_node_B);

  timing_stop(nameFunc);
}

void Pic::apply_external_BC(amrex::MultiFab& mf, const int iStart,
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

Real Pic::calc_E_field_energy() {
  Real sum = 0;
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
    FArrayBox& fab = nodeE[mfi];
    const Box& box = mfi.validbox();
    const Array4<Real>& arr = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    // Do not count the right edges.
    Real sumLoc = 0;
    for (int k = lo.z; k <= hi.z - 1; ++k)
      for (int j = lo.y; j <= hi.y - 1; ++j)
        for (int i = lo.x; i <= hi.x - 1; ++i) {
          sumLoc += pow(arr(i, j, k, ix_), 2) + pow(arr(i, j, k, iy_), 2) +
                    pow(arr(i, j, k, iz_), 2);
        }

    const auto& dx = geom.CellSize();
    const Real coef = 0.5 * dx[ix_] * dx[iy_] * dx[iz_] / fourPI;
    sum += sumLoc * coef;
  }
  ParallelDescriptor::ReduceRealSum(sum,
                                    ParallelDescriptor::IOProcessorNumber());

  if (!ParallelDescriptor::IOProcessor())
    sum = 0;

  return sum;
}

Real Pic::calc_B_field_energy() {
  Real sum = 0;
  for (MFIter mfi(centerB); mfi.isValid(); ++mfi) {
    FArrayBox& fab = centerB[mfi];
    const Box& box = mfi.validbox();
    const Array4<Real>& arr = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Real sumLoc = 0;
    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          sumLoc += pow(arr(i, j, k, ix_), 2) + pow(arr(i, j, k, iy_), 2) +
                    pow(arr(i, j, k, iz_), 2);
        }

    const auto& dx = geom.CellSize();
    const Real coef = 0.5 * dx[ix_] * dx[iy_] * dx[iz_] / fourPI;
    sum += sumLoc * coef;
  }
  ParallelDescriptor::ReduceRealSum(sum,
                                    ParallelDescriptor::IOProcessorNumber());

  if (!ParallelDescriptor::IOProcessor())
    sum = 0;

  return sum;
}

void Pic::compute_cost() {
  average_node_to_cellcenter(costMF, 0, nodePlasma[iTot], iNum_, 1);
}

void Pic::load_balance() {
  if (ParallelDescriptor::NProcs() == 1)
    return;

  if (!tc->loadBalance.is_time_to())
    return;

  std::string nameFunc = "Pic::load_balance";
  timing_start(nameFunc);

  Print() << "--------- Load balancing ------------" << std::endl;

  // iDecomp++;
  Print() << "before dm = " << dm << std::endl;
  compute_cost();
  dm = DistributionMapping::makeSFC(costMF, false);
  Print() << "after dm = " << dm << std::endl;

  redistribute_FabArray(nodeE, dm);
  redistribute_FabArray(nodeEth, dm);
  redistribute_FabArray(nodeB, dm);
  redistribute_FabArray(centerB, dm);

  redistribute_FabArray(centerNetChargeOld, dm);
  redistribute_FabArray(centerNetChargeN, dm);   // false??
  redistribute_FabArray(centerNetChargeNew, dm); // false??

  redistribute_FabArray(centerDivE, dm);
  redistribute_FabArray(centerPhi, dm);

  {
    bool doMoveData = false;

    for (auto& pl : nodePlasma) {
      redistribute_FabArray(pl, dm, doMoveData);
    }

    redistribute_FabArray(nodeMM, dm, doMoveData);
    redistribute_FabArray(costMF, dm, doMoveData);
    redistribute_FabArray(centerMM, dm, doMoveData);

    redistribute_FabArray(tempNode3, dm, doMoveData);
    redistribute_FabArray(tempCenter3, dm, doMoveData);
    redistribute_FabArray(tempCenter1, dm, doMoveData);
    redistribute_FabArray(tempCenter1_1, dm, doMoveData);
  }
  // Load balance particles.
  for (int i = 0; i < nSpecies; i++) {
    parts[i]->SetParticleDistributionMap(0, dm);
    parts[i]->Redistribute();
  }

  fluidInterface->load_balance(dm);

  timing_stop(nameFunc);
}

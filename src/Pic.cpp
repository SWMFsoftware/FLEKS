#include <math.h>

#include <AMReX_MultiFabUtil.H>

#include "Pic.h"
#include "GridUtility.h"
#include "LinearSolver.h"
#include "SWMFDomains.h"
#include "Timing_c.h"
#include "Utility.h"
#include "GridInfo.h"

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

void Pic::set_geom(int nGstIn, const Geometry& geomIn) {
  set_nGst(nGstIn);

  geom = geomIn;

  // centerBA = centerBAIn;
  // nodeBA = convert(centerBA, amrex::IntVect{ AMREX_D_DECL(1, 1, 1) });
  // dm = dmIn;

  // cellStatus.define(centerBA, dm, 1, nGst);

  // costMF.define(centerBA, dm, 1, 0);
  // costMF.setVal(0);
}

void Pic::regrid(const BoxArray& centerBAIn, const DistributionMapping& dmIn) {
  std::string nameFunc = "Pic::regrid";
  Print() << nameFunc << " is runing..." << std::endl;

  if (centerBAIn == centerBA)
    return;

  centerBAOld = centerBA;
  centerBA = centerBAIn;

  nodeBA = convert(centerBA, amrex::IntVect{ AMREX_D_DECL(1, 1, 1) });
  dm = dmIn;

  //===========Move data around begin====================
  distribute_FabArray(nodeE, nodeBA, dm, 3, nGst);
  distribute_FabArray(nodeEth, nodeBA, dm, 3, nGst);
  distribute_FabArray(nodeB, nodeBA, dm, 3, nGst);
  distribute_FabArray(centerB, centerBA, dm, 3, nGst);

  distribute_FabArray(centerNetChargeOld, centerBA, dm, 1, nGst);
  distribute_FabArray(centerNetChargeN, centerBA, dm, 1, nGst);   // false??
  distribute_FabArray(centerNetChargeNew, centerBA, dm, 1, nGst); // false??

  distribute_FabArray(centerDivE, centerBA, dm, 1, nGst);
  distribute_FabArray(centerPhi, centerBA, dm, 1, nGst);

  {
    bool doMoveData = false;

    iTot = nSpecies;
    if (plasmaEnergy.empty()) {
      plasmaEnergy.resize(nSpecies + 1);
    }

    if (nodePlasma.empty()) {
      // The last one is the sum of all species.
      nodePlasma.resize(nSpecies + 1);
      for (auto& pl : nodePlasma) {
        pl.define(nodeBA, dm, nMoments, nGst);
        pl.setVal(0.0);
      }
    } else {
      for (auto& pl : nodePlasma) {
        distribute_FabArray(pl, nodeBA, dm, nMoments, nGst, doMoveData);
      }
    }

    distribute_FabArray(nodeMM, nodeBA, dm, 1, nGst, doMoveData);
    distribute_FabArray(costMF, centerBA, dm, 1, nGst, doMoveData);
    distribute_FabArray(centerMM, centerBA, dm, 1, nGst, doMoveData);

    distribute_FabArray(tempNode3, nodeBA, dm, 3, nGst, doMoveData);
    distribute_FabArray(tempCenter3, centerBA, dm, 3, nGst, doMoveData);
    distribute_FabArray(tempCenter1, centerBA, dm, 1, nGst, doMoveData);
    distribute_FabArray(tempCenter1_1, centerBA, dm, 1, nGst, doMoveData);
  }

  if (parts.empty()) {
    for (int i = 0; i < nSpecies; i++) {
      auto ptr = std::make_unique<Particles>(
          geom, dm, centerBA, tc.get(), i, fluidInterface->getQiSpecies(i),
          fluidInterface->getMiSpecies(i), nPartPerCell);
      parts.push_back(std::move(ptr));
    }
  } else {
    for (int i = 0; i < nSpecies; i++) {
      parts[i]->SetParticleBoxArray(0, centerBA);
      parts[i]->SetParticleDistributionMap(0, dm);
      parts[i]->label_particles_outside_ba();
      parts[i]->Redistribute();
    }
  }
  //===========Move data around end====================

  { //===========Label cellStatus ========================
    distribute_FabArray(cellStatus, centerBA, dm, 1, nGst, false);
    cellStatus.setVal(iBoundary_);
    cellStatus.setVal(iOnNew_, 0);
    for (MFIter mfi(cellStatus); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<int>& cellArr = cellStatus[mfi].array();
      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (centerBAOld.contains({ i, j, k })) {
              cellArr(i, j, k) = iOnOld_;
            }
          }
    }

    cellStatus.FillBoundary(geom.periodicity());

    print_MultiFab(cellStatus, "cellStatus", nGst);
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

// void Pic::make_data() {

//   {
//     // EM field
//     nodeE.define(nodeBA, dm, 3, nGst);
//     nodeE.setVal(0.0);
//     nodeEth.define(nodeBA, dm, 3, nGst);
//     nodeEth.setVal(0.0);

//     nodeB.define(nodeBA, dm, 3, nGst);
//     nodeB.setVal(0.0);

//     centerB.define(centerBA, dm, 3, nGst);
//     centerB.setVal(0.0);

//     centerDivE.define(centerBA, dm, 1, nGst);
//     centerDivE.setVal(0.0);

//     centerPhi.define(centerBA, dm, 1, nGst);
//     centerPhi.setVal(0.0);

//     tempNode3.define(nodeBA, dm, 3, nGst);
//     tempNode3.setVal(0.0);

//     tempCenter3.define(centerBA, dm, 3, nGst);
//     tempCenter3.setVal(0.0);

//     tempCenter1.define(centerBA, dm, 1, nGst);
//     tempCenter1.setVal(0.0);

//     tempCenter1_1.define(centerBA, dm, 1, nGst);
//     tempCenter1_1.setVal(0.0);
//   }

//   {
//     // Plasma
//     // nSpecies = 2;
//     iTot = nSpecies;

//     plasmaEnergy.resize(nSpecies + 1);

//     // The last one is the sum of all species.
//     nodePlasma.resize(nSpecies + 1);
//     for (auto& pl : nodePlasma) {
//       pl.define(nodeBA, dm, nMoments, nGst);
//       pl.setVal(0.0);
//     }

//     // Create particle containers.
//     for (int i = 0; i < nSpecies; i++) {
//       auto ptr = std::make_unique<Particles>(
//           geom, dm, centerBA, tc.get(), i, fluidInterface->getQiSpecies(i),
//           fluidInterface->getMiSpecies(i), nPartPerCell);
//       parts.push_back(std::move(ptr));
//     }

//     nodeMM.define(nodeBA, dm, 1, nGst);
//     const RealMM mm0(0.0);
//     nodeMM.setVal(mm0);

//     //-------divE correction----------------
//     centerNetChargeOld.define(centerBA, dm, 1, nGst);
//     centerNetChargeOld.setVal(0.0);
//     centerNetChargeN.define(centerBA, dm, 1, nGst);
//     centerNetChargeN.setVal(0.0);
//     centerNetChargeNew.define(centerBA, dm, 1, nGst);
//     centerNetChargeNew.setVal(0.0);

//     centerMM.define(centerBA, dm, 1, nGst);
//     centerMM.setVal(0.0);
//     //--------------------------------------
//   }

//   {
//     int nGrid = get_local_node_or_cell_number(nodeE);
//     eSolver.init(nGrid, nDimMax, nDimMax, matvec_E_solver);
//   }

//   {
//     int nGrid = get_local_node_or_cell_number(centerDivE);
//     divESolver.init(nGrid, 1, nDimMax, matvec_divE_accurate);
//   }
// }
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
  // Print() << "nodeMM 1 = " << nodeMM[0].array()(0, 0, 0).data[0] <<
  // std::endl;

  nodeMM.SumBoundary(geom.periodicity());

  // Print() << "nodeMM 2 = " << nodeMM[0].array()(0, 0, 0).data[0] <<
  // std::endl;
  nodeMM.FillBoundary(geom.periodicity());

  // Print() << "nodeMM 3 = " << nodeMM[0].array()(0, 0, 0).data[0] <<
  // std::endl;

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
    // The particles outside the simulation domain is marked for deletion inside
    // divE_correct_particle_position(). Redistribute() deletes these particles.
    // In order to get correct moments, re-inject particles in the ghost cells.
    parts[i]->Redistribute();
    parts[i]->inject_particles_at_boundary(*fluidInterface);
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
                    tempCenter1.nComp(), tempCenter1.nGrow());

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

  // Is this necessary?? --Yuxi
  apply_float_boundary(tempCenter1, geom, 0, tempCenter1.nComp(), -1);

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

  // Real sum = 0;
  for (int i = 0; i < eSolver.get_nSolve(); i++) {
    eSolver.rhs[i] -= eSolver.matvec[i];
    eSolver.xLeft[i] = 0;
    // sum += pow(eSolver.rhs[i], 2);
    // Print() << "i = " << i << " rhs = " << eSolver.rhs[i]
    //         << " matvec = " << eSolver.matvec[i] << " sum = " << sqrt(sum)
    //         << std::endl;
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
  // Use '-1' in order to comapre with old FLEKS.
  apply_float_boundary(nodeE, geom, 0, nodeE.nComp(), -1);
  apply_float_boundary(nodeEth, geom, 0, nodeEth.nComp(), -1);

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

  // print_MultiFab(vecMF, "vecMF1", 0);

  // The right side edges should be filled in.
  vecMF.SumBoundary(geom.periodicity());

  // print_MultiFab(vecMF, "vecMF2", 0);

  // M*E needs ghost cell information.
  vecMF.FillBoundary(geom.periodicity());

  // print_MultiFab(vecMF, "vecMF3", 0);

  if (useZeroBC) {
    // The boundary nodes would not be filled in by convert_1d_3d. So, there is
    // not need to apply zero boundary conditions again here.
  } else {
    // Even after apply_external_BC(), the outmost layer node E is still unknow.
    // See FluidInterface::calc_current for detailed explaniation.
    apply_external_BC(vecMF, 0, nDimMax, &Pic::get_node_E);
  }

  // print_MultiFab(vecMF, "vecMF5", 0);

  lap_node_to_node(vecMF, matvecMF, dm, geom);

  // print_MultiFab(matvecMF, "matvecMF1", 0);

  Real delt2 = pow(fsolver.theta * tc->get_dt(), 2);
  matvecMF.mult(-delt2);

  { // grad(divE)
    div_node_to_center(vecMF, centerDivE, geom.InvCellSize());

    if (fsolver.coefDiff > 0) {
      // Calculate cell center E for center-to-center divE.
      // The outmost boundary layer of tempCenter3 is not accurate.
      average_node_to_cellcenter(tempCenter3, 0, vecMF, 0, 3,
                                 tempCenter3.nGrow());

      // Q: Why apply float BC for all boundary ghost nodes, instead of just the
      // outmost layer?
      // A: For the example described in FluidInterface::calc_current, cell
      // (c+4, c-1) of tempCenter3-block1 is not accurate, so the values at
      // (c+4, c-2)
      // will be wrong if we only apply float BC for the outmost layer.
      apply_float_boundary(tempCenter3, geom, 0, tempCenter3.nComp());

      // print_MultiFab(tempCenter3, "tempcenter2", 2);

      div_center_to_center(tempCenter3, tempCenter1, geom.InvCellSize());
      // tempCenter1.FillBoundary(geom.periodicity());

      // print_MultiFab(tempCenter1, "tempcenter1", 1);

      MultiFab::LinComb(centerDivE, 1 - fsolver.coefDiff, centerDivE, 0,
                        fsolver.coefDiff, tempCenter1, 0, 0, 1, 1);

      // print_MultiFab(centerDivE, "centerDivE", 1);
    }

    // centerDivE.FillBoundary(geom.periodicity());

    grad_center_to_node(centerDivE, tempNode3, geom.InvCellSize());

    tempNode3.mult(delt2);
    MultiFab::Add(matvecMF, tempNode3, 0, 0, matvecMF.nComp(),
                  matvecMF.nGrow());
  }

  // print_MultiFab(matvecMF, "matvecMF2", 0);

  tempNode3.setVal(0);
  update_E_M_dot_E(vecMF, tempNode3);

  MultiFab::Add(matvecMF, tempNode3, 0, 0, matvecMF.nComp(), 0);

  // print_MultiFab(matvecMF, "matvecMF3", 0);

  MultiFab::Add(matvecMF, vecMF, 0, 0, matvecMF.nComp(), 0);

  // print_MultiFab(matvecMF, "matvecMF4", 0);

  convert_3d_to_1d(matvecMF, vecOut, geom);

  // print_MultiFab(matvecMF, "matvecMF5", 0);
}

void Pic::update_E_M_dot_E(const MultiFab& inMF, MultiFab& outMF) {
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

  // print_MultiFab(nodeB, "nodeB_2", geom, 2);
  // print_MultiFab(centerB, "centerB_2", geom, 1);

  const Real* invDx = geom.InvCellSize();
  curl_center_to_node(centerB, tempNode, invDx);

  // print_MultiFab(tempNode, "tempNode",0);

  MultiFab::Saxpy(temp2Node, -fourPI, nodePlasma[iTot], iJhx_, 0,
                  temp2Node.nComp(), temp2Node.nGrow());

  MultiFab::Add(temp2Node, tempNode, 0, 0, tempNode.nComp(), temp2Node.nGrow());

  temp2Node.mult(fsolver.theta * tc->get_dt());

  MultiFab::Add(temp2Node, nodeE, 0, 0, nodeE.nComp(), temp2Node.nGrow());

  convert_3d_to_1d(temp2Node, rhs, geom);
  // print_MultiFab(temp2Node, "temp2node", 0);
}

void Pic::update_B() {
  std::string nameFunc = "Pic::update_B";
  timing_start(nameFunc);

  MultiFab dB(centerBA, dm, 3, nGst);

  curl_node_to_center(nodeEth, dB, geom.InvCellSize());

  MultiFab::Saxpy(centerB, -tc->get_dt(), dB, 0, 0, centerB.nComp(),
                  centerB.nGrow());
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

  BoxArray ba = mf.boxArray();

  const IntVect& ngrow = mf.nGrowVect();
  for (int i = 0; i < nDimMax; ++i) {
    if (geom.isPeriodic(i)) {
      ba.grow(i, ngrow[i]);
    }
  }

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.fabbox();

    //! if there are cells not in the valid + periodic grown box
    //! we need to fill them here
    if (!ba.contains(bx)) {
      amrex::Array4<amrex::Real> const& arr = mf[mfi].array();

      const auto lo = IntVect(bx.loVect());
      const auto hi = IntVect(bx.hiVect());

      IntVect mid = (lo + hi) / 2;

      Vector<BCRec> bcr(1);
      for (int iDim = 0; iDim < nDimMax; iDim++) {
        auto idxLo = mid;
        idxLo[iDim] = lo[iDim];
        bcr[0].setLo(iDim,
                     ba.contains(idxLo) ? BCType::int_dir : BCType::ext_dir);

        auto idxHi = mid;
        idxHi[iDim] = hi[iDim];
        bcr[0].setHi(iDim,
                     ba.contains(idxHi) ? BCType::int_dir : BCType::ext_dir);
      }

      IntVect nGst = mf.nGrowVect();
      int iMin = lo[ix_], iMax = hi[ix_];
      int jMin = lo[iy_], jMax = hi[iy_];
      int kMin = lo[iz_], kMax = hi[iz_];

      // x left
      if (bcr[0].lo(ix_) == BCType::ext_dir)
        for (int iVar = iStart; iVar < nComp; iVar++)
          for (int k = kMin; k <= kMax; k++)
            for (int j = jMin; j <= jMax; j++)
              for (int i = iMin; i <= iMin + nGst[ix_] - 1 + nVirGst; i++) {
                arr(i, j, k, iVar) = (this->*func)(mfi, i, j, k, iVar - iStart);
              }

      // x right
      if (bcr[0].hi(ix_) == BCType::ext_dir)
        for (int iVar = iStart; iVar < nComp; iVar++)
          for (int k = kMin; k <= kMax; k++)
            for (int j = jMin; j <= jMax; j++)
              for (int i = iMax - nGst[ix_] + 1 - nVirGst; i <= iMax; i++) {
                arr(i, j, k, iVar) = (this->*func)(mfi, i, j, k, iVar - iStart);
              }

      // y left
      if (bcr[0].lo(iy_) == BCType::ext_dir)
        for (int iVar = iStart; iVar < nComp; iVar++)
          for (int k = kMin; k <= kMax; k++)
            for (int j = jMin; j <= jMin + nGst[iy_] - 1 + nVirGst; j++)
              for (int i = iMin; i <= iMax; i++) {
                arr(i, j, k, iVar) = (this->*func)(mfi, i, j, k, iVar - iStart);
              }

      // y right
      if (bcr[0].hi(iy_) == BCType::ext_dir)
        for (int iVar = iStart; iVar < nComp; iVar++)
          for (int k = kMin; k <= kMax; k++)
            for (int j = jMax - nGst[iy_] + 1 - nVirGst; j <= jMax; j++)
              for (int i = iMin; i <= iMax; i++) {
                arr(i, j, k, iVar) = (this->*func)(mfi, i, j, k, iVar - iStart);
              }

      // z left
      if (bcr[0].lo(iz_) == BCType::ext_dir)
        for (int iVar = iStart; iVar < nComp; iVar++)
          for (int k = kMin; k <= kMin + nGst[iz_] - 1 + nVirGst; k++)
            for (int j = jMin; j <= jMax; j++)
              for (int i = iMin; i <= iMax; i++) {
                arr(i, j, k, iVar) = (this->*func)(mfi, i, j, k, iVar - iStart);
              }

      // z right
      if (bcr[0].hi(iz_) == BCType::ext_dir)
        for (int iVar = iStart; iVar < nComp; iVar++)
          for (int k = kMax - nGst[iz_] + 1 - nVirGst; k <= kMax; k++)
            for (int j = jMin; j <= jMax; j++)
              for (int i = iMin; i <= iMax; i++) {
                arr(i, j, k, iVar) = (this->*func)(mfi, i, j, k, iVar - iStart);
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
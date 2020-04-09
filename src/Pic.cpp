#include <math.h>

#include <AMReX_MultiFabUtil.H>

#include "GridInfo.h"
#include "GridUtility.h"
#include "LinearSolver.h"
#include "Pic.h"
#include "SWMFDomains.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

//==========================================================
void Pic::init(Real timeIn, const std::string& paramString, int* paramInt,
               double* gridDim, double* paramReal,
               std::shared_ptr<FluidInterface>& fluidIn,
               std::shared_ptr<TimeCtr>& tcIn) {
  tc = tcIn;
  fluidInterface = fluidIn;
  // read_param();
}

//==========================================================
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

//==========================================================
void Pic::fill_new_cells() {
  std::string nameFunc = "Pic::fill_new_cells";

  if (!doNeedFillNewCell)
    return;

  Print() << nameFunc << " begin" << std::endl;

  timing_func(nameFunc);

  fill_E_B_fields();
  fill_particles();

  sum_moments();
  sum_to_center(false);

  doNeedFillNewCell = false;

  Print() << nameFunc << " end" << std::endl;
}

//==========================================================
void Pic::set_geom(int nGstIn, const Geometry& geomIn) {
  set_nGst(nGstIn);
  geom = geomIn;
}

//==========================================================
void Pic::regrid(const BoxArray& picRegionIn, const BoxArray& centerBAIn,
                 const DistributionMapping& dmIn) {
  std::string nameFunc = "Pic::regrid";
  Print() << nameFunc << " is runing..." << std::endl;

  timing_func(nameFunc);

  if (centerBAIn == centerBA)
    return;

  doNeedFillNewCell = true;

  picRegionBA = picRegionIn;

  centerBAOld = centerBA;
  nodeBAOld = convert(centerBAOld, amrex::IntVect{ AMREX_D_DECL(1, 1, 1) });

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

    distribute_FabArray(nodeMM, nodeBA, dm, 1, 1, doMoveData);
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
      // TODO: Combine the following function with the constructor.
      ptr->set_region_ba(picRegionBA);
      parts.push_back(std::move(ptr));
    }
  } else {
    for (int i = 0; i < nSpecies; i++) {
      // Label the particles outside the OLD PIC region.
      parts[i]->label_particles_outside_ba();
      parts[i]->SetParticleBoxArray(0, centerBA);
      parts[i]->set_region_ba(picRegionBA);
      parts[i]->SetParticleDistributionMap(0, dm);
      // Label the particles outside the NEW PIC region.
      parts[i]->label_particles_outside_ba();
      parts[i]->Redistribute();
    }
  }
  //===========Move data around end====================

  { //===========Label cellStatus/nodeStatus ==========

    // The algorithm decides inject particles or not needs at least 2 ghost cell
    // layers.
    distribute_FabArray(cellStatus, centerBA, dm, 1, nGst >= 2 ? nGst : 2,
                        false);
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
            if (cellArr(i, j, k) == iOnNew_ &&
                centerBAOld.contains({ i, j, k })) {
              cellArr(i, j, k) = iOnOld_;
            }
          }
    }

    cellStatus.FillBoundary(geom.periodicity());

    // print_MultiFab(cellStatus, "cellStatus", 1);

    distribute_FabArray(nodeStatus, nodeBA, dm, 1, nGst, false);
    nodeStatus.setVal(iBoundary_);
    nodeStatus.setVal(iOnNew_, 0);

    for (MFIter mfi(nodeStatus); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto& nodeArr = nodeStatus[mfi].array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (nodeArr(i, j, k) == iOnNew_ &&
                nodeBAOld.contains({ i, j, k })) {
              nodeArr(i, j, k) = iOnOld_;
            }
          }
    }

    nodeStatus.FillBoundary(geom.periodicity());

    distribute_FabArray(nodeType, nodeBA, dm, 1, 0, false);
    set_nodeType();

    // print_MultiFab(nodeType, "nodeType");
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

//==========================================================
void Pic::set_nodeType() {
  nodeType.setVal(iNotHandle_);
  const Box& gbx = convert(geom.Domain(), { 0, 0, 0 });

  for (MFIter mfi(nodeType); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();

    const Box& cellBox = convert(box, { 0, 0, 0 });

    const auto& typeArr = nodeType[mfi].array();
    const auto& statusArr = cellStatus[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    bool is2D = false;
    if (geom.isPeriodic(iz_) && (gbx.bigEnd(iz_) == gbx.smallEnd(iz_)))
      is2D = true;

    int iMin = lo.x + 1, iMax = hi.x - 1;
    int jMin = lo.y + 1, jMax = hi.y - 1;
    int kMin = lo.z + 1, kMax = hi.z - 1;
    if (is2D) {
      kMin = lo.z;
      kMax = lo.z;
    }

    int diMax = 0, diMin = -1;
    int djMax = 0, djMin = -1;
    int dkMax = 0, dkMin = -1;
    if (is2D) {
      dkMin = 0;
    }

    auto fHandle = [&](int i, int j, int k) {
      for (int dk = dkMax; dk >= dkMin; dk--)
        for (int dj = djMax; dj >= djMin; dj--)
          for (int di = diMax; di >= diMin; di--) {
            if (statusArr(i + di, j + dj, k + dk) != iBoundary_) {
              // Find the first CELL that shares this node.
              if (cellBox.contains({ i + di, j + dj, k + dk })) {
                return true;
              } else {
                return false;
              }
            }
          }
      Abort("Error: something is wrong here!");
      return false;
    };

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          if (!is2D || k == lo.z) {
            // for 2D (1 cell in the z-direction), only handle the layer of
            // k=0
            if (i == lo.x || i == hi.x || j == lo.y || j == hi.y ||
                (!is2D && (k == lo.z || k == hi.z))) {
              // Block boundary nodes.

              // Check if this block boundary node needs to be handled by this
              // block.

              if (fHandle(i, j, k)) {
                typeArr(i, j, k) = iHandle_;
              }
            } else {
              typeArr(i, j, k) = iHandle_;
            }
          }
        }
  }
}

//==========================================================
void Pic::fill_new_node_E() {

  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
    FArrayBox& fab = nodeE[mfi];
    const Box& box = mfi.validbox();
    const Array4<Real>& arrE = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    const auto& status = nodeStatus[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          if (status(i, j, k) == iOnNew_) {
            arrE(i, j, k, ix_) = fluidInterface->get_ex(mfi, i, j, k);
            arrE(i, j, k, iy_) = fluidInterface->get_ey(mfi, i, j, k);
            arrE(i, j, k, iz_) = fluidInterface->get_ez(mfi, i, j, k);
          }
        }
  }
}

//==========================================================
void Pic::fill_new_node_B() {
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrB = nodeB[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    const auto& status = nodeStatus[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          if (status(i, j, k) == iOnNew_) {
            arrB(i, j, k, ix_) = fluidInterface->get_bx(mfi, i, j, k);
            arrB(i, j, k, iy_) = fluidInterface->get_by(mfi, i, j, k);
            arrB(i, j, k, iz_) = fluidInterface->get_bz(mfi, i, j, k);
          }
        }
  }
}

//==========================================================
void Pic::fill_new_center_B() {
  for (MFIter mfi(centerB); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& centerArr = centerB[mfi].array();
    const auto& nodeArr = nodeB[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    const auto& status = cellStatus[mfi].array();

    for (int iVar = 0; iVar < centerB.nComp(); iVar++)
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (status(i, j, k) != iOnOld_) {
              centerArr(i, j, k, iVar) = 0;
              for (int di = 0; di <= 1; di++)
                for (int dj = 0; dj <= 1; dj++)
                  for (int dk = 0; dk <= 1; dk++) {
                    centerArr(i, j, k, iVar) +=
                        0.125 * nodeArr(i + di, j + dj, k + dk, iVar);
                  }
            }
          }
  }
}

//==========================================================
void Pic::fill_E_B_fields() {

  fill_new_node_E();
  fill_new_node_B();

  nodeE.FillBoundary(geom.periodicity());
  nodeB.FillBoundary(geom.periodicity());

  apply_external_BC(nodeStatus, nodeE, 0, nDimMax, &Pic::get_node_E);
  apply_external_BC(nodeStatus, nodeB, 0, nDimMax, &Pic::get_node_B);

  fill_new_center_B();
  centerB.FillBoundary(geom.periodicity());

  apply_external_BC(cellStatus, centerB, 0, centerB.nComp(),
                    &Pic::get_center_B);
}

//==========================================================
void Pic::fill_particles() {
  inject_particles_for_new_cells();
  inject_particles_for_boundary_cells();
}

//==========================================================
void Pic::particle_mover() {
  std::string nameFunc = "Pic::mover";

  timing_func(nameFunc);

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->mover(nodeEth, nodeB, tc->get_dt());

    if (doReSampling) {
      parts[i]->split_particles(reSamplingLowLimit);
      parts[i]->combine_particles(reSamplingHighLimit);
    }
  }

  inject_particles_for_boundary_cells();
}

//==========================================================
void Pic::sum_moments() {
  std::string nameFunc = "Pic::sum_moments";

  timing_func(nameFunc);

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
  }

  nodeMM.SumBoundary(geom.periodicity());

  nodeMM.FillBoundary(geom.periodicity());
}

//==========================================================
void Pic::divE_correction() {
  std::string nameFunc = "Pic::divE_correction";

  timing_func(nameFunc);

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
  }
  inject_particles_for_boundary_cells();

  sum_to_center(false);
}

//==========================================================
void Pic::divE_correct_particle_position() {
  std::string nameFunc = "Pic::correct_position";

  timing_func(nameFunc);

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->divE_correct_position(centerPhi);
  }
}

//==========================================================
void Pic::calculate_phi(LinearSolver& solver) {
  std::string nameFunc = "Pic::calculate_phi";

  timing_func(nameFunc);

  solver.reset(get_local_node_or_cell_number(centerDivE));

  div_node_to_center(nodeE, tempCenter1, geom.InvCellSize());

  MultiFab::LinComb(tempCenter1, 1.0 / rhoTheta, tempCenter1, 0,
                    -fourPI / rhoTheta, centerNetChargeN, 0, 0,
                    tempCenter1.nComp(), tempCenter1.nGrow());

  convert_3d_to_1d(tempCenter1, solver.rhs, geom);

  solver.solve();

  convert_1d_to_3d(solver.xLeft, centerPhi, geom);
  centerPhi.FillBoundary(geom.periodicity());
}

//==========================================================
void Pic::divE_accurate_matvec(double* vecIn, double* vecOut) {
  std::string nameFunc = "Pic::divE_matvec";

  zero_array(vecOut, divESolver.get_nSolve());

  convert_1d_to_3d(vecIn, tempCenter1, geom);
  tempCenter1.FillBoundary(0, 1, IntVect(1), geom.periodicity());  

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

//==========================================================
void Pic::sum_to_center(bool isBeforeCorrection) {
  std::string nameFunc = "Pic::sum_to_center";

  timing_func(nameFunc);

  centerNetChargeNew.setVal(0.0);

  const RealCMM mm0(0.0);
  centerMM.setVal(mm0);

  bool doNetChargeOnly = !isBeforeCorrection;

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->sum_to_center(centerNetChargeNew, centerMM, doNetChargeOnly);
  }

  if (!doNetChargeOnly) {
    centerMM.SumBoundary(geom.periodicity());
  }

  centerNetChargeNew.SumBoundary(geom.periodicity());

  apply_external_BC(cellStatus, centerNetChargeNew, 0,
                    centerNetChargeNew.nComp(), &Pic::get_zero);

  MultiFab::LinComb(centerNetChargeN, 1 - rhoTheta, centerNetChargeOld, 0,
                    rhoTheta, centerNetChargeNew, 0, 0,
                    centerNetChargeN.nComp(), centerNetChargeN.nGrow());

  if (!isBeforeCorrection) {
    MultiFab::Copy(centerNetChargeOld, centerNetChargeNew, 0, 0,
                   centerNetChargeOld.nComp(), centerNetChargeOld.nGrow());
  }
}

//==========================================================
void Pic::update() {
  std::string nameFunc = "Pic::update";

  timing_func(nameFunc);

  {
    const Real t0 = tc->get_time_si();
    // update time, step number.
    tc->update();
    const Real t1 = tc->get_time_si();
    Print() << "\n================ Begin cycle = " << tc->get_cycle()
            << " from t = " << t0 << " (s) to t = " << t1
            << " (s) ======================" << std::endl;
  }

  update_E();

  particle_mover();
  update_B();

  if (doCorrectDivE) {
    divE_correction();
  }

  // Apply load balance before sum_moments so that the moments and mass matrix
  // on the gird do NOT require MPI communication.
  load_balance();

  sum_moments();
}

//==========================================================
void Pic::update_E() {
  std::string nameFunc = "Pic::update_E";

  timing_func(nameFunc);

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

  apply_external_BC(nodeStatus, nodeE, 0, nDimMax, &Pic::get_node_E);
  apply_external_BC(nodeStatus, nodeEth, 0, nDimMax, &Pic::get_node_E);

}

//==========================================================
void Pic::update_E_matvec(const double* vecIn, double* vecOut,
                          const bool useZeroBC) {
  std::string nameFunc = "Pic::E_matvec";
  timing_func(nameFunc);  

  zero_array(vecOut, eSolver.get_nSolve());

  MultiFab vecMF(nodeBA, dm, 3, nGst);
  vecMF.setVal(0.0);

  MultiFab matvecMF(nodeBA, dm, 3, 1);
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
    // Even after apply_external_BC(), the outmost layer node E is still unknow.
    // See FluidInterface::calc_current for detailed explaniation.
    apply_external_BC(nodeStatus, vecMF, 0, nDimMax, &Pic::get_node_E);
  }

  lap_node_to_node(vecMF, matvecMF, dm, geom, cellStatus);

  Real delt2 = pow(fsolver.theta * tc->get_dt(), 2);
  matvecMF.mult(-delt2);

  { // grad(divE)
    div_node_to_center(vecMF, centerDivE, geom.InvCellSize());

    if (fsolver.coefDiff > 0) {
      // Calculate cell center E for center-to-center divE.
      // The outmost boundary layer of tempCenter3 is not accurate.
      average_node_to_cellcenter(tempCenter3, 0, vecMF, 0, 3,
                                 tempCenter3.nGrow());

      //----The following comments are left here for reference------
      // Q: Why apply float BC for all boundary ghost nodes, instead of just the
      // outmost layer?
      // A: For the example described in FluidInterface::calc_current, cell
      // (c+4, c-1) of tempCenter3-block1 is not accurate, so the values at
      // (c+4, c-2)
      // will be wrong if we only apply float BC for the outmost layer.
      // apply_float_boundary(cellStatus, tempCenter3, geom, 0,
      //                           tempCenter3.nComp());
      //------------------------------------------------------------

      div_center_to_center(tempCenter3, tempCenter1, geom.InvCellSize());

      tempCenter1.FillBoundary(0, 1, IntVect(1), geom.periodicity());

      // 1) The outmost boundary layer of tempCenter3 is not accurate.
      // 2) The 2 outmost boundary layers (all ghosts if there are 2 ghost
      // cells) of tempCenter1 are not accurate
      apply_external_BC(cellStatus, tempCenter1, 0, tempCenter1.nComp(),
                        &Pic::get_zero);

      // print_MultiFab(tempCenter1, "tempCenter1", 1);

      MultiFab::LinComb(centerDivE, 1 - fsolver.coefDiff, centerDivE, 0,
                        fsolver.coefDiff, tempCenter1, 0, 0, 1, 1);
    }

    // print_MultiFab(centerDivE, "centerDivE", 1);
    grad_center_to_node(centerDivE, tempNode3, geom.InvCellSize());
    // print_MultiFab(tempNode3, "tempnode3", geom, 0);

    tempNode3.mult(delt2);
    MultiFab::Add(matvecMF, tempNode3, 0, 0, matvecMF.nComp(),
                  matvecMF.nGrow());
  }

  tempNode3.setVal(0);
  update_E_M_dot_E(vecMF, tempNode3);

  MultiFab::Add(matvecMF, tempNode3, 0, 0, matvecMF.nComp(), 0);

  MultiFab::Add(matvecMF, vecMF, 0, 0, matvecMF.nComp(), 0);

  convert_3d_to_1d(matvecMF, vecOut, geom);
  
}

//==========================================================
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

                Real* const M_I = &(data0[idx0]);

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

//==========================================================
void Pic::update_E_rhs(double* rhs) {
  MultiFab tempNode(nodeBA, dm, 3, nGst);
  tempNode.setVal(0.0);
  MultiFab temp2Node(nodeBA, dm, 3, nGst);
  temp2Node.setVal(0.0);

  apply_external_BC(cellStatus, centerB, 0, centerB.nComp(),
                    &Pic::get_center_B);
  apply_external_BC(nodeStatus, nodeB, 0, nodeB.nComp(), &Pic::get_node_B);

  const Real* invDx = geom.InvCellSize();
  curl_center_to_node(centerB, tempNode, invDx);

  MultiFab::Saxpy(temp2Node, -fourPI, nodePlasma[iTot], iJhx_, 0,
                  temp2Node.nComp(), temp2Node.nGrow());

  MultiFab::Add(temp2Node, tempNode, 0, 0, tempNode.nComp(), temp2Node.nGrow());

  temp2Node.mult(fsolver.theta * tc->get_dt());

  MultiFab::Add(temp2Node, nodeE, 0, 0, nodeE.nComp(), temp2Node.nGrow());

  convert_3d_to_1d(temp2Node, rhs, geom);
}

//==========================================================
void Pic::update_B() {
  std::string nameFunc = "Pic::update_B";
  timing_func(nameFunc);

  MultiFab dB(centerBA, dm, 3, nGst);

  curl_node_to_center(nodeEth, dB, geom.InvCellSize());

  MultiFab::Saxpy(centerB, -tc->get_dt(), dB, 0, 0, centerB.nComp(),
                  centerB.nGrow());
  centerB.FillBoundary(geom.periodicity());

  apply_external_BC(cellStatus, centerB, 0, centerB.nComp(),
                    &Pic::get_center_B);

  average_center_to_node(centerB, nodeB);
  nodeB.FillBoundary(geom.periodicity());

  apply_external_BC(nodeStatus, nodeB, 0, nodeB.nComp(), &Pic::get_node_B);
}

//==========================================================
void Pic::apply_external_BC(const iMultiFab& status, MultiFab& mf,
                            const int iStart, const int nComp, GETVALUE func) {

  if (geom.isAllPeriodic())
    return;
  if (mf.nGrow() == 0)
    return;

  BoxArray ba = convert(picRegionBA, mf.boxArray().ixType()); 

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

      const Array4<const int>& statusArr = status[mfi].array();

      const auto lo = IntVect(bx.loVect());
      const auto hi = IntVect(bx.hiVect());

      for (int iVar = iStart; iVar < nComp; iVar++)
        for (int k = lo[iz_] + 1; k <= hi[iz_] - 1; k++)
          for (int j = lo[iy_]; j <= hi[iy_]; j++)
            for (int i = lo[ix_]; i <= hi[ix_]; i++)
              if (statusArr(i, j, k, 0) == iBoundary_) {
                arr(i, j, k, iVar) = (this->*func)(mfi, i, j, k, iVar - iStart);
              }

      continue;
    }
  }
}

//==========================================================
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

//==========================================================
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

//==========================================================
void Pic::compute_cost() {
  average_node_to_cellcenter(costMF, 0, nodePlasma[iTot], iNum_, 1);
}

//==========================================================
void Pic::load_balance() {
  if (ParallelDescriptor::NProcs() == 1)
    return;

  if (!tc->loadBalance.is_time_to())
    return;

  std::string nameFunc = "Pic::load_balance";
  timing_func(nameFunc);

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
}

//==========================================================
void Pic::convert_1d_to_3d(const double* const p, amrex::MultiFab& MF,
                           amrex::Geometry& geom) {
  bool isCenter = MF.ixType().cellCentered();
  bool isNode = !isCenter;

  MF.setVal(0.0);

  const Box& gbx = convert(geom.Domain(), MF.boxArray().ixType());

  int iCount = 0;
  for (amrex::MFIter mfi(MF, doTiling); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.tilebox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const amrex::Array4<amrex::Real>& arr = MF[mfi].array();

    int iMax = hi.x, jMax = hi.y, kMax = hi.z;
    int iMin = lo.x, jMin = lo.y, kMin = lo.z;

    const auto& nodeArr = nodeType[mfi].array();
    for (int iVar = 0; iVar < MF.nComp(); iVar++)
      for (int k = kMin; k <= kMax; ++k)
        for (int j = jMin; j <= jMax; ++j)
          for (int i = iMin; i <= iMax; ++i)
            if (isCenter || nodeArr(i, j, k) == iHandle_) {
              arr(i, j, k, iVar) = p[iCount++];
            }
  }
}

//==========================================================
void Pic::convert_3d_to_1d(const amrex::MultiFab& MF, double* const p,
                           amrex::Geometry& geom) {

  bool isCenter = MF.ixType().cellCentered();
  bool isNode = !isCenter;

  const Box& gbx = convert(geom.Domain(), MF.boxArray().ixType());

  int iCount = 0;
  for (amrex::MFIter mfi(MF, doTiling); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.tilebox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const amrex::Array4<amrex::Real const>& arr = MF[mfi].array();

    // Avoid double counting the share edges.
    int iMax = hi.x, jMax = hi.y, kMax = hi.z;
    int iMin = lo.x, jMin = lo.y, kMin = lo.z;

    const auto& nodeArr = nodeType[mfi].array();
    for (int iVar = 0; iVar < MF.nComp(); iVar++)
      for (int k = kMin; k <= kMax; ++k)
        for (int j = jMin; j <= jMax; ++j)
          for (int i = iMin; i <= iMax; ++i)
            if (isCenter || nodeArr(i, j, k) == iHandle_) {
              p[iCount++] = arr(i, j, k, iVar);
            }
  }
}

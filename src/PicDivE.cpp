#include <cmath>
#include <vector>

#include <AMReX_MultiFabUtil.H>

#include "GridUtility.h"
#include "LinearSolver.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

//==========================================================
void Pic::divE_correction() {
  std::string nameFunc = "Pic::divE_correction";

  timing_func(nameFunc);

  for (int iIter = 0; iIter < nDivECorrection; iIter++) {

    sum_to_center(true);

    if (doReport)
      Print() << "\n-----" << printPrefix << " div(E) correction at iter "
              << iIter << "----------" << std::endl;

    calculate_phi(divESolver, 0);

    divE_correct_particle_position();
  }

  for (int i = 0; i < nSpecies; ++i) {
    // The particles outside the simulation domain is marked for deletion
    // inside divE_correct_particle_position(). redistribute_particles()
    // deletes these particles. In order to get correct moments, re-inject
    // particles in the ghost cells.
    parts[i]->redistribute_particles();
  }

  inject_particles_for_boundary_cells();

  sum_to_center(false);
}

//==========================================================
void Pic::divE_correct_particle_position() {
  std::string nameFunc = "Pic::correct_position";

  timing_func(nameFunc);
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    for (int i = 0; i < nSpecies; ++i) {
      parts[i]->divE_correct_position(centerPhi, iLev);
    }
  }
}

//==========================================================
void Pic::calculate_phi(LinearSolver& solver, int iLev) {
  std::string nameFunc = "Pic::calculate_phi";

  timing_func(nameFunc);

  {
    MultiFab residual(cGrids[iLev], DistributionMap(iLev), 1, nGst);

    solver.reset(get_local_node_or_cell_number(centerDivE[iLev]));
    // div_node_to_center(nodeE[iLev], residual, Geom(iLev).InvCellSize());
    MultiFab::Copy(residual, centerDivE[iLev], 0, 0, 1, nGst);
    Real coef = 1.0 / rhoTheta;

    MultiFab::LinComb(residual, coef, residual, 0, -fourPI * coef,
                      centerNetChargeN[iLev], 0, 0, residual.nComp(),
                      residual.nGrow());
    if (finest_level > 0) {
      skip_cells_divE_correction(residual, cellStatus[iLev], iLev);
    }

    convert_3d_to_1d(residual, solver.rhs, iLev);

    BL_PROFILE_VAR("Pic::phi_iterate", solve);
    solver.solve(iLev, doReport);
    BL_PROFILE_VAR_STOP(solve);

    convert_1d_to_3d(solver.xLeft, centerPhi[iLev], iLev);
    centerPhi[iLev].FillBoundary(Geom(iLev).periodicity());
  }
}

//==========================================================
void Pic::divE_accurate_matvec(const double* vecIn, double* vecOut, int iLev) {
  std::string nameFunc = "Pic::divE_matvec";
  timing_func(nameFunc);

  // const int iLev = 0;
  zero_array(vecOut, divESolver.get_nSolve());

  MultiFab inMF(cGrids[iLev], DistributionMap(iLev), 1, nGst);

  convert_1d_to_3d(vecIn, inMF, iLev);
  inMF.FillBoundary(0, 1, IntVect(1), Geom(iLev).periodicity());

  MultiFab outMF(cGrids[iLev], DistributionMap(iLev), 1, nGst);
  outMF.setVal(0.0);

  for (MFIter mfi(inMF); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();

    const Array4<Real>& lArr = outMF[mfi].array();
    const Array4<Real const>& rArr = inMF[mfi].array();
    const Array4<RealCMM>& mmArr = centerMM[iLev][mfi].array();

    ParallelFor(box, [&](int i, int j, int k) {
      IntVect ijk = { AMREX_D_DECL(i, j, k) };
      Box subBox(ijk - 1, ijk + 1);

      ParallelFor(subBox, [&](int i2, int j2, int k2) {
        const int gp = (i2 - i + 1) * 9 + (j2 - j + 1) * 3 + k2 - k + 1;
        lArr(i, j, k) += rArr(i2, j2, k2) * mmArr(i, j, k)[gp];
      });
    });
  }
  outMF.mult(fourPI * fourPI);
  convert_3d_to_1d(outMF, vecOut, iLev);
}

//==========================================================
void Pic::sum_to_center(bool isBeforeCorrection) {
  std::string nameFunc = "Pic::sum_to_center";

  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    centerNetChargeNew[iLev].setVal(0.0);

    const RealCMM mm0(0.0);
    centerMM[iLev].setVal(mm0);

    bool doNetChargeOnly = !isBeforeCorrection;

    for (int i = 0; i < nSpecies; ++i) {
      parts[i]->sum_to_center(centerNetChargeNew[iLev], centerMM[iLev],
                              doNetChargeOnly, iLev);
    }
    if (!doNetChargeOnly) {
      centerMM[iLev].SumBoundary(Geom(iLev).periodicity());
    }

    centerNetChargeNew[iLev].SumBoundary(Geom(iLev).periodicity());

    if (iLev == 0) {
      apply_BC(cellStatus[iLev], centerNetChargeNew[iLev], 0,
               centerNetChargeNew[iLev].nComp(), &Pic::get_zero, iLev);
    }

    MultiFab::LinComb(
        centerNetChargeN[iLev], 1 - rhoTheta, centerNetChargeOld[iLev], 0,
        rhoTheta, centerNetChargeNew[iLev], 0, 0,
        centerNetChargeN[iLev].nComp(), centerNetChargeN[iLev].nGrow());

    if (!isBeforeCorrection) {
      MultiFab::Copy(centerNetChargeOld[iLev], centerNetChargeNew[iLev], 0, 0,
                     centerNetChargeOld[iLev].nComp(),
                     centerNetChargeOld[iLev].nGrow());
    }
  }
}

//==========================================================
void Pic::sum_to_center_amr(bool isBeforeCorrection, int iLev) {

  std::string nameFunc = "Pic::sum_to_center";
  timing_func(nameFunc);

  bool doNetChargeOnly = !isBeforeCorrection;

  centerNetChargeNew[iLev].setVal(0.0);
  const RealCMM mm0(0.0);
  centerMM[iLev].setVal(mm0);

  MultiFab jf;
  MultiFab jc;
  int fLev = iLev + 1;
  int cLev = iLev - 1;
  if (iLev == 0) {
    cLev = iLev;
  }
  if (iLev == finest_level) {
    fLev = iLev;
  }
  {
    BoxArray bac = centerB[cLev].boxArray();
    bac.refine(ref_ratio[iLev]);
    jc.define(bac, centerB[cLev].DistributionMap(), 1, 1);
    jc.setVal(0.0);
    BoxArray baf = centerB[fLev].boxArray();
    baf.coarsen(ref_ratio[iLev]);
    baf.grow(1);
    jf.define(baf, centerB[fLev].DistributionMap(), 1, 1);
    jf.setVal(0.0);
  }
  for (int i = 0; i < nSpecies; ++i) {
    parts[i]->sum_to_center_amr(centerNetChargeNew[iLev], jc, jf,
                                centerMM[iLev], doNetChargeOnly, iLev);
  }

  if (!doNetChargeOnly) {
    centerMM[iLev].SumBoundary(Geom(iLev).periodicity());
  }

  centerNetChargeNew[iLev].SumBoundary(Geom(iLev).periodicity());
  jc.SumBoundary();
  jf.SumBoundary();
  amrex::MultiFab tmp;
  tmp.define(centerB[iLev].boxArray(), centerB[iLev].DistributionMap(), 1, 0);
  tmp.setVal(0.0);
  tmp.ParallelCopy(jc);
  MultiFab::Add(centerNetChargeNew[iLev], tmp, 0, 0, 1, 0);
  tmp.setVal(0.0);
  tmp.ParallelCopy(jf);
  MultiFab::Add(centerNetChargeNew[iLev], tmp, 0, 0, 1, 0);

  if (iLev == 0) {
    apply_BC(cellStatus[iLev], centerNetChargeNew[iLev], 0,
             centerNetChargeNew[iLev].nComp(), &Pic::get_zero, iLev);
  }

  MultiFab::LinComb(
      centerNetChargeN[iLev], 1 - rhoTheta, centerNetChargeOld[iLev], 0,
      rhoTheta, centerNetChargeNew[iLev], 0, 0, centerNetChargeN[iLev].nComp(),
      centerNetChargeN[iLev].nGrow());

  if (!isBeforeCorrection) {
    MultiFab::Copy(centerNetChargeOld[iLev], centerNetChargeNew[iLev], 0, 0,
                   centerNetChargeOld[iLev].nComp(),
                   centerNetChargeOld[iLev].nGrow());
  }
}

//==========================================================

//==========================================================
void Pic::amr_divE_correction() {
  std::string nameFunc = "Pic::divE_correction";

  timing_func(nameFunc);

  for (int iIter = 0; iIter < nDivECorrection; iIter++) {
    for (int iLev = finest_level; iLev >= 0; iLev--) {
      sum_to_center_amr(true, iLev);
      skip_cells_divE_correction(centerMM[iLev], cell_status(iLev), iLev);
      calculate_phi(divESolver, iLev);
      for (int i = 0; i < nSpecies; ++i) {
        parts[i]->divE_correct_position(centerPhi, iLev);
      }
      if (finest_level > 0) {
        for (int i = 0; i < nSpecies; ++i) {
          parts[i]->Redistribute();
        }
      }
    }
  }

  inject_particles_for_boundary_cells();
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    sum_to_center_amr(false, iLev);
    if (iLev > 0) {
      skip_cells_divE_correction(centerNetChargeN[iLev], cell_status(iLev),
                                 iLev);
      skip_cells_divE_correction(centerDivE[iLev], cell_status(iLev), iLev);
    }
  }
}

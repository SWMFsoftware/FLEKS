// Composite-grid ECSIM field solve for AMR (see doc/CompositeECSIM.tex and
// Lapenta, JCP 334 (2017) 349).
//
// One implicit solve for E^{n+theta} on all levels. The unknowns are the
// corner nodes of the leaf cells (cells not covered by a finer level), one
// value per physical location:
//  - the nodes of the coarser level own the nodes on a coarse-fine interface;
//  - the fine nodes on the interface (level edge nodes) are not unknowns:
//    they are interpolated along the interface from the coarse nodes
//    (prolongation P; the hanging nodes get the average of their neighbours);
//  - the coarse nodes covered by the finer level are not unknowns either.
// The two conditions of the ECSIM energy proof then hold on the whole grid:
//  1. The current in Ampere's law is deposited with the transpose of the
//     gather: every particle deposits jHat and the mass matrix on the corner
//     nodes of its own cell only, and the fine edge nodes are folded onto
//     the coarse nodes with P^T.
//  2. The curl of B in Ampere's law is the transpose of the curl of E in
//     Faraday's law, both restricted to the leaf cells.
// The conserved energy is sum_c V_c B_c^2 / 8 pi over the leaf cells, plus
// sum_g V_g E_g^2 / 8 pi over the unknown nodes, plus the kinetic energy. The
// node volume V_g is the 1/8 of the volumes of the leaf cells around the node,
// with the volumes of the fine edge nodes folded onto the coarse nodes.
//
// The loops of this file run on the host: the solver is not available on
// GPUs yet (see composite_setup).

#include <AMReX_MultiFabUtil.H>

#include "GridUtility.h"
#include "LinearSolver.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

namespace {
// Coarse parents of the fine node index i for the refinement ratio r: the
// lower parent ic and the weight w1 of the upper parent ic + 1.
AMREX_FORCE_INLINE void coarse_parent(int i, int r, int& ic, Real& w1) {
  ic = amrex::coarsen(i, r);
  w1 = static_cast<Real>(i - ic * r) / r;
}

// The cells around the node (i, j, k) are (i-1..i, j-1..j, k+kLoCell..k):
// 8 cells in 3D, and 4 cells with a 2D AMReX library, where k is always 0.
constexpr int kLoCell = nDim > 2 ? -1 : 0;
constexpr int nCellPerNode = nDim > 2 ? 8 : 4;
} // namespace

//==========================================================
// Fine edge nodes of every level := interpolation from the next coarser level
// (coarsest level first, so that the prolongations compose).
void Pic::composite_prolong(Vector<MultiFab>& mf, int nComp) {
  for (int iLev = 1; iLev < n_lev(); iLev++) {
    const IntVect ratio = refRatio(iLev - 1);
    BoxArray bac = mf[iLev].boxArray();
    bac.coarsen(ratio);
    MultiFab ctmp(bac, mf[iLev].DistributionMap(), nComp, 0);
    ctmp.setVal(0.0);
    ctmp.ParallelCopy(mf[iLev - 1], 0, 0, nComp, IntVect(0), IntVect(0),
                      Geom(iLev - 1).periodicity());

    for (MFIter mfi(mf[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& fine = mf[iLev][mfi].array();
      const Array4<Real const>& crse = ctmp[mfi].array();
      const auto& status = node_status(iLev)[mfi].array();
      amrex::Loop(box, [&](int i, int j, int k) {
        if (!bit::is_lev_edge(status(i, j, k)))
          return;
        int ic, jc, kc;
        Real wx, wy, wz;
        coarse_parent(i, ratio[ix_], ic, wx);
        coarse_parent(j, ratio[iy_], jc, wy);
        coarse_parent(k, nDim > 2 ? ratio[iz_] : 1, kc, wz);
        for (int iVar = 0; iVar < nComp; iVar++) {
          Real val = 0;
          for (int dk = 0; dk <= 1; dk++) {
            const Real wk = dk == 0 ? 1 - wz : wz;
            if (wk == 0)
              continue;
            for (int dj = 0; dj <= 1; dj++) {
              const Real wj = dj == 0 ? 1 - wy : wy;
              if (wj == 0)
                continue;
              for (int di = 0; di <= 1; di++) {
                const Real wi = di == 0 ? 1 - wx : wx;
                if (wi == 0)
                  continue;
                val += wi * wj * wk * crse(ic + di, jc + dj, kc + dk, iVar);
              }
            }
          }
          fine(i, j, k, iVar) = val;
        }
      });
    }
    mf[iLev].FillBoundary(Geom(iLev).periodicity());
  }
}

//==========================================================
// The transpose of composite_prolong: the values on the fine edge nodes are
// added to the coarse nodes with the interpolation weights (finest level
// first). Every physical node is taken once (from its owner box).
void Pic::composite_fold(Vector<MultiFab>& mf, int nComp) {
  for (int iLev = finestLevel(); iLev > 0; iLev--) {
    const IntVect ratio = refRatio(iLev - 1);
    BoxArray bac = mf[iLev].boxArray();
    bac.coarsen(ratio);
    MultiFab ctmp(bac, mf[iLev].DistributionMap(), nComp, 0);
    ctmp.setVal(0.0);

    for (MFIter mfi(mf[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real const>& fine = mf[iLev][mfi].array();
      const Array4<Real>& crse = ctmp[mfi].array();
      const auto& status = node_status(iLev)[mfi].array();
      amrex::Loop(box, [&](int i, int j, int k) {
        if (!bit::is_lev_edge(status(i, j, k)) ||
            !bit::is_owner(status(i, j, k)))
          return;
        int ic, jc, kc;
        Real wx, wy, wz;
        coarse_parent(i, ratio[ix_], ic, wx);
        coarse_parent(j, ratio[iy_], jc, wy);
        coarse_parent(k, nDim > 2 ? ratio[iz_] : 1, kc, wz);
        for (int iVar = 0; iVar < nComp; iVar++) {
          const Real val = fine(i, j, k, iVar);
          for (int dk = 0; dk <= 1; dk++) {
            const Real wk = dk == 0 ? 1 - wz : wz;
            if (wk == 0)
              continue;
            for (int dj = 0; dj <= 1; dj++) {
              const Real wj = dj == 0 ? 1 - wy : wy;
              if (wj == 0)
                continue;
              for (int di = 0; di <= 1; di++) {
                const Real wi = di == 0 ? 1 - wx : wx;
                if (wi == 0)
                  continue;
                crse(ic + di, jc + dj, kc + dk, iVar) += wi * wj * wk * val;
              }
            }
          }
        }
      });
    }
    mf[iLev - 1].ParallelAdd(ctmp, 0, 0, nComp, IntVect(0), IntVect(0),
                             Geom(iLev - 1).periodicity());
  }
}

//==========================================================
// Leaf-cell mask, unknown-node mask and node volumes of the current grid.
void Pic::composite_setup() {
  std::string nameFunc = "Pic::composite_setup";
  timing_func(nameFunc);

  // The loops of this file run on the host (amrex::Loop): several of them
  // sum into a scalar or add to shared coarse nodes.
  if (Gpu::inLaunchRegion())
    Abort("Error: #FIELDSOLVER implicit1 is not available on GPUs yet.");

  if (fsolver.solverBackend != ESolverBackend::LegacyGMRES)
    Abort("Error: #FIELDSOLVER implicit1 does not support #EFIELDSOLVER with "
          "solverType amrex.");

  const int nLev = n_lev();

  // Nothing to do if the grid has not changed.
  bool isCurrent =
      static_cast<int>(compVol.size()) == n_lev_max() && nLevComposite_ == nLev;
  for (int iLev = 0; isCurrent && iLev < nLev; iLev++) {
    isCurrent = !compVol[iLev].empty() &&
                compVol[iLev].boxArray() == node_box_array(iLev) &&
                compVol[iLev].DistributionMap() == DistributionMap(iLev);
  }
  if (isCurrent)
    return;
  nLevComposite_ = nLev;

  compLeaf.resize(n_lev_max());
  compUnknown.resize(n_lev_max());
  compVol.resize(n_lev_max());
  compX.resize(n_lev_max());
  compT.resize(n_lev_max());
  compM.resize(n_lev_max());
  compCurl.resize(n_lev_max());

  // Non-periodic boundaries are treated as in the level-by-level solver: the
  // cells and nodes outside the domain are filled from the boundary conditions
  // and take part in the curl operators.
  isCompositeBounded = !Geom(0).isAllPeriodic();

  // A 'linetied' body drops all its nodes from the linear system and a
  // 'conducting' body the interior ones, as in node_system_layout.
  const bool skipAllBody = is_body_linetied();
  const bool skipInterior = is_body_conducting();

  nCompUnknown = 0;
  for (int iLev = 0; iLev < nLev; iLev++) {
    const BoxArray& ba = boxArray(iLev);
    const BoxArray& nba = node_box_array(iLev);
    const DistributionMapping& dm = DistributionMap(iLev);
    compLeaf[iLev].define(ba, dm, 1, 1);
    compUnknown[iLev].define(nba, dm, 1, 0);
    compVol[iLev].define(nba, dm, 1, 0);
    compX[iLev].define(nba, dm, 3, get_n_ghost());
    compT[iLev].define(nba, dm, 3, 0);
    compM[iLev].define(nba, dm, 3, 0);
    compCurl[iLev].define(ba, dm, 3, 1);
    compCurl[iLev].setVal(0.0);

    // 1 in the leaf cells, 0 in the refined cells and outside the level.
    compLeaf[iLev].setVal(0.0);
    for (MFIter mfi(compLeaf[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& leaf = compLeaf[iLev][mfi].array();
      const auto& status = cell_status(iLev)[mfi].array();
      amrex::Loop(box, [&](int i, int j, int k) {
        leaf(i, j, k) = bit::is_refined(status(i, j, k)) ? 0.0 : 1.0;
      });
    }
    compLeaf[iLev].FillBoundary(Geom(iLev).periodicity());

    if (isCompositeBounded) {
      // A cell outside the domain is a leaf cell if the cell of the domain
      // next to it is, so that the boundary nodes get the full stencil and
      // volume.
      const Box& domain = Geom(iLev).Domain();
      int nTouch = 0;
      for (MFIter mfi(compLeaf[iLev]); mfi.isValid(); ++mfi) {
        const Box box = amrex::grow(mfi.validbox(), 1);
        const Array4<Real>& leaf = compLeaf[iLev][mfi].array();
        const auto& status = cell_status(iLev)[mfi].array();
        amrex::Loop(box, [&](int i, int j, int k) {
          if (bit::is_domain_boundary(status(i, j, k))) {
            IntVect iv(AMREX_D_DECL(i, j, k));
            for (int iDim = 0; iDim < nDim; iDim++)
              iv[iDim] = std::min(std::max(iv[iDim], domain.smallEnd(iDim)),
                                  domain.bigEnd(iDim));
            leaf(i, j, k) = leaf(iv);
            if (iLev > 0)
              nTouch++;
          }
        });
      }
      // The boundary conditions are applied on all levels below, but a
      // refined level at a non-periodic boundary could not be tested: the
      // initial magnetic field is not a number there with any field solver.
      ParallelDescriptor::ReduceIntSum(nTouch);
      if (nTouch > 0)
        Abort("Error: #FIELDSOLVER implicit1 does not support refined levels "
              "that touch a non-periodic boundary.");
    }

    // Volume of the leaf cells of this level around the node.
    const Real vol = get_cell_volume(iLev);
    for (MFIter mfi(compVol[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real const>& leaf = compLeaf[iLev][mfi].array();
      const Array4<Real>& volArr = compVol[iLev][mfi].array();
      const Array4<int>& unk = compUnknown[iLev][mfi].array();
      const auto& status = node_status(iLev)[mfi].array();
      amrex::Loop(box, [&](int i, int j, int k) {
        Real sum = 0;
        for (int dk = kLoCell; dk <= 0; dk++)
          for (int dj = -1; dj <= 0; dj++)
            for (int di = -1; di <= 0; di++)
              sum += leaf(i + di, j + dj, k + dk);
        volArr(i, j, k) = sum * vol / nCellPerNode;
        const bool isEdge = iLev > 0 && bit::is_lev_edge(status(i, j, k));
        const int st = status(i, j, k);
        const bool isBody = (skipAllBody && bit::is_body(st)) ||
                            (skipInterior && bit::is_body_interior(st));
        unk(i, j, k) =
            (bit::is_owner(st) && sum > 0 && !isEdge && !isBody) ? 1 : 0;
        if (unk(i, j, k) == 1)
          nCompUnknown++;
      });
    }
  }
  composite_fold(compVol, 1);

  if (isCompositeBounded) {
    // Node volumes without the cells outside the domain, for the electric
    // energy. They agree with the cell based sum of calc_E_field_energy.
    compVolE.resize(n_lev_max());
    for (int iLev = 0; iLev < nLev; iLev++) {
      compVolE[iLev].define(node_box_array(iLev), DistributionMap(iLev), 1, 0);
      const Real vol = get_cell_volume(iLev);
      for (MFIter mfi(compVolE[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real const>& leaf = compLeaf[iLev][mfi].array();
        const Array4<Real>& volArr = compVolE[iLev][mfi].array();
        const auto& status = cell_status(iLev)[mfi].array();
        amrex::Loop(box, [&](int i, int j, int k) {
          Real sum = 0;
          for (int dk = kLoCell; dk <= 0; dk++)
            for (int dj = -1; dj <= 0; dj++)
              for (int di = -1; di <= 0; di++)
                if (!bit::is_domain_boundary(status(i + di, j + dj, k + dk)))
                  sum += leaf(i + di, j + dj, k + dk);
          volArr(i, j, k) = sum * vol / nCellPerNode;
        });
      }
    }
    composite_fold(compVolE, 1);
  }

  if (isCompositeBounded) {
    // Rows that are replaced by the identity, as in update_E_matvec: 1 for
    // all components on the faces with a prescribed electric field (inflow,
    // fixed, wave), 2 for the tangential components on conducting walls.
    compFixed.resize(n_lev_max());
    compIn.resize(n_lev_max());
    for (int iLev = 0; iLev < nLev; iLev++) {
      compFixed[iLev].define(node_box_array(iLev), DistributionMap(iLev), 3, 0);
      compIn[iLev].define(node_box_array(iLev), DistributionMap(iLev), 3, 0);
      const BoundaryBounds bnd(Geom(iLev), node_box_array(iLev).ixType(),
                               &bcField);
      for (MFIter mfi(compFixed[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<int>& fixed = compFixed[iLev][mfi].array();
        amrex::Loop(box, 3, [&](int i, int j, int k, int iVar) {
          const int ijk[3] = { i, j, k };
          fixed(i, j, k, iVar) = 0;
          if (hasConductingBC_ && is_conducting_wall_Et(bnd, ijk, iVar))
            fixed(i, j, k, iVar) = 2;
          for (int d = 0; d < nDim; ++d) {
            if (!bnd.isNode[d])
              continue;
            if ((FieldBC::is_dirichlet(bnd.bcLo[d]) &&
                 ijk[d] == bnd.loBnd[d]) ||
                (FieldBC::is_dirichlet(bnd.bcHi[d]) && ijk[d] == bnd.hiBnd[d]))
              fixed(i, j, k, iVar) = 1;
          }
        });
      }
    }
  }
}

//==========================================================
// Replace the rows of a level that are fixed by the boundary conditions:
// 'dirichlet' on the faces with a prescribed electric field and 'wall' for
// the tangential components on conducting walls.
void Pic::composite_set_fixed_rows(MultiFab& out, const MultiFab& dirichlet,
                                   const MultiFab& wall, const int iLev) {
  for (MFIter mfi(out); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arr = out[mfi].array();
    const Array4<Real const>& dir = dirichlet[mfi].array();
    const Array4<Real const>& wal = wall[mfi].array();
    const Array4<int const>& fixed = compFixed[iLev][mfi].array();
    amrex::Loop(box, 3, [&](int i, int j, int k, int iVar) {
      if (fixed(i, j, k, iVar) == 1) {
        arr(i, j, k, iVar) = dir(i, j, k, iVar);
      } else if (fixed(i, j, k, iVar) == 2) {
        arr(i, j, k, iVar) = wal(i, j, k, iVar);
      }
    });
  }
}

//==========================================================
// The 1D vector of the unknowns -> level arrays holding the composite field
// on all the corner nodes of the leaf cells (and zero elsewhere).
void Pic::composite_scatter(const double* const p, Vector<MultiFab>& mf) {
  int iCount = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    mf[iLev].setVal(0.0);
    for (MFIter mfi(mf[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arr = mf[iLev][mfi].array();
      const Array4<int const>& unk = compUnknown[iLev][mfi].array();
      const auto lo = lbound(box);
      const auto hi = ubound(box);
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i)
            if (unk(i, j, k) == 1)
              for (int iVar = 0; iVar < 3; iVar++)
                arr(i, j, k, iVar) = p[iCount++];
    }
    // Copies of the nodes shared by the boxes.
    mf[iLev].SumBoundary(Geom(iLev).periodicity());
    mf[iLev].FillBoundary(Geom(iLev).periodicity());
    if (is_fake_2d()) {
      Periodicity period(IntVect(AMREX_D_DECL(0, 0, 1)));
      mf[iLev].FillBoundary(period);
    }
  }
  composite_prolong(mf, 3);
}

//==========================================================
void Pic::composite_gather(const Vector<MultiFab>& mf, double* const p) {
  int iCount = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi(mf[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real const>& arr = mf[iLev][mfi].array();
      const Array4<int const>& unk = compUnknown[iLev][mfi].array();
      const auto lo = lbound(box);
      const auto hi = ubound(box);
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i)
            if (unk(i, j, k) == 1)
              for (int iVar = 0; iVar < 3; iVar++)
                p[iCount++] = arr(i, j, k, iVar);
    }
  }
}

//==========================================================
// out = V_lev * C^T (leaf * in), for a cell centered field on one level. The
// ghost cells of 'center' are overwritten.
void Pic::composite_curl_transpose(MultiFab& center, MultiFab& out, int iLev) {
  center.FillBoundary(Geom(iLev).periodicity());
  for (int iVar = 0; iVar < 3; iVar++)
    MultiFab::Multiply(center, compLeaf[iLev], 0, iVar, 1, 1);
  curl_center_to_node(center, out, Geom(iLev).InvCellSize());
  out.mult(get_cell_volume(iLev), 0, 3, 0);
}

//==========================================================
// vecOut = A vecIn for the composite system
// A = I + V^-1 [(theta dt)^2 C^T W_B C + 4 pi theta dt P^T M P],
// where W_B holds the leaf cell volumes and M is the unnormalised mass matrix.
// With non-periodic boundaries the nodes outside the domain are filled as in
// update_E_matvec: with the homogeneous part of the boundary conditions for
// a correction of E (useZeroBC), or with the boundary values for E itself.
void Pic::composite_matvec(const double* vecIn, double* vecOut,
                           const bool useZeroBC) {
  std::string nameFunc = "Pic::composite_matvec";
  timing_func(nameFunc);

  composite_scatter(vecIn, compX);

  // The wider stencils of ratioDivC2C and #UPWINDE also use the coarse nodes
  // covered by a finer level and the ghost nodes of the finer levels.
  const bool doDivC2C = fsolver.coefDiff > 0;
  if (isCompositeBounded) {
    for (int iLev = 0; iLev < n_lev(); iLev++)
      MultiFab::Copy(compIn[iLev], compX[iLev], 0, 0, 3, 0);
    if (useZeroBC) {
      // The correction is zero where the electric field is prescribed.
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        for (MFIter mfi(compX[iLev]); mfi.isValid(); ++mfi) {
          const Box& box = mfi.validbox();
          const Array4<Real>& arr = compX[iLev][mfi].array();
          const Array4<int const>& fixed = compFixed[iLev][mfi].array();
          amrex::Loop(box, 3, [&](int i, int j, int k, int iVar) {
            if (fixed(i, j, k, iVar) == 1)
              arr(i, j, k, iVar) = 0.0;
          });
        }
      }
    }
    composite_fill_ghost(compX, useZeroBC ? bcHomogeneous_ : bcMatvec_);
  } else if (doDivC2C || useUpwindE) {
    composite_fill_ghost(compX, bcNone_);
  }

  const Real delt2 = pow(fsolver.theta * tc->get_dt(), 2);
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    curl_node_to_center(compX[iLev], compCurl[iLev], Geom(iLev).InvCellSize());
    composite_curl_transpose(compCurl[iLev], compT[iLev], iLev);

    if (doDivC2C) {
      // The level solver uses -lap + grad(div) for curl curl, with div =
      // (1 - coefDiff) * divN2C + coefDiff * divC2C. Here curl curl is
      // C^T W_B C, so the difference coefDiff * grad(divC2C - divN2C) of the
      // leaf cells is added.
      MultiFab& divE = centerDivE[iLev];
      MultiFab& divN2C = compCurl[iLev];
      div_node_to_center(compX[iLev], divE, Geom(iLev).InvCellSize());
      MultiFab::Copy(divN2C, divE, 0, 0, 1, 1);
      blend_div_E_c2c(compX[iLev], iLev);
      MultiFab::Subtract(divE, divN2C, 0, 0, 1, 1);
      MultiFab::Multiply(divE, compLeaf[iLev], 0, 0, 1, 1);
      grad_center_to_node(divE, compM[iLev], Geom(iLev).InvCellSize());
      MultiFab::Saxpy(compT[iLev], get_cell_volume(iLev), compM[iLev], 0, 0, 3,
                      0);
    }
    compT[iLev].mult(delt2, 0, 3, 0);

    // nodeMM is divided by the cell volume, and update_E_M_dot_E multiplies
    // with 4 pi theta dt.
    update_E_M_dot_E(compX[iLev], compM[iLev], iLev);
    MultiFab::Saxpy(compT[iLev], get_cell_volume(iLev), compM[iLev], 0, 0, 3,
                    0);
  }
  composite_fold(compT, 3);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (int iVar = 0; iVar < 3; iVar++)
      MultiFab::Divide(compT[iLev], compVol[iLev], 0, iVar, 1, 0);
    MultiFab::Add(compT[iLev], compX[iLev], 0, 0, 3, 0);
    if (useUpwindE)
      add_upwind_E_matvec(compX[iLev], compT[iLev], iLev);
  }
  if (isCompositeBounded)
    for (int iLev = 0; iLev < n_lev(); iLev++)
      composite_set_fixed_rows(compT[iLev], compX[iLev], compIn[iLev], iLev);
  // 'conducting' body: keep the Krylov vectors in the constrained space.
  if (is_body_conducting())
    for (int iLev = 0; iLev < n_lev(); iLev++)
      project_body_E(compT[iLev], iLev);
  composite_gather(compT, vecOut);
}

//==========================================================
// rhs = E^n + V^-1 theta dt (C^T W_B B^n - 4 pi P^T jHat)
void Pic::composite_rhs(double* rhs) {
  const Real delt = fsolver.theta * tc->get_dt();
  if (doSmoothJ)
    composite_smooth_J();
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (isCompositeBounded || useCurlUpwindB) {
      // The magnetic field in the cells and on the nodes outside the domain,
      // as in update_E_rhs.
      nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
      if (useCurlUpwindB) {
        apply_centerB_BC(iLev);
      } else {
        centerB[iLev].FillBoundary(Geom(iLev).periodicity());
        if (iLev == 0) {
          apply_field_bc(cell_status(0), centerB[0], 0, centerB[0].nComp(),
                         &Pic::get_center_B, 0, true, -2.0);
        } else {
          fill_fine_lev_bny_from_coarse(
              centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
              refRatio(iLev - 1), Geom(iLev - 1), Geom(iLev), cell_status(iLev),
              *get_cell_interp());
        }
      }
      if (iLev == 0) {
        apply_field_bc(node_status(0), nodeB[0], 0, nodeB[0].nComp(),
                       &Pic::get_node_B, 0, true, 0.0);
      } else {
        fill_fine_lev_bny_from_coarse(
            nodeB[iLev - 1], nodeB[iLev], 0, nodeB[iLev - 1].nComp(),
            refRatio(iLev - 1), Geom(iLev - 1), Geom(iLev), node_status(iLev),
            node_bilinear_interp);
      }
      MultiFab::Copy(compCurl[iLev], centerB[iLev], 0, 0, 3, 1);
    } else {
      MultiFab::Copy(compCurl[iLev], centerB[iLev], 0, 0, 3, 0);
    }
    composite_curl_transpose(compCurl[iLev], compT[iLev], iLev);
    // jHat is divided by the cell volume. With #SMOOTHJ the smoothed current
    // is added below.
    if (!doSmoothJ)
      MultiFab::Saxpy(compT[iLev], -fourPI * get_cell_volume(iLev), jHat[iLev],
                      0, 0, 3, 0);
    compT[iLev].mult(delt, 0, 3, 0);

    // The currents of the background and of the inductive electric fields,
    // as in update_E_rhs. update_E_M_dot_E multiplies with 4 pi theta dt.
    if (solveFieldInCoMov) {
      update_E_M_dot_E(eBg[iLev], compM[iLev], iLev);
      MultiFab::Saxpy(compT[iLev], get_cell_volume(iLev), compM[iLev], 0, 0, 3,
                      0);
    }
    if (has_active_inductive_E0() && !nodeE0[iLev].empty()) {
      update_E_M_dot_E(nodeE0[iLev], compM[iLev], iLev);
      MultiFab::Saxpy(compT[iLev], -get_cell_volume(iLev), compM[iLev], 0, 0, 3,
                      0);
    }
  }
  composite_fold(compT, 3);
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (int iVar = 0; iVar < 3; iVar++)
      MultiFab::Divide(compT[iLev], compVol[iLev], 0, iVar, 1, 0);
    MultiFab::Add(compT[iLev], nodeE[iLev], 0, 0, 3, 0);
  }
  if (doSmoothJ) {
    for (int iLev = 0; iLev < n_lev(); iLev++)
      MultiFab::Saxpy(compT[iLev], -fourPI * delt, compJ[iLev], 0, 0, 3, 0);
  }
  if (isCompositeBounded) {
    // The right-hand side of the fixed rows is the electric field itself, so
    // that the residual rhs - A E^n is zero there: the prescribed field at
    // the time of the solve on the inflow faces, E^n otherwise.
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      const MultiFab* dirichletE = &nodeE[iLev];
      if (hasInflowBC_) {
        MultiFab::Copy(compX[iLev], nodeE[iLev], 0, 0, 3,
                       std::min(compX[iLev].nGrow(), nodeE[iLev].nGrow()));
        apply_field_bc(node_status(iLev), compX[iLev], 0, nDim3,
                       &Pic::get_node_E, iLev, false, 0.0);
        dirichletE = &compX[iLev];
      }
      composite_set_fixed_rows(compT[iLev], *dirichletE, nodeE[iLev], iLev);
    }
  }
  if (is_body_conducting())
    for (int iLev = 0; iLev < n_lev(); iLev++)
      project_body_E(compT[iLev], iLev);
  composite_gather(compT, rhs);
}

//==========================================================
// The level arrays of E from the unknowns: prolongation onto the fine edge
// nodes; the covered coarse nodes (not used by any leaf cell or particle) are
// injected from the finer level, and the ghost nodes outside a fine level are
// interpolated from the coarser one, both for output only.
void Pic::composite_fill_E(const double* const p, Vector<MultiFab>& mf) {
  composite_scatter(p, mf);
  composite_fill_ghost(mf, bcFinal_);
}

//==========================================================
// The nodes that are not unknowns and are not set by composite_scatter: the
// covered coarse nodes are injected from the finer level and the ghost nodes
// outside a fine level are interpolated from the coarser one. The boundary
// conditions of E are applied on every level as in the level-by-level solver:
// bcFinal_ for E after the solve (update_E_impl), bcMatvec_ for E in the
// matvec, bcHomogeneous_ for a correction of E (update_E_matvec), bcNone_
// for other fields.
void Pic::composite_fill_ghost(Vector<MultiFab>& mf, const int bcMode) {
  const bool doBC = isCompositeBounded && bcMode != bcNone_;
  for (int iLev = finestLevel(); iLev > 0; iLev--) {
    average_down_nodal(mf[iLev], mf[iLev - 1], refRatio(iLev - 1));
  }
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (iLev > 0 && doBC) {
      // The values set on the boundary faces of the coarser level are
      // interpolated to the fine edge nodes.
      composite_prolong(mf, 3);
    }
    mf[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev > 0) {
      fill_fine_lev_bny_from_coarse(
          mf[iLev - 1], mf[iLev], 0, 3, refRatio(iLev - 1), Geom(iLev - 1),
          Geom(iLev), node_status(iLev), node_bilinear_interp);
      mf[iLev].FillBoundary(Geom(iLev).periodicity());
    }
    if (!doBC)
      continue;
    if (bcMode == bcHomogeneous_) {
      apply_field_bc_homogeneous(node_status(iLev), mf[iLev], iLev);
    } else if (bcMode == bcMatvec_) {
      apply_field_bc(node_status(iLev), mf[iLev], 0, nDim3, &Pic::get_node_E,
                     iLev, false, 0.0);
      // The boundary values written on the faces are copied to the ghost
      // nodes of the neighbouring boxes.
      mf[iLev].FillBoundary(Geom(iLev).periodicity());
    } else {
      // The wave source is applied inside apply_field_bc().
      apply_field_bc(node_status(iLev), mf[iLev], 0, nDim3, &Pic::get_node_E,
                     iLev, false);
    }
  }
}

//==========================================================
// compJ := the current of #SMOOTHJ: jHat with the fine edge nodes folded
// onto the coarse nodes, as a field on all nodes, smoothed level by level.
void Pic::composite_smooth_J() {
  const int nLev = n_lev();
  compJ.resize(n_lev_max());
  for (int iLev = 0; iLev < nLev; iLev++) {
    if (compJ[iLev].empty() || compJ[iLev].boxArray() != node_box_array(iLev) ||
        compJ[iLev].DistributionMap() != DistributionMap(iLev))
      compJ[iLev].define(node_box_array(iLev), DistributionMap(iLev), 3,
                         jHat[iLev].nGrow());
    MultiFab::Copy(compT[iLev], jHat[iLev], 0, 0, 3, 0);
    compT[iLev].mult(get_cell_volume(iLev), 0, 3, 0);
  }
  composite_fold(compT, 3);
  for (int iLev = 0; iLev < nLev; iLev++)
    for (int iVar = 0; iVar < 3; iVar++)
      MultiFab::Divide(compT[iLev], compVol[iLev], 0, iVar, 1, 0);

  // The folded current on the unknown nodes and the fine edge nodes. The
  // other nodes (inside a body, outside the domain) keep the deposited
  // values; compM marks the nodes set by composite_scatter.
  std::vector<double> jVec(3 * static_cast<size_t>(nCompUnknown));
  std::vector<double> ones(jVec.size(), 1.0);
  composite_scatter(ones.data(), compX);
  for (int iLev = 0; iLev < nLev; iLev++)
    MultiFab::Copy(compM[iLev], compX[iLev], 0, 0, 3, 0);
  composite_gather(compT, jVec.data());
  composite_scatter(jVec.data(), compX);
  for (int iLev = 0; iLev < nLev; iLev++) {
    MultiFab::Copy(compJ[iLev], jHat[iLev], 0, 0, 3, compJ[iLev].nGrow());
    for (MFIter mfi(compJ[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& j = compJ[iLev][mfi].array();
      const Array4<Real const>& x = compX[iLev][mfi].array();
      const Array4<Real const>& mark = compM[iLev][mfi].array();
      amrex::Loop(box, 3, [&](int i, int jj, int k, int iVar) {
        if (mark(i, jj, k, iVar) > 0.5)
          j(i, jj, k, iVar) = x(i, jj, k, iVar);
      });
    }
  }
  composite_fill_ghost(compJ, bcNone_);

  for (int iLev = 0; iLev < nLev; iLev++) {
    for (int icount = 0; icount < nSmoothJ; icount++) {
      if (iLev > 0) {
        fill_fine_lev_bny_from_coarse(compJ[iLev - 1], compJ[iLev], 0, 3,
                                      refRatio(iLev - 1), Geom(iLev - 1),
                                      Geom(iLev), node_status(iLev),
                                      node_bilinear_interp);
        compJ[iLev].FillBoundary(Geom(iLev).periodicity());
      }
      smooth_multifab(compJ[iLev], iLev, icount % 2 + 1, coefSmoothJ);
    }
  }
}

//==========================================================
void Pic::update_E_composite() {
  std::string nameFunc = "Pic::update_E_composite";
  timing_func(nameFunc);

  composite_setup();

  if (doTestComposite && !isCompositeTested) {
    composite_self_test();
    isCompositeTested = true;
  }

  eSolver.reset(nCompUnknown);
  const int nSolve = eSolver.get_nSolve();

  // Energy budget diagnostics in the first steps only (they loop over the
  // particles).
  const bool doBudget = doTestComposite && tc->get_cycle() < 3;

  composite_rhs(eSolver.rhs);

  // E^n on the unknowns. The values on the other nodes do not enter.
  std::vector<double> eOld(nSolve);
  composite_gather(nodeE, eOld.data());

  // Solve A dE = rhs - A E^n for the change dE = E^{n+theta} - E^n.
  // E^n is used with its boundary values, the change dE with the homogeneous
  // boundary conditions.
  // The solvers pass their level index to update_E_matvec. The value
  // iLevComposite_ (-1) is not an AMR level: it selects composite_matvec,
  // which acts on the unknowns of all levels.
  if (fsolver.mode == FieldSolverMode::NewtonKrylov) {
    solve_E_newton_krylov(iLevComposite_);
  } else {
    composite_matvec(eOld.data(), eSolver.matvec, false);
    for (int i = 0; i < nSolve; ++i) {
      eSolver.rhs[i] -= eSolver.matvec[i];
      eSolver.xLeft[i] = 0;
    }

    if (doReport)
      Print() << "\n-------" << printPrefix
              << " GMRES composite E solver ------------------" << std::endl;

    BL_PROFILE_VAR("Pic::E_iterate", eSolverComp);
    eSolver.solve(iLevComposite_, doReport);
    BL_PROFILE_VAR_STOP(eSolverComp);
  }

  // E^{n+theta} and E^{n+1} = (E^{n+theta} - (1 - theta) E^n) / theta
  std::vector<double> eNew(nSolve);
  for (int i = 0; i < nSolve; ++i) {
    const double eTh = eOld[i] + eSolver.xLeft[i];
    eNew[i] = (eTh - (1.0 - fsolver.theta) * eOld[i]) / fsolver.theta;
    eOld[i] = eTh;
  }
  if (doBudget) {
    // Energy budget of the step from the grid side: the work on the
    // particles dt E_theta.(jHat + M E_theta) and the change of the electric
    // energy. eOld holds E_theta here.
    const MPI_Comm iComm = ParallelDescriptor::Communicator();
    std::vector<double> work(nSolve), e0(nSolve), vol(nSolve);
    composite_gather(nodeE, e0.data());
    composite_scatter(eOld.data(), compX);
    const Real c0 = fourPI * fsolver.theta * tc->get_dt();
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      update_E_M_dot_E(compX[iLev], compM[iLev], iLev);
      MultiFab::LinComb(compT[iLev], get_cell_volume(iLev) / c0, compM[iLev], 0,
                        get_cell_volume(iLev), jHat[iLev], 0, 0, 3, 0);
    }
    // The same work summed level by level over all nodes (owner copies).
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      Real sumLev = 0;
      for (MFIter mfi(compT[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real const>& t = compT[iLev][mfi].array();
        const Array4<Real const>& x = compX[iLev][mfi].array();
        const auto& status = node_status(iLev)[mfi].array();
        amrex::Loop(box, 3, [&](int i, int j, int k, int iVar) {
          if (bit::is_owner(status(i, j, k)))
            sumLev += t(i, j, k, iVar) * x(i, j, k, iVar);
        });
      }
      ParallelDescriptor::ReduceRealSum(sumLev);
      Print() << printPrefix << "Composite budget: level " << iLev
              << " dt sum E.J = " << std::setprecision(12)
              << sumLev * tc->get_dt() << std::endl;
    }
    composite_fold(compT, 3);
    composite_gather(compT, work.data());
    const Real jDotE =
        tc->get_dt() * dot_product_mpi(eOld.data(), work.data(), nSolve, iComm);
    for (int iLev = 0; iLev < n_lev(); iLev++)
      for (int iVar = 0; iVar < 3; iVar++)
        MultiFab::Copy(compT[iLev], compVol[iLev], 0, iVar, 1, 0);
    composite_gather(compT, vol.data());
    Real wE0 = 0, wE1 = 0;
    for (int i = 0; i < nSolve; ++i) {
      wE0 += vol[i] * e0[i] * e0[i];
      wE1 += vol[i] * eNew[i] * eNew[i];
    }
    ParallelDescriptor::ReduceRealSum(wE0);
    ParallelDescriptor::ReduceRealSum(wE1);
    Print() << printPrefix
            << "Composite budget: dt E.J = " << std::setprecision(12) << jDotE
            << " dW_E = " << (wE1 - wE0) * 0.5 / fourPI << std::endl;
  }

  if (doSmoothE) {
    // As in update_E_impl: the change of E is smoothed first, then E at the
    // two time levels after the boundary conditions.
    composite_scatter(eSolver.xLeft, compX);
    composite_fill_ghost(compX, bcNone_);
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      smooth_E(compX[iLev], iLev);
      MultiFab::Copy(nodeEth[iLev], compX[iLev], 0, 0, 3, get_n_ghost());
      MultiFab::Add(nodeEth[iLev], nodeE[iLev], 0, 0, 3, get_n_ghost());
      MultiFab::LinComb(nodeE[iLev], -(1.0 - fsolver.theta) / fsolver.theta,
                        nodeE[iLev], 0, 1. / fsolver.theta, nodeEth[iLev], 0, 0,
                        3, get_n_ghost());
    }
    for (int iStage = 0; iStage < 2; iStage++) {
      // The fine edge nodes follow the coarse nodes again.
      composite_prolong(nodeE, 3);
      composite_prolong(nodeEth, 3);
      composite_fill_ghost(nodeE, bcFinal_);
      composite_fill_ghost(nodeEth, bcFinal_);
      if (iStage == 0) {
        for (int iLev = 0; iLev < n_lev(); iLev++) {
          smooth_E(nodeEth[iLev], iLev);
          smooth_E(nodeE[iLev], iLev);
        }
      }
      if (n_lev() == 1)
        break;
    }
  } else {
    // Same order as in update_E_impl: the boundary conditions may depend on
    // it.
    composite_fill_E(eNew.data(), nodeE);
    composite_fill_E(eOld.data(), nodeEth);
  }

  // Electric field boundary on the inner body (see #BODYBOUNDARY).
  if (use_body()) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      apply_body_E_bc(nodeE[iLev], iLev);
      apply_body_E_bc(nodeEth[iLev], iLev);
    }
  }

  if (doBudget) {
    // Change of the kinetic energy that the mover will produce with nodeEth
    // and nodeB, level by level.
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      Real dK = 0, kin = 0, qSum = 0;
      Long nPart = 0;
      for (int iSp = 0; iSp < nSpecies; iSp++) {
        nPart += parts[iSp]->NumberOfParticlesAtLevel(iLev, true, true);
        const Real qom = parts[iSp]->get_charge() / parts[iSp]->get_mass();
        const Real qdto2mc = qom * 0.5 * tc->get_dt();
        for (PicParticles::PIter pti(*parts[iSp], iLev); pti.isValid(); ++pti) {
          const Array4<Real const>& eArr = nodeEth[iLev][pti].array();
          const Array4<Real const>& bArr = nodeB[iLev][pti].array();
          for (const auto& p : pti.GetArrayOfStructs()) {
            if (p.id() < 0)
              continue;
            IntVect loIdx;
            RealVect dShift;
            Real coef[2][2][2];
            find_node_interpolation(p.pos(), Geom(iLev).ProbLo(),
                                    Geom(iLev).InvCellSize(), loIdx, dShift,
                                    coef);
            Real ep[3] = { 0, 0, 0 }, bp[3] = { 0, 0, 0 };
            const int kLo = nDim > 2 ? loIdx[nDim - 1] : 0;
            for (int kk = 0; kk <= (nDim > 2 ? 1 : 0); kk++)
              for (int jj = 0; jj <= 1; jj++)
                for (int ii = 0; ii <= 1; ii++)
                  for (int iVar = 0; iVar < 3; iVar++) {
                    ep[iVar] +=
                        coef[ii][jj][kk] *
                        eArr(loIdx[ix_] + ii, loIdx[iy_] + jj, kLo + kk, iVar);
                    bp[iVar] +=
                        coef[ii][jj][kk] *
                        bArr(loIdx[ix_] + ii, loIdx[iy_] + jj, kLo + kk, iVar);
                  }
            const Real v[3] = { p.rdata(PicParticles::iup_),
                                p.rdata(PicParticles::iup_ + 1),
                                p.rdata(PicParticles::iup_ + 2) };
            const Real ve[3] = { v[0] + qdto2mc * ep[0], v[1] + qdto2mc * ep[1],
                                 v[2] + qdto2mc * ep[2] };
            const Real om[3] = { qdto2mc * bp[0], qdto2mc * bp[1],
                                 qdto2mc * bp[2] };
            Real v1[3];
            boris_push_nonrelativistic(v, ve, om, v1);
            const Real mp = p.rdata(PicParticles::iqp_) / qom;
            kin += 0.5 * mp * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            qSum += mp;
            dK += 0.5 * mp *
                  (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] - v[0] * v[0] -
                   v[1] * v[1] - v[2] * v[2]);
          }
        }
      }
      ParallelDescriptor::ReduceRealSum(dK);
      ParallelDescriptor::ReduceRealSum(kin);
      ParallelDescriptor::ReduceRealSum(qSum);
      ParallelDescriptor::ReduceLongSum(nPart);
      Print() << printPrefix << "Composite budget: level " << iLev
              << " dK of particles = " << std::setprecision(12) << dK
              << " K = " << kin << " mass = " << qSum << " nPart = " << nPart
              << std::endl;
    }
  }

  if (doBudget) {
    // nodeEth must equal the composite field on all corners of leaf cells.
    composite_scatter(eOld.data(), compX);
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      Real diffMax = 0;
      for (MFIter mfi(compX[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real const>& a = nodeEth[iLev][mfi].array();
        const Array4<Real const>& x = compX[iLev][mfi].array();
        const Array4<Real const>& leaf = compLeaf[iLev][mfi].array();
        amrex::Loop(box, 3, [&](int i, int j, int k, int iVar) {
          Real sum = 0;
          for (int dk = kLoCell; dk <= 0; dk++)
            for (int dj = -1; dj <= 0; dj++)
              for (int di = -1; di <= 0; di++)
                sum += leaf(i + di, j + dj, k + dk);
          if (sum > 0)
            diffMax = std::max(diffMax,
                               std::abs(a(i, j, k, iVar) - x(i, j, k, iVar)));
        });
      }
      ParallelDescriptor::ReduceRealMax(diffMax);
      Print() << printPrefix << "Composite budget: level " << iLev
              << " max|nodeEth - P E| = " << diffMax << std::endl;
    }
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    div_node_to_center(nodeE[iLev], centerDivE[iLev], Geom(iLev).InvCellSize());
  }
}

//==========================================================
// sum_g V_g E_g^2 / 8 pi over the unknown nodes.
Real Pic::calc_E_field_energy_composite() {
  composite_setup();

  // The inductive electric field of a time dependent intrinsic magnetic
  // field is part of the total field, as in calc_E_field_energy.
  const bool hasActiveE0 = has_active_inductive_E0();
  Real sum = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    const bool addE0 = hasActiveE0 && !nodeE0[iLev].empty();
    for (MFIter mfi(nodeE[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real const>& arr = nodeE[iLev][mfi].array();
      const Array4<Real const> arrE0 =
          addE0 ? nodeE0[iLev][mfi].const_array() : Array4<Real const>();
      const Array4<Real const>& vol = isCompositeBounded
                                          ? compVolE[iLev][mfi].array()
                                          : compVol[iLev][mfi].array();
      const Array4<int const>& unk = compUnknown[iLev][mfi].array();
      amrex::Loop(box, [&](int i, int j, int k) {
        if (unk(i, j, k) == 1) {
          Real e2 = 0;
          for (int iVar = 0; iVar < 3; iVar++) {
            const Real e =
                arr(i, j, k, iVar) + (addE0 ? arrE0(i, j, k, iVar) : 0.0);
            e2 += e * e;
          }
          sum += vol(i, j, k) * e2;
        }
      });
    }
  }
  sum *= 0.5 / fourPI;
  ParallelDescriptor::ReduceRealSum(sum,
                                    ParallelDescriptor::IOProcessorNumber());
  if (!ParallelDescriptor::IOProcessor())
    sum = 0;
  return sum;
}

//==========================================================
// Checks of the composite operators with random vectors:
//  1. sum of the node volumes = volume of the domain;
//  2. P^T is the transpose of P and the Ampere curl is the transpose of the
//     Faraday curl: u.(C^T W_B C v) = v.(C^T W_B C u) = sum_leaf V (Cu).(Cv);
//  3. the mass matrix part is the transpose pair of gather and deposit only
//     through the particles, which the energy of the run checks.
void Pic::composite_self_test() {
  const MPI_Comm iComm = ParallelDescriptor::Communicator();
  const int n = 3 * nCompUnknown;
  std::vector<double> u(n), v(n), su(n), sv(n);
  for (int i = 0; i < n; ++i) {
    u[i] = amrex::Random() - 0.5;
    v[i] = amrex::Random() - 0.5;
  }

  // Volume
  Real volSum = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi(compVol[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real const>& vol = compVol[iLev][mfi].array();
      const Array4<int const>& unk = compUnknown[iLev][mfi].array();
      amrex::Loop(box, [&](int i, int j, int k) {
        if (unk(i, j, k) == 1)
          volSum += vol(i, j, k);
      });
    }
  }
  ParallelDescriptor::ReduceRealSum(volSum);
  Real volDomain = 1;
  for (int iDim = 0; iDim < nDim; iDim++)
    volDomain *= Geom(0).ProbLength(iDim);

  // S w = fold(V C^T leaf C scatter(w)); also returns sum_leaf V (C w1).(C w2)
  Vector<MultiFab> curlSave(n_lev());
  auto apply_s = [&](const std::vector<double>& in, std::vector<double>& out,
                     bool doSave, Real& cDotC) {
    composite_scatter(in.data(), compX);
    cDotC = 0;
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      curl_node_to_center(compX[iLev], compCurl[iLev],
                          Geom(iLev).InvCellSize());
      if (doSave) {
        curlSave[iLev].define(boxArray(iLev), DistributionMap(iLev), 3, 0);
        MultiFab::Copy(curlSave[iLev], compCurl[iLev], 0, 0, 3, 0);
      } else {
        for (MFIter mfi(curlSave[iLev]); mfi.isValid(); ++mfi) {
          const Box& box = mfi.validbox();
          const Array4<Real const>& a = curlSave[iLev][mfi].array();
          const Array4<Real const>& b = compCurl[iLev][mfi].array();
          const Array4<Real const>& leaf = compLeaf[iLev][mfi].array();
          const Real vol = get_cell_volume(iLev);
          amrex::Loop(box, 3, [&](int i, int j, int k, int iVar) {
            cDotC += vol * leaf(i, j, k) * a(i, j, k, iVar) * b(i, j, k, iVar);
          });
        }
      }
      composite_curl_transpose(compCurl[iLev], compT[iLev], iLev);
    }
    composite_fold(compT, 3);
    composite_gather(compT, out.data());
  };

  Real dummy = 0, cuDotCv = 0;
  apply_s(u, su, true, dummy);
  apply_s(v, sv, false, cuDotCv);
  ParallelDescriptor::ReduceRealSum(cuDotCv);
  const Real uSv = dot_product_mpi(u.data(), sv.data(), n, iComm);
  const Real vSu = dot_product_mpi(v.data(), su.data(), n, iComm);

  Long nUnknownAll = nCompUnknown;
  ParallelDescriptor::ReduceLongSum(nUnknownAll);

  Print() << printPrefix
          << "Composite solver test: unknown nodes = " << nUnknownAll << "\n"
          << printPrefix
          << "  sum(V_g)/V_domain - 1    = " << volSum / volDomain - 1 << "\n"
          << printPrefix << "  u.S v = " << uSv << ", v.S u = " << vSu
          << ", sum_leaf V (Cu).(Cv) = " << cuDotCv << "\n"
          << printPrefix
          << "  (u.Sv - v.Su)/|u.Sv|     = " << (uSv - vSu) / std::abs(uSv)
          << "\n"
          << printPrefix
          << "  (u.Sv - (Cu).(Cv))/|u.Sv| = " << (uSv - cuDotCv) / std::abs(uSv)
          << std::endl;
}

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PhysBCFunct.H>

#include "GridUtility.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

namespace {

struct AbsorbWeights {
  Real decay[3] = { 0.0, 0.0, 0.0 };
  Real drive[3] = { 0.0, 0.0, 0.0 };
};

inline BoxArray get_boundary_active_ba(const BoxArray& activeRegion,
                                       const MultiFab& mf, const Geometry& geom,
                                       int nDimVal, int iz) {
  BoxArray ba = convert(activeRegion, mf.boxArray().ixType());
  const IntVect& ngrow = mf.nGrowVect();
  if (nDimVal > 2 && geom.Domain().bigEnd(iz) == geom.Domain().smallEnd(iz)) {
    ba.grow(iz, ngrow[iz]);
  }
  return ba;
}

// The active region of the level of geom in the index type of mf: activeRegion
// is in level-0 index space (geom0), so it is refined to that level first.
inline BoxArray get_level_active_ba(const BoxArray& activeRegion,
                                    const MultiFab& mf, const Geometry& geom,
                                    const Geometry& geom0, int nDimVal,
                                    int iz) {
  IntVect ratio(1);
  for (int d = 0; d < AMREX_SPACEDIM; ++d)
    ratio[d] = geom.Domain().length(d) / geom0.Domain().length(d);
  return get_boundary_active_ba(amrex::refine(activeRegion, ratio), mf, geom,
                                nDimVal, iz);
}

// True if the level-boundary (ghost) cell/node iv of a level lies outside
// the active region ba of that level (see get_level_active_ba), once it is
// mapped across the periodic boundaries onto the cell/node it is an image of.
// Only those ghosts are outer boundary ghosts. On a refined level the
// coarse-fine interface ghosts are level boundary too, but they lie inside the
// active region: they belong to the coarse solution. (The domain-boundary
// status bit cannot tell them apart: it is also set on ghost nodes lying
// exactly on a domain face, e.g. on a periodic seam or on both z layers of a
// fake-2D grid.)
AMREX_GPU_HOST_DEVICE inline bool is_outside_active(const BoxListView& ba,
                                                    const int* periodic,
                                                    const int* domainLen,
                                                    const BoundaryBounds& bnd,
                                                    int nDimVal, IntVect iv) {
  for (int d = 0; d < nDimVal; ++d) {
    if (!periodic[d])
      continue;
    if (iv[d] < bnd.loBnd[d])
      iv[d] += domainLen[d];
    else if (iv[d] > bnd.hiBnd[d])
      iv[d] -= domainLen[d];
  }
  return !ba.contains(iv);
}

struct ActiveComp {
  int dir;
  int side;
  Real effAmp;
  Real omega_t;
  Real k_vec[3];
  Real pol[3];
  int profile;
  Real tCenter;
  Real tInvWidth;
  const WaveComponent* orig;
};

} // namespace

//==========================================================
void Pic::apply_field_bc(const iMultiFab& status, MultiFab& mf,
                         const int iStart, const int nComp, GETVALUE func,
                         const int iLev, const bool isB, Real dtOverride,
                         Real timeOverride) {
  if (Geom(iLev).isAllPeriodic() || mf.nGrow() == 0)
    return;

  std::string nameFunc = "Pic::apply_field_bc";
  timing_func(nameFunc);

  const auto qty =
      isB ? FieldBC::Quantity::Magnetic : FieldBC::Quantity::Electric;
  const Vector<BCRec> bcr =
      FieldBC::create_bcrec(bcField, qty, nComp, Geom(iLev));

  // Native AMReX physical boundary fill for foextrap and reflections:
  GpuBndryFuncFab<FabFillNoOp> bfunc(FabFillNoOp{});
  PhysBCFunct<GpuBndryFuncFab<FabFillNoOp> > physbcf(Geom(iLev), bcr, bfunc);
  physbcf(mf, iStart, nComp, mf.nGrowVect(), 0.0, 0);

  const Real t =
      (timeOverride >= 0.0) ? timeOverride : (tc ? tc->get_time() : 0.0);

  // Fill external Dirichlet (coupled / fixed) faces from func if provided:
  // Hybrid PIC fills (or, for the first B layer, Faraday-advances) the
  // coarse-fine ghosts of a refined level itself, so the prescribed state must
  // stay out of them. The full-PIC solver still relies on this fill there.
  if (func != nullptr) {
    fill_ext_dir(status, mf, iStart, nComp, func, iLev, &bcField,
                 iLev > 0 || useHybridPIC || useCurlUpwindB, t);
  }

  // Node-centered conducting wall requires zeroing boundary wall nodes
  // and mirroring nodal ghost nodes.
  if (hasConductingBC_ && mf.boxArray().ixType().nodeCentered()) {
    apply_conducting_wall(status, mf, iStart, nComp, iLev, bcField, isB);
  }

  if (hasInflowBC_) {
    apply_inflow_wall(status, mf, iStart, nComp, iLev, bcField, isB, func, t);
  }

  if (hasAbsorbBC_)
    apply_absorbing_wall(status, mf, iStart, nComp, iLev, bcField, isB,
                         dtOverride);

  // Wave boundary condition overwrites faces where active.
  if (waveBC.active) {
    apply_wave_field(status, mf, iStart, nComp, iLev, bcField, isB ? 0 : 1, t,
                     func);
  }
}

//==========================================================
void Pic::apply_field_bc_homogeneous(const iMultiFab& status, MultiFab& mf,
                                     const int iLev) {
  if (Geom(iLev).isAllPeriodic() || mf.nGrow() == 0)
    return;

  std::string nameFunc = "Pic::apply_field_bc_homogeneous";
  timing_func(nameFunc);

  const int nComp = mf.nComp();

  // Zero-gradient ghosts on outflow/vacuum faces: a ghost takes the value of
  // the boundary cell/node along the face normal, as the foextrap fill of
  // apply_field_bc() does. The ghosts of the other faces are left alone: the
  // caller zeroes the correction on inflow/fixed faces, conducting faces are
  // handled below, and absorbing and wave faces are not treated here (the
  // absorbing weights are unstable above c*dt/dx = 1 anyway).
  if (bcField.has(FieldBC::outflow) || bcField.has(FieldBC::vacuum)) {
    const BoxArray ba =
        (iLev > 0) ? get_level_active_ba(active_region_ref(), mf, Geom(iLev),
                                         Geom(0), nDim, iz_)
                   : get_boundary_active_ba(active_region_ref(), mf, Geom(iLev),
                                            nDim, iz_);
    const BoundaryBounds bnd(Geom(iLev), mf.boxArray().ixType(), &bcField);
    auto isOpen = [] AMREX_GPU_HOST_DEVICE(int type) {
      return type == FieldBC::outflow || type == FieldBC::vacuum;
    };

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bxFab = mfi.fabbox();
      if (ba.contains(bxFab))
        continue;

      Array4<Real> const& arr = mf[mfi].array();
      const Array4<const int>& statusArr = status[mfi].array();

      parallel_for_omp(bxFab, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (!bit::is_lev_boundary(statusArr(i, j, k, 0)))
          return;

        const int ijk[3] = { i, j, k };
        int m[3] = { i, j, k };
        bool touched = false;
        for (int d = 0; d < nDim; ++d) {
          if (isOpen(bnd.bcLo[d]) && ijk[d] < bnd.loBnd[d]) {
            m[d] = bnd.loBnd[d];
            touched = true;
          } else if (isOpen(bnd.bcHi[d]) && ijk[d] > bnd.hiBnd[d]) {
            m[d] = bnd.hiBnd[d];
            touched = true;
          }
        }
        if (!touched)
          return;

        for (int iVar = 0; iVar < nComp; ++iVar)
          arr(i, j, k, iVar) = arr(m[0], m[1], m[2], iVar);
      });
    }
  }

  if (hasConductingBC_)
    apply_conducting_wall(status, mf, 0, nComp, iLev, bcField, false);
}

//==========================================================
void Pic::fill_ext_dir(const iMultiFab& status, MultiFab& mf, const int iStart,
                       const int nComp, GETVALUE func, const int iLev,
                       const BoxBC<FieldBC::Type>* bc,
                       bool skipCoarseFineGhosts, Real timeOverride) {
  if (!func)
    return;

  // With skipCoarseFineGhosts only the ghosts outside the active region of the
  // level are filled. Otherwise every level-boundary ghost of a fab that is
  // not inside the level-0 active region is, which on a refined level includes
  // the coarse-fine interface ghosts.
  const BoxArray ba = skipCoarseFineGhosts
                          ? get_level_active_ba(active_region_ref(), mf,
                                                Geom(iLev), Geom(0), nDim, iz_)
                          : get_boundary_active_ba(active_region_ref(), mf,
                                                   Geom(iLev), nDim, iz_);
  const BoundaryBounds bnd(Geom(iLev), mf.boxArray().ixType(), bc);
  const Geometry& geom = Geom(iLev);
  const Real* probLo = geom.ProbLo();
  const Real* dx = geom.CellSize();
  const Real tEval =
      (timeOverride >= 0.0) ? timeOverride : (tc ? tc->get_time() : 0.0);
  const GETVALUE fnNodeB = &Pic::get_node_B;
  const GETVALUE fnCenterB = &Pic::get_center_B;
  const GETVALUE fnNodeE = &Pic::get_node_E;
  const GETVALUE fnCenterE = &Pic::get_center_E;
  const bool isDiscontinuityB = discontinuityBC.active && bc != nullptr &&
                                nComp == 3 &&
                                (func == fnNodeB || func == fnCenterB);
  const bool isDiscontinuityE = discontinuityBC.active && bc != nullptr &&
                                nComp == 3 &&
                                (func == fnNodeE || func == fnCenterE);
  const bool hasBC = (bc != nullptr);
  int periodic[3] = { 0, 0, 0 }, domainLen[3] = { 0, 0, 0 };
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    periodic[d] = geom.isPeriodic(d);
    domainLen[d] = geom.Domain().length(d);
  }
  std::unique_ptr<BoxListStorage> activeBoxes;
  if (skipCoarseFineGhosts)
    activeBoxes = std::make_unique<BoxListStorage>(ba);
  const BoxListView activeView =
      skipCoarseFineGhosts ? activeBoxes->view() : BoxListView();

  HostFabCache statusCache;
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    if (ba.contains(bxFab))
      continue;

    const Array4<const int> statusArr = statusCache.array(status, mfi);

    parallel_host_fill_fab(
        mf, mfi, bxFab, [&](int i, int j, int k, const Array4<Real>& arr) {
          if (!bit::is_lev_boundary(statusArr(i, j, k, 0)))
            return;

          if (skipCoarseFineGhosts &&
              !is_outside_active(activeView, periodic, domainLen, bnd, nDim,
                                 IntVect(AMREX_D_DECL(i, j, k))))
            return;

          // If at an outer physical boundary with non-Dirichlet condition,
          // PhysBCFunct and dedicated wall operators handle it.
          if (hasBC) {
            const int ijk[3] = { i, j, k };
            bool skipForPhysWall = false;
            bool onInflowGhost = false;
            for (int d = 0; d < nDim; ++d) {
              const auto isNonDirichletWall = [](int t) {
                return t == FieldBC::outflow || t == FieldBC::vacuum ||
                       t == FieldBC::conducting || t == FieldBC::absorb;
              };
              if ((ijk[d] < bnd.loBnd[d] && isNonDirichletWall(bnd.bcLo[d])) ||
                  (ijk[d] > bnd.hiBnd[d] && isNonDirichletWall(bnd.bcHi[d]))) {
                skipForPhysWall = true;
                break;
              }
              if ((ijk[d] < bnd.loBnd[d] && bnd.bcLo[d] == FieldBC::inflow) ||
                  (ijk[d] > bnd.hiBnd[d] && bnd.bcHi[d] == FieldBC::inflow)) {
                onInflowGhost = true;
              }
            }
            if (skipForPhysWall)
              return;

            if (onInflowGhost && (isDiscontinuityB || isDiscontinuityE)) {
              Real pos[3] = { 0.0, 0.0, 0.0 };
              for (int d = 0; d < nDim; ++d) {
                pos[d] = probLo[d] + (static_cast<Real>(ijk[d]) +
                                      (bnd.isNode[d] ? 0.0 : 0.5)) *
                                         dx[d];
              }
              Real val[3] = { 0.0, 0.0, 0.0 };
              if (isDiscontinuityB) {
                discontinuityBC.eval_B(pos, tEval, val);
              } else {
                discontinuityBC.eval_E(pos, tEval, val);
              }
              for (int iVar = 0; iVar < nComp; ++iVar) {
                arr(i, j, k, iStart + iVar) = val[iVar];
              }
              return;
            }
          }

          for (int iVar = 0; iVar < nComp; ++iVar) {
            arr(i, j, k, iStart + iVar) = (this->*func)(
                mfi, IntVect{ AMREX_D_DECL(i, j, k) }, iVar, iLev);
          }
        });
  }
}

//==========================================================
void Pic::apply_BC(const iMultiFab& status, MultiFab& mf, const int iStart,
                   const int nComp, GETVALUE func, const int iLev,
                   const BoxBC<FieldBC::Type>* bc, bool skipCoarseFineGhosts) {
  if (Geom(iLev).isAllPeriodic() || mf.nGrow() == 0)
    return;

  std::string nameFunc = "Pic::apply_BC";
  timing_func(nameFunc);

  if (bc != nullptr) {
    const Vector<BCRec> bcr = FieldBC::create_bcrec(
        *bc, FieldBC::Quantity::Scalar, nComp, Geom(iLev));
    GpuBndryFuncFab<FabFillNoOp> bfunc(FabFillNoOp{});
    PhysBCFunct<GpuBndryFuncFab<FabFillNoOp> > physbcf(Geom(iLev), bcr, bfunc);
    physbcf(mf, iStart, nComp, mf.nGrowVect(), 0.0, 0);

    if (func != nullptr) {
      fill_ext_dir(status, mf, iStart, nComp, func, iLev, bc,
                   skipCoarseFineGhosts);
    }
  } else if (func == nullptr) {
    // Float / extrapolation BC on physical and embedded level boundaries:
    // 1. Native AMReX foextrap at physical domain boundaries:
    Vector<BCRec> bcr(nComp);
    for (int c = 0; c < nComp; ++c) {
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (Geom(iLev).isPeriodic(d)) {
          bcr[c].setLo(d, BCType::int_dir);
          bcr[c].setHi(d, BCType::int_dir);
        } else {
          bcr[c].setLo(d, BCType::foextrap);
          bcr[c].setHi(d, BCType::foextrap);
        }
      }
    }
    GpuBndryFuncFab<FabFillNoOp> bfunc(FabFillNoOp{});
    PhysBCFunct<GpuBndryFuncFab<FabFillNoOp> > physbcf(Geom(iLev), bcr, bfunc);
    physbcf(mf, iStart, nComp, mf.nGrowVect(), 0.0, 0);

    // 2. Extrapolate from nearest valid neighbor for embedded active-region
    // boundaries. When skipCoarseFineGhosts is true, coarse-fine interface
    // ghosts (which lie inside the active region) are preserved.
    const BoxArray ba =
        skipCoarseFineGhosts
            ? get_level_active_ba(active_region_ref(), mf, Geom(iLev), Geom(0),
                                  nDim, iz_)
            : get_boundary_active_ba(active_region_ref(), mf, Geom(iLev), nDim,
                                     iz_);
    const BoundaryBounds bnd(Geom(iLev), mf.boxArray().ixType(), bc);
    const Geometry& geom = Geom(iLev);
    int periodic[3] = { 0, 0, 0 }, domainLen[3] = { 0, 0, 0 };
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      periodic[d] = geom.isPeriodic(d);
      domainLen[d] = geom.Domain().length(d);
    }
    std::unique_ptr<BoxListStorage> activeBoxes;
    if (skipCoarseFineGhosts)
      activeBoxes = std::make_unique<BoxListStorage>(ba);
    const BoxListView activeView =
        skipCoarseFineGhosts ? activeBoxes->view() : BoxListView();
    const int nDimLoc = nDim;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bxFab = mfi.fabbox();
      const Box& bxValid = mfi.validbox();
      if (skipCoarseFineGhosts && ba.contains(bxFab))
        continue;
      if (!skipCoarseFineGhosts && iLev == 0 && ba.contains(bxFab))
        continue;

      Array4<Real> const arr = mf[mfi].array();
      const Array4<const int> statusArr = status[mfi].const_array();
      Box box = bxValid;
      box.grow(1);
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (bit::is_lev_boundary(statusArr(i, j, k, 0))) {
          if (skipCoarseFineGhosts &&
              !is_outside_active(activeView, periodic, domainLen, bnd, nDimLoc,
                                 IntVect(AMREX_D_DECL(i, j, k)))) {
            return;
          }
          bool isNeiFound = false;
          const int kmin = (nDimLoc > 2) ? -1 : 0;
          const int kmax = (nDimLoc > 2) ? 1 : 0;
          for (int kk = kmin; kk <= kmax && !isNeiFound; ++kk) {
            for (int jj = -1; jj <= 1 && !isNeiFound; ++jj) {
              for (int ii = -1; ii <= 1 && !isNeiFound; ++ii) {
                if (!bit::is_lev_boundary(
                        statusArr(i + ii, j + jj, k + kk, 0))) {
                  isNeiFound = true;
                  for (int iVar = iStart; iVar < iStart + nComp; ++iVar) {
                    arr(i, j, k, iVar) = arr(i + ii, j + jj, k + kk, iVar);
                  }
                }
              }
            }
          }
        }
      });
    }
  } else {
    fill_ext_dir(status, mf, iStart, nComp, func, iLev, nullptr,
                 skipCoarseFineGhosts);
  }
}

//==========================================================
void Pic::apply_conducting_wall(const iMultiFab& status, MultiFab& mf,
                                const int iStart, const int nComp,
                                const int iLev, const BoxBC<FieldBC::Type>& bc,
                                bool isB) {
  std::string nameFunc = "Pic::apply_conducting_wall";
  timing_func(nameFunc);

  const BoxArray ba =
      get_boundary_active_ba(active_region_ref(), mf, Geom(iLev), nDim, iz_);
  const BoundaryBounds bnd(Geom(iLev), mf.boxArray().ixType(), &bc);

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    if (ba.contains(bxFab))
      continue;

    Array4<Real> const& arr = mf[mfi].array();
    const Array4<const int>& statusArr = status[mfi].array();
    const Box& bxValid = mfi.validbox();
    const Dim3 vLo = bxValid.smallEnd().dim3();
    const Dim3 vHi = bxValid.bigEnd().dim3();

    parallel_for_omp(bxFab, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const int ijk[3] = { i, j, k };
      const int vLoArr[3] = { vLo.x, vLo.y, vLo.z };
      const int vHiArr[3] = { vHi.x, vHi.y, vHi.z };

      // 1. Boundary nodes on physical wall (node-centered only).
      for (int d = 0; d < nDim; ++d) {
        if (!bnd.isNode[d])
          continue;

        const bool onLoWall =
            (bnd.bcLo[d] == FieldBC::conducting) && (ijk[d] == bnd.loBnd[d]);
        const bool onHiWall =
            (bnd.bcHi[d] == FieldBC::conducting) && (ijk[d] == bnd.hiBnd[d]);
        if (onLoWall || onHiWall) {
          bool inValid = true;
          for (int od = 0; od < nDim; ++od) {
            if (od != d && (ijk[od] < vLoArr[od] || ijk[od] > vHiArr[od])) {
              inValid = false;
              break;
            }
          }
          if (inValid) {
            for (int iVar = 0; iVar < nComp; ++iVar) {
              const int comp = iStart + iVar;
              if (isB) {
                if (iVar == d)
                  arr(i, j, k, comp) = 0.0;
              } else {
                if (iVar != d)
                  arr(i, j, k, comp) = 0.0;
              }
            }
          }
        }
      }

      // 2. Ghost cells/nodes outside domain.
      if (!bit::is_lev_boundary(statusArr(i, j, k, 0)))
        return;

      int m[3] = { i, j, k };
      bool isCondLow[3] = { false, false, false };
      bool isCondHigh[3] = { false, false, false };
      bool touched = false;

      for (int d = 0; d < nDim; ++d) {
        if (bnd.bcLo[d] == FieldBC::conducting && ijk[d] < bnd.loBnd[d]) {
          isCondLow[d] = true;
          m[d] = bnd.isNode[d] ? (2 * bnd.loBnd[d] - ijk[d])
                               : (2 * bnd.loBnd[d] - 1 - ijk[d]);
          touched = true;
        } else if (bnd.bcHi[d] == FieldBC::conducting &&
                   ijk[d] > bnd.hiBnd[d]) {
          isCondHigh[d] = true;
          m[d] = bnd.isNode[d] ? (2 * bnd.hiBnd[d] - ijk[d])
                               : (2 * bnd.hiBnd[d] + 1 - ijk[d]);
          touched = true;
        }
      }

      if (!touched)
        return;

      for (int iVar = 0; iVar < nComp; ++iVar) {
        const int comp = iStart + iVar;
        if (isB) {
          bool isNormal = false;
          for (int d = 0; d < nDim; ++d) {
            if ((isCondLow[d] || isCondHigh[d]) && iVar == d) {
              isNormal = true;
              break;
            }
          }
          if (isNormal) {
            arr(i, j, k, comp) = -arr(m[0], m[1], m[2], comp);
          } else {
            arr(i, j, k, comp) = arr(m[0], m[1], m[2], comp);
          }
        } else {
          bool isTangential = false;
          for (int d = 0; d < nDim; ++d) {
            if ((isCondLow[d] || isCondHigh[d]) && iVar != d) {
              isTangential = true;
              break;
            }
          }
          if (isTangential) {
            arr(i, j, k, comp) = -arr(m[0], m[1], m[2], comp);
          } else {
            arr(i, j, k, comp) = arr(m[0], m[1], m[2], comp);
          }
        }
      }
    });
  }
}

//==========================================================
void Pic::apply_absorbing_wall(const iMultiFab& status, MultiFab& mf,
                               const int iStart, const int nComp,
                               const int iLev, const BoxBC<FieldBC::Type>& bc,
                               bool isB, Real dtOverride) {
  // dtOverride < -1.5 signals a non-advancing refresh: keep stateful ghosts.
  if (dtOverride < -1.5)
    return;

  const Real dt = (dtOverride > 0.0) ? dtOverride : (tc ? tc->get_dt() : 0.0);
  const bool zeroGradInit = (dtOverride == 0.0) || (dt <= 0.0);

  std::string nameFunc = "Pic::apply_absorbing_wall";
  timing_func(nameFunc);

  const BoxArray ba =
      get_boundary_active_ba(active_region_ref(), mf, Geom(iLev), nDim, iz_);
  const BoundaryBounds bnd(Geom(iLev), mf.boxArray().ixType(), &bc);
  const Real* dx = Geom(iLev).CellSize();

  if (zeroGradInit) {
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bxFab = mfi.fabbox();
      if (ba.contains(bxFab))
        continue;

      Array4<Real> const& arr = mf[mfi].array();
      const Array4<const int>& statusArr = status[mfi].array();

      parallel_for_omp(bxFab, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (!bit::is_lev_boundary(statusArr(i, j, k, 0)))
          return;

        const int ijk[3] = { i, j, k };
        int m[3] = { i, j, k };
        bool touched = false;

        for (int d = 0; d < nDim; ++d) {
          if (bnd.bcLo[d] == FieldBC::absorb && ijk[d] < bnd.loBnd[d]) {
            m[d] = bnd.loBnd[d];
            touched = true;
          } else if (bnd.bcHi[d] == FieldBC::absorb && ijk[d] > bnd.hiBnd[d]) {
            m[d] = bnd.hiBnd[d];
            touched = true;
          }
        }

        if (!touched)
          return;

        for (int iVar = 0; iVar < nComp; ++iVar) {
          const int comp = iStart + iVar;
          arr(i, j, k, comp) = arr(m[0], m[1], m[2], comp);
        }
      });
    }
    return;
  }

  // Characteristic speed: explicit #ABSORB value when given; otherwise c=1 in
  // full-PIC, or local fast magnetosonic speed in hybrid-PIC (defaulting to 1
  // in vacuum).
  const Real csFixed = (absorbCharSpeed > 0.0) ? absorbCharSpeed : 1.0;
  const bool autoHybridSpeed = (useHybridPIC && absorbCharSpeed <= 0.0);
  const bool hasHybridSources =
      autoHybridSpeed && !nodePlasma.empty() &&
      nSpecies < static_cast<int>(nodePlasma.size()) &&
      !nodePlasma[nSpecies][iLev].empty() && !centerB[iLev].empty() &&
      (mf.boxArray().size() == boxArray(iLev).size());
  const bool hasB0 =
      hasHybridSources && use_intrinsic_B() && !centerB0[iLev].empty();
  const MultiFab& hybridBMF =
      (isB && mf.boxArray().ixType().cellCentered() && nComp >= 3)
          ? mf
          : centerB[iLev];
  const Real cs2 =
      (electronTemperature > 0.0)
          ? ((electronGamma > 0.0) ? electronGamma : 1.0) * electronTemperature
          : 0.0;
  const Real rhoFloor = amrex::max(rhoMinOhm, static_cast<Real>(1e-10));
  const Real invFourPI = 1.0 / fourPI;
  const Box domBox = Geom(iLev).Domain();
  const int domHi[3] = { domBox.bigEnd(0), domBox.bigEnd(1),
                         (nDim > 2) ? domBox.bigEnd(2) : 0 };
  const int domLo[3] = { domBox.smallEnd(0), domBox.smallEnd(1),
                         (nDim > 2) ? domBox.smallEnd(2) : 0 };

  AbsorbWeights weights;
  for (int d = 0; d < nDim; ++d) {
    const Real drive0 = csFixed * dt / dx[d];
    weights.decay[d] = (1.0 - drive0) / (1.0 + drive0);
    weights.drive[d] = 2.0 * drive0 / (1.0 + drive0);
  }
  const Real invDx[3] = { 1.0 / dx[0], (nDim > 1) ? (1.0 / dx[1]) : 0.0,
                          (nDim > 2) ? (1.0 / dx[2]) : 0.0 };
  const int rhoComp = iRho_;

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    if (ba.contains(bxFab))
      continue;

    Array4<Real> const& arr = mf[mfi].array();
    const Array4<const int>& statusArr = status[mfi].array();
    const Array4<const Real> bArr =
        hasHybridSources ? hybridBMF[mfi].const_array() : Array4<const Real>();
    const Array4<const Real> b0Arr =
        hasB0 ? centerB0[iLev][mfi].const_array() : Array4<const Real>();
    const Array4<const Real> rhoArr =
        hasHybridSources ? nodePlasma[nSpecies][iLev][mfi].const_array()
                         : Array4<const Real>();

    parallel_for_omp(bxFab, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      if (!bit::is_lev_boundary(statusArr(i, j, k, 0)))
        return;

      const int ijk[3] = { i, j, k };
      int m[3] = { i, j, k };
      bool isAbsorb[3] = { false, false, false };
      int nAbsorb = 0;

      for (int d = 0; d < nDim; ++d) {
        const bool isLow =
            (bnd.bcLo[d] == FieldBC::absorb) && (ijk[d] < bnd.loBnd[d]);
        const bool isHigh =
            (bnd.bcHi[d] == FieldBC::absorb) && (ijk[d] > bnd.hiBnd[d]);
        if (isLow) {
          m[d] = bnd.isNode[d] ? (2 * bnd.loBnd[d] - ijk[d])
                               : (2 * bnd.loBnd[d] - 1 - ijk[d]);
          isAbsorb[d] = true;
          nAbsorb++;
        } else if (isHigh) {
          m[d] = bnd.isNode[d] ? (2 * bnd.hiBnd[d] - ijk[d])
                               : (2 * bnd.hiBnd[d] + 1 - ijk[d]);
          isAbsorb[d] = true;
          nAbsorb++;
        }
      }

      if (nAbsorb == 0)
        return;

      Real cellDecay = 1.0;
      Real cellDrive = 0.0;

      if (hasHybridSources) {
        const int cm0 = std::clamp(m[0], domLo[0], domHi[0]);
        const int cm1 = std::clamp(m[1], domLo[1], domHi[1]);
        const int cm2 = (nDim > 2) ? std::clamp(m[2], domLo[2], domHi[2]) : 0;
        Real csLocal = 1.0;
        if (bArr.contains(cm0, cm1, cm2) && rhoArr.contains(cm0, cm1, cm2)) {
          const Real rhoCell = rhoArr(cm0, cm1, cm2, rhoComp);
          if (rhoCell > 0.0) {
            Real bx = bArr(cm0, cm1, cm2, ix_);
            Real by = bArr(cm0, cm1, cm2, iy_);
            Real bz = bArr(cm0, cm1, cm2, iz_);
            if (hasB0 && b0Arr.contains(cm0, cm1, cm2)) {
              bx += b0Arr(cm0, cm1, cm2, ix_);
              by += b0Arr(cm0, cm1, cm2, iy_);
              bz += b0Arr(cm0, cm1, cm2, iz_);
            }
            const Real b2 = bx * bx + by * by + bz * bz;
            const Real rhoEff = amrex::max(rhoCell, rhoFloor);
            csLocal = std::sqrt(b2 * invFourPI / rhoEff + cs2);
          }
        }
        for (int d = 0; d < nDim; ++d) {
          if (!isAbsorb[d])
            continue;
          const Real drive0 =
              amrex::min(csLocal * dt * invDx[d], static_cast<Real>(1.0));
          cellDecay *= (1.0 - drive0) / (1.0 + drive0);
          cellDrive += 2.0 * drive0 / (1.0 + drive0);
        }
      } else {
        for (int d = 0; d < nDim; ++d) {
          if (!isAbsorb[d])
            continue;
          cellDecay *= weights.decay[d];
          cellDrive += weights.drive[d];
        }
      }

      cellDrive /= static_cast<Real>(nAbsorb);

      for (int iVar = 0; iVar < nComp; ++iVar) {
        const int comp = iStart + iVar;
        arr(i, j, k, comp) = cellDecay * arr(i, j, k, comp) +
                             cellDrive * arr(m[0], m[1], m[2], comp);
      }
    });
  }
}

//==========================================================
void Pic::apply_inflow_wall(const iMultiFab& status, MultiFab& mf,
                            const int iStart, const int nComp, const int iLev,
                            const BoxBC<FieldBC::Type>& bc, bool isB,
                            GETVALUE func, Real timeOverride) {
  // isB indicates whether mf is magnetic (true) or electric (false).

  std::string nameFunc = "Pic::apply_inflow_wall";
  timing_func(nameFunc);

  const BoxArray ba =
      get_boundary_active_ba(active_region_ref(), mf, Geom(iLev), nDim, iz_);
  const BoundaryBounds bnd(Geom(iLev), mf.boxArray().ixType(), &bc);
  const Geometry& geom = Geom(iLev);
  const Real* probLo = geom.ProbLo();
  const Real* dx = geom.CellSize();
  const Real tEval =
      (timeOverride >= 0.0) ? timeOverride : (tc ? tc->get_time() : 0.0);
  const GETVALUE fnNodeB = &Pic::get_node_B;
  const GETVALUE fnCenterB = &Pic::get_center_B;
  const GETVALUE fnNodeE = &Pic::get_node_E;
  const GETVALUE fnCenterE = &Pic::get_center_E;
  const bool useDiscOnWall = discontinuityBC.active && nComp == 3 &&
                             func != nullptr &&
                             (func == fnNodeB || func == fnCenterB ||
                              func == fnNodeE || func == fnCenterE);

  const bool hasFunc = (func != nullptr);
  const int nDimLoc = nDim;

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    if (ba.contains(bxFab))
      continue;

    const Box& bxValid = mfi.validbox();
    const Dim3 vLo = bxValid.smallEnd().dim3();
    const Dim3 vHi = bxValid.bigEnd().dim3();
    const int vLoArr[3] = { vLo.x, vLo.y, vLo.z };
    const int vHiArr[3] = { vHi.x, vHi.y, vHi.z };

    // 1. Boundary nodes on physical inflow wall (node-centered only).
    if (hasFunc || useDiscOnWall) {
      parallel_host_fill_fab(
          mf, mfi, bxFab, [&](int i, int j, int k, const Array4<Real>& arr) {
            const int ijk[3] = { i, j, k };
            for (int d = 0; d < nDimLoc; ++d) {
              if (!bnd.isNode[d])
                continue;

              const bool onLoWall = (FieldBC::is_dirichlet(bnd.bcLo[d]) &&
                                     (ijk[d] == bnd.loBnd[d]));
              const bool onHiWall = (FieldBC::is_dirichlet(bnd.bcHi[d]) &&
                                     (ijk[d] == bnd.hiBnd[d]));
              if (onLoWall || onHiWall) {
                bool inValid = true;
                for (int od = 0; od < nDimLoc; ++od) {
                  if (od != d &&
                      (ijk[od] < vLoArr[od] || ijk[od] > vHiArr[od])) {
                    inValid = false;
                    break;
                  }
                }
                if (inValid) {
                  const bool isInflowFace =
                      (onLoWall && bnd.bcLo[d] == FieldBC::inflow) ||
                      (onHiWall && bnd.bcHi[d] == FieldBC::inflow);
                  if (isInflowFace && useDiscOnWall) {
                    Real pos[3] = { 0.0, 0.0, 0.0 };
                    for (int od = 0; od < nDimLoc; ++od) {
                      pos[od] = probLo[od] + (static_cast<Real>(ijk[od]) +
                                              (bnd.isNode[od] ? 0.0 : 0.5)) *
                                                 dx[od];
                    }
                    Real val[3] = { 0.0, 0.0, 0.0 };
                    if (isB) {
                      discontinuityBC.eval_B(pos, tEval, val);
                    } else {
                      discontinuityBC.eval_E(pos, tEval, val);
                    }
                    for (int iVar = 0; iVar < nComp; ++iVar) {
                      arr(i, j, k, iStart + iVar) = val[iVar];
                    }
                  } else {
                    for (int iVar = 0; iVar < nComp; ++iVar) {
                      const int comp = iStart + iVar;
                      if (hasFunc) {
                        arr(i, j, k, comp) = (this->*func)(
                            mfi, IntVect{ AMREX_D_DECL(i, j, k) }, iVar, iLev);
                      } else {
                        int m[3] = { i, j, k };
                        m[d] =
                            onLoWall ? (bnd.loBnd[d] + 1) : (bnd.hiBnd[d] - 1);
                        arr(i, j, k, comp) = arr(m[0], m[1], m[2], comp);
                      }
                    }
                  }
                }
              }
            }
          });
    } else {
      Array4<Real> const arr = mf[mfi].array();
      parallel_for_omp(bxFab, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const int ijk[3] = { i, j, k };
        for (int d = 0; d < nDimLoc; ++d) {
          if (!bnd.isNode[d])
            continue;

          const bool onLoWall =
              (FieldBC::is_dirichlet(bnd.bcLo[d]) && (ijk[d] == bnd.loBnd[d]));
          const bool onHiWall =
              (FieldBC::is_dirichlet(bnd.bcHi[d]) && (ijk[d] == bnd.hiBnd[d]));
          if (onLoWall || onHiWall) {
            bool inValid = true;
            for (int od = 0; od < nDimLoc; ++od) {
              if (od != d && (ijk[od] < vLoArr[od] || ijk[od] > vHiArr[od])) {
                inValid = false;
                break;
              }
            }
            if (inValid) {
              for (int iVar = 0; iVar < nComp; ++iVar) {
                const int comp = iStart + iVar;
                int m[3] = { i, j, k };
                m[d] = onLoWall ? (bnd.loBnd[d] + 1) : (bnd.hiBnd[d] - 1);
                arr(i, j, k, comp) = arr(m[0], m[1], m[2], comp);
              }
            }
          }
        }
      });
    }

    // Cell-centered data is not pinned here: for cell-centered B the cell at
    // loBnd/hiBnd is the last *physical* cell, advanced by Faraday's law.
    // Overwriting it with the prescribed state every call breaks div(B) = 0
    // against its evolving neighbour, so Bx grows in a layer next to the
    // face. The prescribed state enters through the ghost cells (filled by
    // fill_ext_dir) and the pinned boundary nodes above.

    // 2. Ghost cells/nodes outside domain if func == nullptr
    // (When func != nullptr, fill_ext_dir handles them). This is a second
    // pass, so every ghost copies the updated wall node regardless of the
    // loop (thread) order.
    if (hasFunc)
      continue;

    Array4<Real> const arr = mf[mfi].array();
    const Array4<const int> statusArr = status[mfi].const_array();
    parallel_for_omp(bxFab, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      if (!bit::is_lev_boundary(statusArr(i, j, k, 0)))
        return;
      const int ijk[3] = { i, j, k };
      int m[3] = { i, j, k };
      bool touched = false;
      for (int d = 0; d < nDimLoc; ++d) {
        if (FieldBC::is_dirichlet(bnd.bcLo[d]) && ijk[d] < bnd.loBnd[d]) {
          m[d] = bnd.loBnd[d];
          touched = true;
        } else if (FieldBC::is_dirichlet(bnd.bcHi[d]) &&
                   ijk[d] > bnd.hiBnd[d]) {
          m[d] = bnd.hiBnd[d];
          touched = true;
        }
      }
      if (touched) {
        for (int iVar = 0; iVar < nComp; ++iVar) {
          arr(i, j, k, iStart + iVar) = arr(m[0], m[1], m[2], iStart + iVar);
        }
      }
    });
  }
}

//==========================================================
void Pic::ensure_wave_bc_initialized(bool allowUniformFallback) {
  if (!waveBC.turb.active || waveBC.turb.initialized)
    return;
  bool isWaveFace[6] = { false, false, false, false, false, false };
  bool anyWaveFace = false;
  for (int d = 0; d < nDim; ++d) {
    isWaveFace[2 * d + 0] = (bcField.face(d, 0) == FieldBC::wave);
    isWaveFace[2 * d + 1] = (bcField.face(d, 1) == FieldBC::wave);
    if (isWaveFace[2 * d + 0] || isWaveFace[2 * d + 1])
      anyWaveFace = true;
  }
  if (!anyWaveFace)
    return;
  if (!fi)
    return;
  const bool hasNodeFluid =
      fi->is_node_fluid_ready() && !fi->get_node_fluid(0).empty();
  const bool hasUniformFallback =
      allowUniformFallback && !fi->is_init_from_swmf() && (fi->get_nS() > 0) &&
      (static_cast<int>(fi->get_uniform_state().size()) >=
       5 * fi->get_nS() + 3) &&
      (static_cast<int>(fi->inflowState.size()) >= fi->get_nS());
  if (!hasNodeFluid && !hasUniformFallback)
    return;
  waveBC.ensure_turbulence_initialized(fi, Geom(0), nDim, is_fake_2d(),
                                       isWaveFace);
}

//==========================================================
void Pic::apply_wave_field(const iMultiFab& status, MultiFab& mf,
                           const int iStart, const int nComp, const int iLev,
                           const BoxBC<FieldBC::Type>& bc, int iField, Real t,
                           GETVALUE func) {
  bool isWaveFace[6] = { false, false, false, false, false, false };
  bool anyWaveFace = false;
  for (int d = 0; d < nDim; ++d) {
    isWaveFace[2 * d + 0] = (bc.face(d, 0) == FieldBC::wave);
    isWaveFace[2 * d + 1] = (bc.face(d, 1) == FieldBC::wave);
    if (isWaveFace[2 * d + 0] || isWaveFace[2 * d + 1])
      anyWaveFace = true;
  }
  if (!anyWaveFace)
    return;

  if (waveBC.turb.active && !waveBC.turb.initialized) {
    waveBC.ensure_turbulence_initialized(fi, Geom(0), nDim, is_fake_2d(),
                                         isWaveFace);
  }

  const Real turbEnv = waveBC.turbulent_envelope(t);
  const bool hasTurb = waveBC.turb.active && waveBC.turb.initialized &&
                       (turbEnv > 0.0) && (iField == 0 || iField == 1);
  std::vector<Real> turbPhaseT;
  if (hasTurb) {
    const int N = waveBC.turb.nMode;
    turbPhaseT.resize(N);
    for (int m = 0; m < N; ++m)
      turbPhaseT[m] = waveBC.turb.phase[m] - waveBC.turb.omega[m] * t;
  }
  const Real* turbPhaseTPtr = hasTurb ? turbPhaseT.data() : nullptr;

  // Pre-filter legacy #WAVEBC components active at time t for iField.
  // Hoisting envelope and time-phase evaluation out of the cell loops.
  // (ActiveComp is defined at namespace scope above.)
  std::vector<ActiveComp> activeComps;
  for (const auto& f : waveBC.faces) {
    if (f.direction >= nDim || bc.face(f.direction, f.side) != FieldBC::wave)
      continue;
    for (const auto& c : f.comps) {
      if (c.iField != iField || c.amplitude == 0.0)
        continue;
      const Real env = waveBC.envelope(c, t);
      if (env <= 0.0)
        continue;
      ActiveComp ac;
      ac.dir = f.direction;
      ac.side = f.side;
      ac.effAmp = c.amplitude * env;
      ac.omega_t = c.frequency * t - c.phase;
      for (int d = 0; d < 3; ++d) {
        ac.k_vec[d] = c.k_vec[d];
        ac.pol[d] = c.pol[d];
      }
      ac.profile = c.profile;
      ac.tCenter = c.tCenter;
      ac.tInvWidth = (c.tWidth > 0.0) ? (1.0 / c.tWidth) : 1.0;
      ac.orig = &c;
      activeComps.push_back(ac);
    }
  }
  if (!hasTurb && activeComps.empty())
    return;

  std::string nameFunc = "Pic::apply_wave_field";
  timing_func(nameFunc);

  const BoxArray ba =
      get_boundary_active_ba(active_region_ref(), mf, Geom(iLev), nDim, iz_);
  const BoundaryBounds bnd(Geom(iLev), mf.boxArray().ixType(), &bc);

  const Real* plo = Geom(iLev).ProbLo();
  const Real* dx = Geom(iLev).CellSize();
  const Real offset[3] = { bnd.isNode[0] ? 0.0 : 0.5, bnd.isNode[1] ? 0.0 : 0.5,
                           bnd.isNode[2] ? 0.0 : 0.5 };

  HostFabCache statusCache;
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    const Box& bxValid = mfi.validbox();
    if (ba.contains(bxFab))
      continue;

    // Box culling: check if this FAB touches any face active for wave
    // injection.
    bool boxTouchesWaveFace = false;
    if (hasTurb) {
      for (int d = 0; d < nDim && !boxTouchesWaveFace; ++d) {
        if (isWaveFace[2 * d + 0] &&
            (bxFab.smallEnd(d) < bnd.loBnd[d] ||
             (bnd.isNode[d] && bxFab.smallEnd(d) <= bnd.loBnd[d]))) {
          boxTouchesWaveFace = true;
        }
        if (isWaveFace[2 * d + 1] &&
            (bxFab.bigEnd(d) > bnd.hiBnd[d] ||
             (bnd.isNode[d] && bxFab.bigEnd(d) >= bnd.hiBnd[d]))) {
          boxTouchesWaveFace = true;
        }
      }
    }
    for (const auto& ac : activeComps) {
      if (boxTouchesWaveFace)
        break;
      const int d = ac.dir;
      const int side = ac.side;
      if (side == 0) {
        if (bxFab.smallEnd(d) < bnd.loBnd[d] ||
            (bnd.isNode[d] && bxFab.smallEnd(d) <= bnd.loBnd[d])) {
          boxTouchesWaveFace = true;
          break;
        }
      } else {
        if (bxFab.bigEnd(d) > bnd.hiBnd[d] ||
            (bnd.isNode[d] && bxFab.bigEnd(d) >= bnd.hiBnd[d])) {
          boxTouchesWaveFace = true;
          break;
        }
      }
    }
    if (!boxTouchesWaveFace)
      continue;

    const Array4<const int> statusArr = statusCache.array(status, mfi);
    const Dim3 vLo = bxValid.smallEnd().dim3();
    const Dim3 vHi = bxValid.bigEnd().dim3();
    const int vLoArr[3] = { vLo.x, vLo.y, vLo.z };
    const int vHiArr[3] = { vHi.x, vHi.y, vHi.z };

    parallel_host_fill_fab(
        mf, mfi, bxFab, [&](int i, int j, int k, const Array4<Real>& arr) {
          const int ijk[3] = { i, j, k };
          const bool isGhost = bit::is_lev_boundary(statusArr(i, j, k, 0));

          Real pos[3] = { plo[0] + dx[0] * (i + offset[0]),
                          plo[1] + dx[1] * (j + offset[1]),
                          (nDim > 2) ? plo[2] + dx[2] * (k + offset[2]) : 0.0 };

          Real waveVal[3] = { 0.0, 0.0, 0.0 };
          bool hasWave = false;
          bool isBndNode = false;
          int matchedDir = -1;
          int matchedSide = -1;

          for (const auto& ac : activeComps) {
            const int d = ac.dir;
            const int side = ac.side;
            if (matchedDir >= 0 && (d != matchedDir || side != matchedSide))
              continue;
            const int idx = ijk[d];
            const bool onGhostFace =
                isGhost && ((side == 0 && idx < bnd.loBnd[d]) ||
                            (side == 1 && idx > bnd.hiBnd[d]));

            // Boundary node on physical wall (node-centered only).
            bool onNodeWall = false;
            if (bnd.isNode[d]) {
              const bool onLoNode = (side == 0 && idx == bnd.loBnd[d]);
              const bool onHiNode = (side == 1 && idx == bnd.hiBnd[d]);
              if (onLoNode || onHiNode) {
                bool inTangential = true;
                for (int od = 0; od < nDim; ++od) {
                  if (od != d &&
                      (ijk[od] < vLoArr[od] || ijk[od] > vHiArr[od])) {
                    inTangential = false;
                    break;
                  }
                }
                if (inTangential)
                  onNodeWall = true;
              }
            }

            if (!onGhostFace && !onNodeWall)
              continue;

            if (matchedDir < 0) {
              matchedDir = d;
              matchedSide = side;
            }

            if (onNodeWall)
              isBndNode = true;

            // Evaluate wave profile value with precomputed parameters.
            Real val = 0.0;
            if (ac.profile == WaveComponent::kCustom && ac.orig &&
                ac.orig->custom) {
              val = ac.orig->custom(*ac.orig, t, pos) *
                    (ac.effAmp / ac.orig->amplitude);
            } else {
              const Real kdotx = ac.k_vec[0] * pos[0] + ac.k_vec[1] * pos[1] +
                                 ac.k_vec[2] * pos[2];
              val = ac.effAmp * std::sin(kdotx - ac.omega_t);
              if (ac.profile == WaveComponent::kPacket) {
                const Real tau = (t - ac.tCenter) * ac.tInvWidth;
                val *= std::exp(-tau * tau);
              }
            }

            if (val == 0.0)
              continue;

            hasWave = true;
            if (iField == 0 || iField == 1) {
              for (int iVar = 0; iVar < std::min(nComp, 3); ++iVar) {
                waveVal[iVar] += val * ac.pol[iVar];
              }
            } else {
              waveVal[0] += val;
            }
          }

          if (hasWave) {
            const Real w = waveBC.window_weight(pos);
            for (int iVar = 0; iVar < 3; ++iVar)
              waveVal[iVar] *= w;
          }

          if (hasTurb) {
            int turbFaceMask = 0;
            for (int d = 0; d < nDim; ++d) {
              for (int side = 0; side < 2; ++side) {
                const int fIdx = 2 * d + side;
                if (!isWaveFace[fIdx])
                  continue;
                const int idx = ijk[d];
                const bool onGhostFace =
                    isGhost && ((side == 0 && idx < bnd.loBnd[d]) ||
                                (side == 1 && idx > bnd.hiBnd[d]));
                bool onNodeWall = false;
                if (bnd.isNode[d]) {
                  const bool onLoNode = (side == 0 && idx == bnd.loBnd[d]);
                  const bool onHiNode = (side == 1 && idx == bnd.hiBnd[d]);
                  if (onLoNode || onHiNode) {
                    bool inTangential = true;
                    for (int od = 0; od < nDim; ++od) {
                      if (od != d &&
                          (ijk[od] < vLoArr[od] || ijk[od] > vHiArr[od])) {
                        inTangential = false;
                        break;
                      }
                    }
                    if (inTangential)
                      onNodeWall = true;
                  }
                }
                if (onGhostFace || onNodeWall) {
                  turbFaceMask |= (1 << fIdx);
                  if (onNodeWall)
                    isBndNode = true;
                }
              }
            }
            if (turbFaceMask != 0) {
              Real turbVal[3] = { 0.0, 0.0, 0.0 };
              waveBC.eval_turbulent_field(pos, turbPhaseTPtr, turbEnv, iField,
                                          turbVal, turbFaceMask);
              if (turbVal[0] != 0.0 || turbVal[1] != 0.0 || turbVal[2] != 0.0) {
                hasWave = true;
                for (int iVar = 0; iVar < 3; ++iVar)
                  waveVal[iVar] += turbVal[iVar];
              }
            }
          }

          if (hasWave) {
            if (isBndNode && func) {
              // Reset nodal physical boundary from base state before
              // superimposing wave.
              for (int iVar = 0; iVar < nComp; ++iVar) {
                const Real baseVal = (this->*func)(
                    mfi, IntVect{ AMREX_D_DECL(i, j, k) }, iVar, iLev);
                arr(i, j, k, iStart + iVar) = baseVal + waveVal[iVar % 3];
              }
            } else {
              // Ghost cells already contain base state from apply_BC(); add
              // wave perturbation.
              if (iField == 0 || iField == 1) {
                for (int iVar = 0; iVar < nComp; ++iVar) {
                  arr(i, j, k, iStart + iVar) += waveVal[iVar % 3];
                }
              } else if (nComp > 0) {
                arr(i, j, k, iStart) += waveVal[0];
              }
            }
          }
        });
  }
}

//==========================================================
void Pic::wave_velocity_kick(const Real* pos, Real t, Real& dvx, Real& dvy,
                             Real& dvz, int waveFaceMask, int iSpecies,
                             Real charge) {
  dvx = dvy = dvz = 0.0;
  if (!waveBC.active)
    return;

  if (!waveBC.faces.empty()) {
    int matchedDir = -1;
    int matchedSide = -1;
    for (const auto& f : waveBC.faces) {
      if (f.direction >= nDim ||
          bcField.face(f.direction, f.side) != FieldBC::wave)
        continue;
      const int faceBit = 1 << (2 * f.direction + f.side);
      if (waveFaceMask >= 0 && (waveFaceMask & faceBit) == 0)
        continue;
      if (matchedDir >= 0 &&
          (f.direction != matchedDir || f.side != matchedSide))
        continue;
      for (const auto& c : f.comps) {
        if (c.iField != 2 || c.amplitude == 0.0) // velocity kick
          continue;
        if (c.iSpecies >= 0 && iSpecies >= 0 && c.iSpecies != iSpecies)
          continue;
        if (waveBC.envelope(c, t) <= 0.0)
          continue;
        if (matchedDir < 0) {
          matchedDir = f.direction;
          matchedSide = f.side;
        }
        const Real val = waveBC.value(c, t, pos);
        if (val == 0.0)
          continue;
        dvx += val * c.pol[0];
        dvy += val * c.pol[1];
        dvz += val * c.pol[2];
      }
    }
    const Real w = waveBC.window_weight(pos);
    dvx *= w;
    dvy *= w;
    dvz *= w;
  }

  if (waveBC.turb.active && waveBC.turb.initialized) {
    waveBC.eval_turbulent_velocity_kick(pos, t, charge, dvx, dvy, dvz,
                                        waveFaceMask);
  }
}

//==========================================================
// Inner body field boundary (see #BODY / #BODYBOUNDARY).
//
// All these operators work on the nodes/cells flagged with bit::iBody_ and use
// the radial direction from the body center as the surface normal, which is
// the smooth spherical normal rather than the staircase face normal.
//==========================================================

void Pic::zero_body_E(amrex::MultiFab& mf, const int iLev) {
  // linetied: the body nodes are not part of the plasma domain and carry no
  // electric field.
  mask_body(mf, node_status(iLev));
}

//==========================================================
void Pic::zero_body_interior_E(amrex::MultiFab& mf, const int iLev) {
  // The interior of the body is a cavity: it is shielded by the surface layer
  // and carries no electric field, whatever the field boundary condition is.
  mask_body_interior(mf, node_status(iLev));
}

//==========================================================
void Pic::project_body_E(amrex::MultiFab& mf, const int iLev) {
  // conducting: E <- (E.n) n, i.e. the tangential electric field vanishes and
  // the radial component is kept. This is a surface condition: it is applied
  // to the one-cell-thick surface layer of the body, while the interior is
  // field free.
  if (mf.nComp() < 3)
    return;

  const auto plo = Geom(iLev).ProbLo();
  const auto dx = Geom(iLev).CellSize();
  const Real cx = get_body_center()[ix_];
  const Real cy = get_body_center()[iy_];
  const Real cz = get_body_center()[iz_];
  const auto& status = node_status(iLev);

  const int activeDim = get_dim();

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& box = mfi.fabbox();
    auto arr = mf[mfi].array();
    const auto statusArr = status[mfi].array();

    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (!bit::is_body(statusArr(i, j, k)) ||
          bit::is_body_interior(statusArr(i, j, k)))
        return;

      const Real x = plo[ix_] + i * dx[ix_] - cx;
      const Real y = (activeDim > 1) ? plo[iy_] + j * dx[iy_] - cy : 0.0;
      const Real z = (activeDim > 2) ? plo[iz_] + k * dx[iz_] - cz : 0.0;
      const Real r = std::sqrt(x * x + y * y + z * z);
      if (r <= 0.0)
        return;

      const Real nx = x / r;
      const Real ny = y / r;
      const Real nz = z / r;

      const Real ex = arr(i, j, k, ix_);
      const Real ey = arr(i, j, k, iy_);
      const Real ez = arr(i, j, k, iz_);
      const Real er = ex * nx + ey * ny + ez * nz;

      arr(i, j, k, ix_) = er * nx;
      arr(i, j, k, iy_) = er * ny;
      arr(i, j, k, iz_) = er * nz;
    });
  }
}

//==========================================================
void Pic::project_body_B(amrex::MultiFab& mf, const int iLev) {
  // conducting: the TOTAL field is made tangential, (B1 + B0).n = 0, so the
  // evolved field picks up B1.n = -B0.n and the surface current
  // n x (B1 + B0) closes on the intrinsic field instead of being shorted out.
  // Without B0 this is the plain B1 <- B1 - (B1.n) n. BATSRUS splits it too.
  if (mf.nComp() < 3)
    return;

  const bool isCell = mf.ixType().cellCentered();
  const auto& status = isCell ? cell_status(iLev) : node_status(iLev);
  const auto plo = Geom(iLev).ProbLo();
  const auto dx = Geom(iLev).CellSize();
  const Real cx = get_body_center()[ix_];
  const Real cy = get_body_center()[iy_];
  const Real cz = get_body_center()[iz_];
  // Cell centers sit half a cell above the lower node.
  const Real shift = isCell ? 0.5 : 0.0;
  const int activeDim = get_dim();

  const amrex::MultiFab* b0 = nullptr;
  if (use_intrinsic_B()) {
    const auto& src = isCell ? centerB0 : nodeB0;
    if (iLev < static_cast<int>(src.size()) && !src[iLev].empty()) {
      b0 = &src[iLev];
    }
  }
  const bool hasB0 = (b0 != nullptr);

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& box = mfi.fabbox();
    auto arr = mf[mfi].array();
    const auto statusArr = status[mfi].array();
    const Array4<amrex::Real const> b0Arr =
        b0 ? (*b0)[mfi].array() : Array4<amrex::Real const>();

    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (!bit::is_body(statusArr(i, j, k)) ||
          bit::is_body_interior(statusArr(i, j, k)))
        return;

      const Real x = plo[ix_] + (i + shift) * dx[ix_] - cx;
      const Real y =
          (activeDim > 1) ? plo[iy_] + (j + shift) * dx[iy_] - cy : 0.0;
      const Real z =
          (activeDim > 2) ? plo[iz_] + (k + shift) * dx[iz_] - cz : 0.0;
      const Real r = std::sqrt(x * x + y * y + z * z);
      if (r <= 0.0)
        return;

      const Real nx = x / r;
      const Real ny = y / r;
      const Real nz = z / r;

      const Real bx = arr(i, j, k, ix_);
      const Real by = arr(i, j, k, iy_);
      const Real bz = arr(i, j, k, iz_);
      // Normal component of the total field, so that B1.n = -B0.n.
      Real br = bx * nx + by * ny + bz * nz;
      if (hasB0) {
        br += b0Arr(i, j, k, ix_) * nx + b0Arr(i, j, k, iy_) * ny +
              b0Arr(i, j, k, iz_) * nz;
      }

      arr(i, j, k, ix_) = bx - br * nx;
      arr(i, j, k, iy_) = by - br * ny;
      arr(i, j, k, iz_) = bz - br * nz;
    });
  }
}

//==========================================================
void Pic::fill_body_E_insulating(amrex::MultiFab& mf, const int iLev) {
  // In hybrid PIC, no plasma exists inside the body. An uncharged dielectric
  // obstacle does not shield electric fields; E satisfies Laplace's equation
  // nabla^2(E) = 0 inside the body, matching the surrounding plasma field.
  if (mf.nComp() < 3)
    return;

  const auto& status = node_status(iLev);
  const int activeDim = get_dim();

  if (smoothScratchMF[iLev].empty() ||
      smoothScratchMF[iLev].boxArray() != mf.boxArray() ||
      smoothScratchMF[iLev].DistributionMap() != mf.DistributionMap() ||
      smoothScratchMF[iLev].nComp() != mf.nComp() ||
      smoothScratchMF[iLev].nGrow() != mf.nGrow()) {
    smoothScratchMF[iLev].define(mf.boxArray(), mf.DistributionMap(),
                                 mf.nComp(), mf.nGrow());
  }

  MultiFab& mfOld = smoothScratchMF[iLev];

  for (int iter = 0; iter < 30; ++iter) {
    MultiFab::Copy(mfOld, mf, 0, 0, mf.nComp(), mf.nGrow());
    mfOld.FillBoundary(Geom(iLev).periodicity());

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      auto arr = mf[mfi].array();
      const auto arrOld = mfOld[mfi].array();
      const auto statusArr = status[mfi].array();

      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (!bit::is_body(statusArr(i, j, k)))
          return;

        Real numNeighbors = 0.0;
        Real sumE[3] = { 0.0, 0.0, 0.0 };

        sumE[0] += arrOld(i + 1, j, k, 0) + arrOld(i - 1, j, k, 0);
        sumE[1] += arrOld(i + 1, j, k, 1) + arrOld(i - 1, j, k, 1);
        sumE[2] += arrOld(i + 1, j, k, 2) + arrOld(i - 1, j, k, 2);
        numNeighbors += 2.0;

        if (activeDim > 1) {
          sumE[0] += arrOld(i, j + 1, k, 0) + arrOld(i, j - 1, k, 0);
          sumE[1] += arrOld(i, j + 1, k, 1) + arrOld(i, j - 1, k, 1);
          sumE[2] += arrOld(i, j + 1, k, 2) + arrOld(i, j - 1, k, 2);
          numNeighbors += 2.0;
        }
        if (activeDim > 2) {
          sumE[0] += arrOld(i, j, k + 1, 0) + arrOld(i, j, k - 1, 0);
          sumE[1] += arrOld(i, j, k + 1, 1) + arrOld(i, j, k - 1, 1);
          sumE[2] += arrOld(i, j, k + 1, 2) + arrOld(i, j, k - 1, 2);
          numNeighbors += 2.0;
        }

        if (numNeighbors > 0.0) {
          arr(i, j, k, 0) = sumE[0] / numNeighbors;
          arr(i, j, k, 1) = sumE[1] / numNeighbors;
          arr(i, j, k, 2) = sumE[2] / numNeighbors;
        }
      });
    }
    mf.FillBoundary(Geom(iLev).periodicity());
  }
}

//==========================================================
void Pic::apply_body_E_bc(amrex::MultiFab& mf, const int iLev) {
  if (!use_body())
    return;

  if (bodyFieldBC == BodyFieldBC::linetied) {
    // E = 0 on every body node: the surface is line-tied and the interior is a
    // field-free cavity.
    zero_body_E(mf, iLev);
  } else if (bodyFieldBC == BodyFieldBC::conducting) {
    // E = 0 in the interior, E_t = 0 on the surface layer.
    zero_body_interior_E(mf, iLev);
    project_body_E(mf, iLev);
  } else if (bodyFieldBC == BodyFieldBC::insulating) {
    // In hybrid PIC, harmonic relaxation fills E into the body cavity.
    if (useHybridPIC) {
      fill_body_E_insulating(mf, iLev);
    }
  } else if (bodyFieldBC == BodyFieldBC::ionosphere) {
#ifdef _MI_COUPLING_
    apply_body_E_ionosphere(mf, iLev);
#else
    Abort("BodyFieldBC::ionosphere requires ./Config.pl -mi.");
#endif
  }
}

#ifdef _MI_COUPLING_
namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_curl_at_node(
    const Array4<const Real>& c, int i, int j, int k, int nDimVal,
    const GpuArray<Real, AMREX_SPACEDIM>& invDx, Real out[3]) noexcept {
  constexpr int ix = 0, iy = 1, iz = 2;
  const int km1 = (nDimVal > 2) ? (k - 1) : k;
  const Real invDxX = invDx[0];
  const Real invDxY = (AMREX_SPACEDIM > 1) ? invDx[1] : Real(0.0);
#if AMREX_SPACEDIM > 2
  const Real invDxZ = invDx[2];
#else
  const Real invDxZ = Real(0.0);
#endif

  const Real cZDY =
      Real(0.25) * invDxY *
      (c(i, j, k, iz) - c(i, j - 1, k, iz) + c(i, j, km1, iz) -
       c(i, j - 1, km1, iz) + c(i - 1, j, k, iz) - c(i - 1, j - 1, k, iz) +
       c(i - 1, j, km1, iz) - c(i - 1, j - 1, km1, iz));

  const Real cYDZ =
      (nDimVal > 2) ? Real(0.25) * invDxZ *
                          (c(i, j, k, iy) - c(i, j, km1, iy) +
                           c(i - 1, j, k, iy) - c(i - 1, j, km1, iy) +
                           c(i, j - 1, k, iy) - c(i, j - 1, km1, iy) +
                           c(i - 1, j - 1, k, iy) - c(i - 1, j - 1, km1, iy))
                    : Real(0.0);

  const Real cXDZ =
      (nDimVal > 2) ? Real(0.25) * invDxZ *
                          (c(i, j, k, ix) - c(i, j, km1, ix) +
                           c(i - 1, j, k, ix) - c(i - 1, j, km1, ix) +
                           c(i, j - 1, k, ix) - c(i, j - 1, km1, ix) +
                           c(i - 1, j - 1, k, ix) - c(i - 1, j - 1, km1, ix))
                    : Real(0.0);

  const Real cZDX =
      Real(0.25) * invDxX *
      (c(i, j, k, iz) - c(i - 1, j, k, iz) + c(i, j, km1, iz) -
       c(i - 1, j, km1, iz) + c(i, j - 1, k, iz) - c(i - 1, j - 1, k, iz) +
       c(i, j - 1, km1, iz) - c(i - 1, j - 1, km1, iz));

  const Real cYDX =
      Real(0.25) * invDxX *
      (c(i, j, k, iy) - c(i - 1, j, k, iy) + c(i, j, km1, iy) -
       c(i - 1, j, km1, iy) + c(i, j - 1, k, iy) - c(i - 1, j - 1, k, iy) +
       c(i, j - 1, km1, iy) - c(i - 1, j - 1, km1, iy));

  const Real cXDY =
      Real(0.25) * invDxY *
      (c(i, j, k, ix) - c(i, j - 1, k, ix) + c(i, j, km1, ix) -
       c(i, j - 1, km1, ix) + c(i - 1, j, k, ix) - c(i - 1, j - 1, k, ix) +
       c(i - 1, j, km1, ix) - c(i - 1, j - 1, km1, ix));

  out[0] = cZDY - cYDZ;
  out[1] = cXDZ - cZDX;
  out[2] = cYDX - cXDY;
}

} // namespace

//==========================================================
void Pic::init_ionosphere() {
  if (!is_body_ionosphere() || !intrinsicB_ || !intrinsicB_->use_dipole())
    return;

  if (!ionosphere_)
    ionosphere_ = std::make_unique<Ionosphere>(ionoParams_);
  else
    ionosphere_->set_params(ionoParams_);

  const auto& b0p = intrinsicB_->get_params();
  const Real rPlanetSI = (fi->get_rPlanet_SI() > 0.0)
                             ? fi->get_rPlanet_SI()
                             : (get_body_radius() * fi->get_No2SiL());
  const Real si2NoL = fi->get_Si2NoL();
  const Real si2NoB = fi->get_Si2NoB();
  const Real si2NoE =
      (fi->get_No2SiE() > 0.0) ? (Real(1.0) / fi->get_No2SiE()) : Real(1.0);
  const Real bEqSI = intrinsicB_->get_dipole_strength_code() * fi->get_No2SiB();

  // Default FAC extraction shell: rBody + 4*dx0 if not explicitly specified.
  if (ionosphere_->params().rFacCode <= 0.0) {
    const Real dx0 = Geom(0).CellSize(ix_);
    ionosphere_->params().rFacCode = get_body_radius() + Real(4.0) * dx0;
  }

  ionosphere_->init(b0p.moment, b0p.center, rPlanetSI, si2NoL, si2NoE, si2NoB,
                    bEqSI);
  lastIonoSolveTimeSI_ = -1.0;
  lastIonoSolveCycle_ = -1;
  facShellCache_.clear();
  bodyIonoCache_.clear();

  Print() << "  Ionosphere initialized: rIono=" << ionosphere_->frame().rIonoSI
          << " m (" << ionosphere_->frame().rIonoCode
          << " code), rFac=" << ionosphere_->get_params().rFacCode
          << " code, grid=" << ionosphere_->n_theta() << "x"
          << ionosphere_->n_phi() << "\n";
}

//==========================================================
void Pic::ensure_fac_shell_cache() {
  const int nLev = n_lev();
  const int nTh = ionosphere_->n_theta();
  const int nPh = ionosphere_->n_phi();
  const int nTot = nTh * nPh;
  const Real rFacCode = ionosphere_->get_params().rFacCode;

  bool needRebuild = !facShellCache_.initialized ||
                     facShellCache_.nTotal != nTot ||
                     facShellCache_.lastFacRadius != rFacCode ||
                     static_cast<int>(facShellCache_.levels.size()) != nLev;

  if (!needRebuild) {
    for (int iLev = 0; iLev < nLev; ++iLev) {
      const auto& lv = facShellCache_.levels[iLev];
      if (!lv.valid || lv.ba != centerB[iLev].boxArray() ||
          lv.dm != centerB[iLev].DistributionMap()) {
        needRebuild = true;
        break;
      }
    }
  }

  if (!needRebuild)
    return;

  facShellCache_.clear();
  facShellCache_.nTotal = nTot;
  facShellCache_.lastFacRadius = rFacCode;
  facShellCache_.jrFullHost.assign(static_cast<size_t>(nTot), Real(0.0));
  facShellCache_.levels.resize(nLev);

  const auto& frame = ionosphere_->frame();
  const Real latBndRad =
      ionosphere_->get_params().latBoundaryDeg * (dPI / Real(180.0));
  const auto& b0p = intrinsicB_->get_params();
  const auto& thArr = ionosphere_->theta();
  const auto& phArr = ionosphere_->phi();
  const int finest = finestLevel();

  // Pre-build per-level, per-box host vectors of static FacSample descriptors.
  std::vector<std::vector<std::vector<FacSample> > > hostSamples(nLev);
  std::vector<HostFabCache> statusCaches(nLev);
  for (int iLev = 0; iLev < nLev; ++iLev) {
    auto& lv = facShellCache_.levels[iLev];
    lv.ba = centerB[iLev].boxArray();
    lv.dm = centerB[iLev].DistributionMap();
    lv.valid = true;
    const int nLocal = centerB[iLev].local_size();
    lv.boxSamples.resize(nLocal);
    hostSamples[iLev].resize(nLocal);
  }

  int localSampleIdx = 0;
  std::vector<int> localIonoIdx;

  // Map each ionospheric (iTheta, iPhi) point once per grid/regrid and assign
  // it to its unique finest uncovered local valid box.
  for (int i = 0; i < nTh; ++i) {
    const Real th = thArr[i];
    if (std::abs(Real(0.5) * dPI - th) <= latBndRad)
      continue;
    const Real sn = std::sin(th);
    const Real cs = std::cos(th);

    for (int j = 0; j < nPh - 1; ++j) {
      const Real ph = phArr[j];
      const Real xIonoSmg[3] = {
        frame.rIonoCode * sn * std::cos(ph),
        frame.rIonoCode * sn * std::sin(ph),
        frame.rIonoCode * cs,
      };
      Real xFacSmg[3], thFac = 0.0, phFac = 0.0, jDum[2][3];
      const int iHem = map_dipole_smg(xIonoSmg, rFacCode, xFacSmg, thFac, phFac,
                                      false, jDum);
      if (iHem == 0)
        continue;

      Real xFacSim[3], xIonoSim[3];
      for (int b = 0; b < 3; ++b) {
        xFacSim[b] = frame.center[b] + frame.rot[0][b] * xFacSmg[0] +
                     frame.rot[1][b] * xFacSmg[1] +
                     frame.rot[2][b] * xFacSmg[2];
        xIonoSim[b] = frame.center[b] + frame.rot[0][b] * xIonoSmg[0] +
                      frame.rot[1][b] * xIonoSmg[1] +
                      frame.rot[2][b] * xIonoSmg[2];
      }

      Real bIono[3] = { 0.0, 0.0, 0.0 };
      eval_dipole_b(b0p, xIonoSim[0], xIonoSim[1], xIonoSim[2], bIono);
      const Real bIonoMag = std::sqrt(
          bIono[0] * bIono[0] + bIono[1] * bIono[1] + bIono[2] * bIono[2]);
      if (bIonoMag <= Real(1.0e-30))
        continue;

      const Real rHat[3] = {
        (xIonoSim[0] - frame.center[0]) / frame.rIonoCode,
        (xIonoSim[1] - frame.center[1]) / frame.rIonoCode,
        (xIonoSim[2] - frame.center[2]) / frame.rIonoCode,
      };
      const Real bIonoR =
          bIono[0] * rHat[0] + bIono[1] * rHat[1] + bIono[2] * rHat[2];

      Real b0Fac[3] = { 0.0, 0.0, 0.0 };
      eval_dipole_b(b0p, xFacSim[0], xFacSim[1], xFacSim[2], b0Fac);

      bool claimed = false;
      for (int iLev = finest; iLev >= 0 && !claimed; --iLev) {
        const auto plo = Geom(iLev).ProbLo();
        const auto invDx = Geom(iLev).InvCellSize();
        const Real gx = (xFacSim[0] - plo[0]) * invDx[0];
        const Real gy = (xFacSim[1] - plo[1]) * invDx[1];
        const Real gz = (xFacSim[2] - plo[2]) * invDx[2];
        const int ix = static_cast<int>(std::floor(gx));
        const int iy = static_cast<int>(std::floor(gy));
        const int iz = static_cast<int>(std::floor(gz));
        const IntVect iv(AMREX_D_DECL(ix, iy, iz));
        const Box curlNeed(IntVect(AMREX_D_DECL(ix - 1, iy - 1, iz - 1)),
                           IntVect(AMREX_D_DECL(ix + 1, iy + 1, iz + 1)));

        const bool checkRefined = (iLev < finest);
        const auto& statusMF = cell_status(iLev);

        for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
          if (!mfi.validbox().contains(iv) || !mfi.fabbox().contains(curlNeed))
            continue;
          if (checkRefined) {
            const auto stArr = statusCaches[iLev].array(statusMF, mfi);
            if (bit::is_refined(stArr(iv)))
              continue;
          }

          FacSample s;
          s.ix = ix;
          s.iy = iy;
          s.iz = iz;
          s.wx = gx - ix;
          s.wy = gy - iy;
          s.wz = gz - iz;
          s.ionoIdx = localSampleIdx++;
          s.b0x = b0Fac[0];
          s.b0y = b0Fac[1];
          s.b0z = b0Fac[2];
          s.bIonoR = bIonoR;

          hostSamples[iLev][mfi.LocalIndex()].push_back(s);
          localIonoIdx.push_back(i * nPh + j);
          claimed = true;
          break;
        }
      }
    }
  }

  for (int iLev = 0; iLev < nLev; ++iLev) {
    auto& lv = facShellCache_.levels[iLev];
    for (int li = 0; li < static_cast<int>(hostSamples[iLev].size()); ++li) {
      const auto& src = hostSamples[iLev][li];
      if (src.empty())
        continue;
      lv.boxSamples[li].resize(src.size());
      Gpu::copyAsync(Gpu::hostToDevice, src.begin(), src.end(),
                     lv.boxSamples[li].begin());
    }
  }
  Gpu::streamSynchronize();

  facShellCache_.nLocalActive = localSampleIdx;
  facShellCache_.localIonoIdx = std::move(localIonoIdx);
  facShellCache_.jrDevBuf.resize(static_cast<size_t>(localSampleIdx));
  facShellCache_.jrCompactHost.resize(static_cast<size_t>(localSampleIdx));
  facShellCache_.initialized = true;
}

//==========================================================
void Pic::extract_ionosphere_fac(std::vector<Real>& jrOutwardSI) {
  ensure_fac_shell_cache();

  const int nTh = ionosphere_->n_theta();
  const int nPh = ionosphere_->n_phi();
  const int nTot = facShellCache_.nTotal;
  const int nLocal = facShellCache_.nLocalActive;

  constexpr Real mu0 = Real(4.0) * dPI * Real(1.0e-7);
  const Real curlBToJSI = fi->get_No2SiB() / (mu0 * fi->get_No2SiL());
  const int nDimVal = get_dim();

  if (nLocal > 0) {
    Real* jrDevPtr = facShellCache_.jrDevBuf.data();
    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      const auto& lv = facShellCache_.levels[iLev];
      const auto invDx = Geom(iLev).InvCellSizeArray();

      for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
        const auto& devSamples = lv.boxSamples[mfi.LocalIndex()];
        const int nS = static_cast<int>(devSamples.size());
        if (nS == 0)
          continue;

        const FacSample* sPtr = devSamples.data();
        const auto b1Arr = centerB[iLev][mfi].const_array();

        ParallelFor(nS, [=] AMREX_GPU_DEVICE(int idx) noexcept {
          const FacSample s = sPtr[idx];

          // Inline curl(B1) on the 8 corners of the containing cell.
          Real c000[3], c100[3], c010[3], c110[3];
          Real c001[3], c101[3], c011[3], c111[3];
          eval_curl_at_node(b1Arr, s.ix, s.iy, s.iz, nDimVal, invDx, c000);
          eval_curl_at_node(b1Arr, s.ix + 1, s.iy, s.iz, nDimVal, invDx, c100);
          eval_curl_at_node(b1Arr, s.ix, s.iy + 1, s.iz, nDimVal, invDx, c010);
          eval_curl_at_node(b1Arr, s.ix + 1, s.iy + 1, s.iz, nDimVal, invDx,
                            c110);
          eval_curl_at_node(b1Arr, s.ix, s.iy, s.iz + 1, nDimVal, invDx, c001);
          eval_curl_at_node(b1Arr, s.ix + 1, s.iy, s.iz + 1, nDimVal, invDx,
                            c101);
          eval_curl_at_node(b1Arr, s.ix, s.iy + 1, s.iz + 1, nDimVal, invDx,
                            c011);
          eval_curl_at_node(b1Arr, s.ix + 1, s.iy + 1, s.iz + 1, nDimVal, invDx,
                            c111);

          Real jSI[3];
          const Real omx = Real(1.0) - s.wx;
          const Real omy = Real(1.0) - s.wy;
          const Real omz = Real(1.0) - s.wz;
          for (int c = 0; c < 3; ++c) {
            const Real v00 = omx * c000[c] + s.wx * c100[c];
            const Real v10 = omx * c010[c] + s.wx * c110[c];
            const Real v01 = omx * c001[c] + s.wx * c101[c];
            const Real v11 = omx * c011[c] + s.wx * c111[c];
            const Real v0 = omy * v00 + s.wy * v10;
            const Real v1 = omy * v01 + s.wy * v11;
            jSI[c] = (omz * v0 + s.wz * v1) * curlBToJSI;
          }

          const Real bFac0 = b1Arr(s.ix, s.iy, s.iz, 0) + s.b0x;
          const Real bFac1 = b1Arr(s.ix, s.iy, s.iz, 1) + s.b0y;
          const Real bFac2 = b1Arr(s.ix, s.iy, s.iz, 2) + s.b0z;
          const Real bFacMag2 = bFac0 * bFac0 + bFac1 * bFac1 + bFac2 * bFac2;

          Real jrVal = Real(0.0);
          if (bFacMag2 > Real(1.0e-60)) {
            const Real jDotB = jSI[0] * bFac0 + jSI[1] * bFac1 + jSI[2] * bFac2;
            jrVal = (jDotB / bFacMag2) * s.bIonoR;
          }
          jrDevPtr[s.ionoIdx] = jrVal;
        });
      }
    }

    Gpu::copy(Gpu::deviceToHost, facShellCache_.jrDevBuf.begin(),
              facShellCache_.jrDevBuf.end(),
              facShellCache_.jrCompactHost.begin());
  }

  jrOutwardSI.assign(static_cast<size_t>(nTot), Real(0.0));
  for (int idx = 0; idx < nLocal; ++idx) {
    jrOutwardSI[static_cast<size_t>(facShellCache_.localIonoIdx[idx])] =
        facShellCache_.jrCompactHost[idx];
  }

  ParallelDescriptor::ReduceRealSum(jrOutwardSI.data(), nTot,
                                    ParallelDescriptor::IOProcessorNumber());

  if (ParallelDescriptor::IOProcessor()) {
    for (int i = 0; i < nTh; ++i) {
      jrOutwardSI[static_cast<size_t>(i * nPh + (nPh - 1))] =
          jrOutwardSI[static_cast<size_t>(i * nPh)];
    }
  }
}

//==========================================================
void Pic::update_ionosphere(bool force) {
  if (!is_body_ionosphere() || !ionosphere_ || !ionosphere_->is_initialized())
    return;

  const int cycle = tc ? tc->get_cycle() : 0;
  const Real timeSI = tc ? (tc->get_time() * fi->get_No2SiT()) : Real(0.0);
  const auto& p = ionosphere_->get_params();

  if (!force) {
    if (p.dtCoupleSI > Real(0.0)) {
      if (lastIonoSolveTimeSI_ >= Real(0.0) &&
          (timeSI - lastIonoSolveTimeSI_) < p.dtCoupleSI)
        return;
    } else if (p.dnCouple > 1) {
      if (lastIonoSolveCycle_ >= 0 &&
          (cycle - lastIonoSolveCycle_) < p.dnCouple)
        return;
    }
  }

  std::string nameFunc = "Pic::update_ionosphere";
  timing_func(nameFunc);

  extract_ionosphere_fac(facShellCache_.jrFullHost);
  const bool ok = ionosphere_->solve(facShellCache_.jrFullHost, true);

  lastIonoSolveTimeSI_ = timeSI;
  lastIonoSolveCycle_ = cycle;

  if (!ok) {
    Print() << "  [MI WARNING] Ionosphere solver did not reach tolerance "
            << "(iterN=" << ionosphere_->last_iter_north()
            << ", residN=" << ionosphere_->last_resid_north()
            << ", iterS=" << ionosphere_->last_iter_south()
            << ", residS=" << ionosphere_->last_resid_south() << ")\n";
  }

  if (p.dnOutput > 0 && (force || (cycle % p.dnOutput == 0))) {
    char buf[256];
    std::snprintf(buf, sizeof(buf), "%s/plots/ionosphere_n%07d.dat",
                  component.c_str(), cycle);
    ionosphere_->write_diagnostics(buf, timeSI, cycle);
  }
}

//==========================================================
void Pic::ensure_body_iono_cache(const int iLev) {
  const int nLev = n_lev();
  if (static_cast<int>(bodyIonoCache_.size()) != nLev)
    bodyIonoCache_.resize(nLev);

  const auto& status = node_status(iLev);
  auto& cache = bodyIonoCache_[iLev];
  if (cache.valid && cache.ba == status.boxArray() &&
      cache.dm == status.DistributionMap())
    return;

  cache.clear();
  cache.ba = status.boxArray();
  cache.dm = status.DistributionMap();
  cache.localToBodyIdx.assign(status.local_size(), -1);

  const auto plo = Geom(iLev).ProbLo();
  const auto dx = Geom(iLev).CellSize();
  const int activeDim = get_dim();
  const IonosphereDeviceView view = ionosphere_->device_view();

  HostFabCache statusCache;
  for (MFIter mfi(status); mfi.isValid(); ++mfi) {
    const Box& fbx = mfi.fabbox();
    const auto stArr = statusCache.array(status, mfi);

    IntVect intLo(AMREX_D_DECL(0, 0, 0)), intHi(AMREX_D_DECL(-1, -1, -1));
    bool hasInt = false;
    std::vector<BodySurfaceNode> surfNodes;

    const auto lo = lbound(fbx);
    const auto hi = ubound(fbx);
    for (int k = lo.z; k <= hi.z; ++k) {
      const Real z = (activeDim > 2) ? plo[iz_] + k * dx[iz_] : Real(0.0);
      for (int j = lo.y; j <= hi.y; ++j) {
        const Real y = (activeDim > 1) ? plo[iy_] + j * dx[iy_] : Real(0.0);
        for (int i = lo.x; i <= hi.x; ++i) {
          const int st = stArr(i, j, k);
          if (bit::is_body_interior(st)) {
            const IntVect iv(AMREX_D_DECL(i, j, k));
            if (!hasInt) {
              intLo = iv;
              intHi = iv;
              hasInt = true;
            } else {
              intLo.min(iv);
              intHi.max(iv);
            }
            continue;
          }
          if (!bit::is_body(st))
            continue;

          BodySurfaceNode node{};
          node.i = i;
          node.j = j;
          node.k = k;
          node.idx00 = -1;

          const Real xSim = plo[ix_] + i * dx[ix_];
          const Real dVec[3] = {
            xSim - view.frame.center[0],
            y - view.frame.center[1],
            z - view.frame.center[2],
          };
          Real xSmg[3] = { 0.0, 0.0, 0.0 };
          for (int a = 0; a < 3; ++a) {
            xSmg[a] = view.frame.rot[a][0] * dVec[0] +
                      view.frame.rot[a][1] * dVec[1] +
                      view.frame.rot[a][2] * dVec[2];
          }

          Real xMap[3] = { 0.0, 0.0, 0.0 };
          Real theta = 0.0, phi = 0.0;
          Real jacSmg[2][3];
          const int iHem = map_dipole_smg(xSmg, view.frame.rIonoCode, xMap,
                                          theta, phi, true, jacSmg);
          const Real magLat = std::abs(Real(0.5) * dPI - theta);
          if (iHem != 0 && magLat > view.latBoundaryRad && view.nTheta >= 3 &&
              view.nPhi >= 3) {
            const Real tNorm = theta * view.invDTheta;
            const Real pNorm = phi * view.invDPhi;
            int i0 = static_cast<int>(std::floor(tNorm));
            int j0 = static_cast<int>(std::floor(pNorm));
            if (i0 < 0)
              i0 = 0;
            if (i0 >= view.nTheta - 1)
              i0 = view.nTheta - 2;
            if (j0 < 0)
              j0 = 0;
            if (j0 >= view.nPhi - 1)
              j0 = view.nPhi - 2;

            const Real wt = tNorm - i0;
            const Real wp = pNorm - j0;

            node.idx00 = i0 * view.nPhi + j0;
            node.idx01 = i0 * view.nPhi + (j0 + 1);
            node.idx10 = (i0 + 1) * view.nPhi + j0;
            node.idx11 = (i0 + 1) * view.nPhi + (j0 + 1);

            node.w00 = (Real(1.0) - wt) * (Real(1.0) - wp);
            node.w01 = (Real(1.0) - wt) * wp;
            node.w10 = wt * (Real(1.0) - wp);
            node.w11 = wt * wp;

            // Pre-rotate and scale SMG dipole Jacobian to simulation E units:
            // jTh = eScale * R^T * jacSmg[0], jPh = eScale * R^T * jacSmg[1].
            node.jTh0 = view.eScale * (view.frame.rot[0][0] * jacSmg[0][0] +
                                       view.frame.rot[1][0] * jacSmg[0][1] +
                                       view.frame.rot[2][0] * jacSmg[0][2]);
            node.jTh1 = view.eScale * (view.frame.rot[0][1] * jacSmg[0][0] +
                                       view.frame.rot[1][1] * jacSmg[0][1] +
                                       view.frame.rot[2][1] * jacSmg[0][2]);
            node.jTh2 = view.eScale * (view.frame.rot[0][2] * jacSmg[0][0] +
                                       view.frame.rot[1][2] * jacSmg[0][1] +
                                       view.frame.rot[2][2] * jacSmg[0][2]);

            node.jPh0 = view.eScale * (view.frame.rot[0][0] * jacSmg[1][0] +
                                       view.frame.rot[1][0] * jacSmg[1][1] +
                                       view.frame.rot[2][0] * jacSmg[1][2]);
            node.jPh1 = view.eScale * (view.frame.rot[0][1] * jacSmg[1][0] +
                                       view.frame.rot[1][1] * jacSmg[1][1] +
                                       view.frame.rot[2][1] * jacSmg[1][2]);
            node.jPh2 = view.eScale * (view.frame.rot[0][2] * jacSmg[1][0] +
                                       view.frame.rot[1][2] * jacSmg[1][1] +
                                       view.frame.rot[2][2] * jacSmg[1][2]);
          }

          surfNodes.push_back(node);
        }
      }
    }

    const Box intSubBox = hasInt ? Box(intLo, intHi, fbx.ixType()) : Box();
    if (!intSubBox.isEmpty() || !surfNodes.empty()) {
      const int bi = static_cast<int>(cache.bodyBoxes.size());
      cache.localToBodyIdx[mfi.LocalIndex()] = bi;
      BodyBoxCache bCache;
      bCache.localIndex = mfi.LocalIndex();
      bCache.interiorSubBox = intSubBox;
      if (!surfNodes.empty()) {
        bCache.surfaceNodes.resize(surfNodes.size());
        Gpu::copyAsync(Gpu::hostToDevice, surfNodes.begin(), surfNodes.end(),
                       bCache.surfaceNodes.begin());
      }
      cache.bodyBoxes.push_back(std::move(bCache));
    }
  }

  Gpu::streamSynchronize();
  cache.cachedGeneration = 0;
  cache.valid = true;
}

//==========================================================
void Pic::refresh_body_iono_E(const int iLev) {
  auto& cache = bodyIonoCache_[iLev];
  const IonosphereDeviceView view = ionosphere_->device_view();
  if (!view.dPhiDTheta || !view.dPhiDPhi)
    return;

  for (auto& bCache : cache.bodyBoxes) {
    const int nSurf = static_cast<int>(bCache.surfaceNodes.size());
    if (nSurf == 0)
      continue;

    BodySurfaceNode* nodes = bCache.surfaceNodes.data();
    ParallelFor(nSurf, [=] AMREX_GPU_DEVICE(int idx) noexcept {
      BodySurfaceNode n = nodes[idx];
      if (n.idx00 < 0)
        return;

      const Real dPdT =
          n.w00 * view.dPhiDTheta[n.idx00] + n.w01 * view.dPhiDTheta[n.idx01] +
          n.w10 * view.dPhiDTheta[n.idx10] + n.w11 * view.dPhiDTheta[n.idx11];
      const Real dPdP =
          n.w00 * view.dPhiDPhi[n.idx00] + n.w01 * view.dPhiDPhi[n.idx01] +
          n.w10 * view.dPhiDPhi[n.idx10] + n.w11 * view.dPhiDPhi[n.idx11];

      n.ex = dPdT * n.jTh0 + dPdP * n.jPh0;
      n.ey = dPdT * n.jTh1 + dPdP * n.jPh1;
      n.ez = dPdT * n.jTh2 + dPdP * n.jPh2;
      nodes[idx] = n;
    });
  }
  Gpu::streamSynchronize();
  cache.cachedGeneration = ionosphere_->generation();
}

//==========================================================
void Pic::apply_body_E_ionosphere(amrex::MultiFab& mf, const int iLev) {
  if (!ionosphere_ || !ionosphere_->is_initialized() || mf.nComp() < 3)
    return;

  ensure_body_iono_cache(iLev);
  auto& cache = bodyIonoCache_[iLev];
  if (cache.bodyBoxes.empty())
    return;

  if (cache.cachedGeneration != ionosphere_->generation())
    refresh_body_iono_E(iLev);

  const auto& status = node_status(iLev);
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const int bi = cache.localToBodyIdx[mfi.LocalIndex()];
    if (bi < 0)
      continue;

    const auto& bCache = cache.bodyBoxes[bi];
    const Box& fbx = mfi.fabbox();
    auto arr = mf[mfi].array();

    if (!bCache.interiorSubBox.isEmpty()) {
      const Box zbox = fbx & bCache.interiorSubBox;
      if (!zbox.isEmpty()) {
        const auto stArr = status[mfi].const_array();
        parallel_for_omp(zbox,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                           if (bit::is_body_interior(stArr(i, j, k))) {
                             arr(i, j, k, ix_) = Real(0.0);
                             arr(i, j, k, iy_) = Real(0.0);
                             arr(i, j, k, iz_) = Real(0.0);
                           }
                         });
      }
    }

    const int nSurf = static_cast<int>(bCache.surfaceNodes.size());
    if (nSurf > 0) {
      const BodySurfaceNode* nodes = bCache.surfaceNodes.data();
      ParallelFor(nSurf, [=] AMREX_GPU_DEVICE(int idx) noexcept {
        const BodySurfaceNode n = nodes[idx];
        if (!fbx.contains(IntVect(AMREX_D_DECL(n.i, n.j, n.k))))
          return;
        arr(n.i, n.j, n.k, ix_) = n.ex;
        arr(n.i, n.j, n.k, iy_) = n.ey;
        arr(n.i, n.j, n.k, iz_) = n.ez;
      });
    }
  }
}
#endif

//==========================================================
void Pic::warn_absorb_cfl() {
  const Real cs = (absorbCharSpeed > 0.0) ? absorbCharSpeed : 1.0;
  const Real* dx = Geom(0).CellSize();
  Real cfl = 0.0;
  for (int d = 0; d < nDim; ++d) {
    if (bcField.face(d, 0) == FieldBC::absorb ||
        bcField.face(d, 1) == FieldBC::absorb)
      cfl = std::max(cfl, cs * tc->get_dt() / dx[d]);
  }
  if (cfl <= 1.0)
    return;

  isAbsorbCflWarned_ = true;
  amrex::Print() << printPrefix
                 << " [CFL warning] absorbing boundary: charSpeed*dt/dx = "
                 << cfl
                 << " (> 1). The absorbing boundary is unstable above 1; "
                    "use a smaller time step or 'outflow' faces.\n";
}

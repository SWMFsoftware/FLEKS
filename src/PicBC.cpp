#include <algorithm>
#include <cmath>
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

} // namespace

//==========================================================
void Pic::apply_field_bc(const iMultiFab& status, MultiFab& mf,
                         const int iStart, const int nComp, GETVALUE func,
                         const int iLev, const bool isB) {
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

  // Fill external Dirichlet (coupled / fixed) faces from func if provided:
  if (func != nullptr) {
    fill_ext_dir(status, mf, iStart, nComp, func, iLev, &bcField);
  }

  // Node-centered conducting wall requires zeroing boundary wall nodes
  // and mirroring nodal ghost nodes.
  if (hasConductingBC_ && mf.boxArray().ixType().nodeCentered()) {
    apply_conducting_wall(status, mf, iStart, nComp, iLev, bcField, isB);
  }

  if (hasAbsorbBC_)
    apply_absorbing_wall(status, mf, iStart, nComp, iLev, bcField, isB);

  // Wave boundary condition overwrites faces where active.
  if (waveBC.active) {
    const Real t = tc ? tc->get_time() : 0.0;
    apply_wave_field(status, mf, iStart, nComp, iLev, bcField, isB ? 0 : 1, t,
                     func);
  }
}

//==========================================================
void Pic::fill_ext_dir(const iMultiFab& status, MultiFab& mf, const int iStart,
                       const int nComp, GETVALUE func, const int iLev,
                       const BoxBC<FieldBC::Type>* bc) {
  if (!func)
    return;

  const BoxArray ba =
      get_boundary_active_ba(activeRegion, mf, Geom(iLev), nDim, iz_);
  const BoundaryBounds bnd(Geom(iLev), mf.boxArray().ixType(), bc);

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    if (ba.contains(bxFab))
      continue;

    Array4<Real> const arr = mf.array(mfi);
    Array4<const int> const statusArr = status.array(mfi);

    // Note: [&mfi] is retained because host member-function pointer func
    // requires mfi. All other data are captured by value for GPU readiness.
    ParallelFor(bxFab, [=, &mfi](int i, int j, int k) {
      if (!bit::is_lev_boundary(statusArr(i, j, k, 0)))
        return;

      // If at an outer physical boundary with non-Dirichlet condition,
      // PhysBCFunct and dedicated wall operators handle it.
      if (bc != nullptr) {
        const int ijk[3] = { i, j, k };
        bool skipForPhysWall = false;
        for (int d = 0; d < nDim; ++d) {
          if ((ijk[d] < bnd.loBnd[d] && (bnd.bcLo[d] == FieldBC::outflow ||
                                         bnd.bcLo[d] == FieldBC::inflow ||
                                         bnd.bcLo[d] == FieldBC::conducting)) ||
              (ijk[d] > bnd.hiBnd[d] && (bnd.bcHi[d] == FieldBC::outflow ||
                                         bnd.bcHi[d] == FieldBC::inflow ||
                                         bnd.bcHi[d] == FieldBC::conducting))) {
            skipForPhysWall = true;
            break;
          }
        }
        if (skipForPhysWall)
          return;
      }

      for (int iVar = 0; iVar < nComp; ++iVar) {
        arr(i, j, k, iStart + iVar) =
            (this->*func)(mfi, IntVect{ AMREX_D_DECL(i, j, k) }, iVar, iLev);
      }
    });
  }
}

//==========================================================
void Pic::apply_BC(const iMultiFab& status, MultiFab& mf, const int iStart,
                   const int nComp, GETVALUE func, const int iLev,
                   const BoxBC<FieldBC::Type>* bc) {
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
      fill_ext_dir(status, mf, iStart, nComp, func, iLev, bc);
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
    // boundaries:
    const BoxArray ba =
        get_boundary_active_ba(activeRegion, mf, Geom(iLev), nDim, iz_);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bxFab = mfi.fabbox();
      const Box& bxValid = mfi.validbox();
      if (!ba.contains(bxFab)) {
        Array4<Real> const& arr = mf[mfi].array();
        const Array4<const int>& statusArr = status[mfi].array();
        Box box = bxValid;
        box.grow(1);
        ParallelFor(box, [&](int i, int j, int k) {
          if (bit::is_lev_boundary(statusArr(i, j, k, 0))) {
            bool isNeiFound = false;
            const int kmin = (nDim > 2) ? -1 : 0;
            const int kmax = (nDim > 2) ? 1 : 0;
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
    }
  } else {
    fill_ext_dir(status, mf, iStart, nComp, func, iLev, nullptr);
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
      get_boundary_active_ba(activeRegion, mf, Geom(iLev), nDim, iz_);
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

    ParallelFor(bxFab, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const int ijk[3] = { i, j, k };
      const int vLoArr[3] = { vLo.x, vLo.y, vLo.z };
      const int vHiArr[3] = { vHi.x, vHi.y, vHi.z };

      // 1. Boundary nodes on physical wall (node-centred only).
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
                               bool isB) {
  (void)isB;
  const Real dt = tc ? tc->get_dt() : 0.0;
  if (dt <= 0.0)
    return;

  std::string nameFunc = "Pic::apply_absorbing_wall";
  timing_func(nameFunc);

  // Characteristic speed; default c=1, override via #ABSORB.
  const Real cs = (absorbCharSpeed > 0.0) ? absorbCharSpeed : 1.0;

  const BoxArray ba =
      get_boundary_active_ba(activeRegion, mf, Geom(iLev), nDim, iz_);
  const BoundaryBounds bnd(Geom(iLev), mf.boxArray().ixType(), &bc);
  const Real* dx = Geom(iLev).CellSize();

  AbsorbWeights weights;
  for (int d = 0; d < nDim; ++d) {
    const Real drive0 = cs * dt / dx[d];
    weights.decay[d] = (1.0 - drive0) / (1.0 + drive0);
    weights.drive[d] = 2.0 * drive0 / (1.0 + drive0);
  }

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    if (ba.contains(bxFab))
      continue;

    Array4<Real> const& arr = mf[mfi].array();
    const Array4<const int>& statusArr = status[mfi].array();

    ParallelFor(bxFab, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      if (!bit::is_lev_boundary(statusArr(i, j, k, 0)))
        return;

      const int ijk[3] = { i, j, k };
      int m[3] = { i, j, k };
      Real cellDecay = 1.0;
      Real cellDrive = 0.0;
      int nAbsorb = 0;

      for (int d = 0; d < nDim; ++d) {
        bool isLow =
            (bnd.bcLo[d] == FieldBC::absorb) && (ijk[d] < bnd.loBnd[d]);
        bool isHigh =
            (bnd.bcHi[d] == FieldBC::absorb) && (ijk[d] > bnd.hiBnd[d]);
        if (isLow) {
          m[d] = bnd.isNode[d] ? (2 * bnd.loBnd[d] - ijk[d])
                               : (2 * bnd.loBnd[d] - 1 - ijk[d]);
          cellDecay *= weights.decay[d];
          cellDrive += weights.drive[d];
          nAbsorb++;
        } else if (isHigh) {
          m[d] = bnd.isNode[d] ? (2 * bnd.hiBnd[d] - ijk[d])
                               : (2 * bnd.hiBnd[d] + 1 - ijk[d]);
          cellDecay *= weights.decay[d];
          cellDrive += weights.drive[d];
          nAbsorb++;
        }
      }

      if (nAbsorb == 0)
        return;

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
                            const BoxBC<FieldBC::Type>& bc, bool isB) {
  (void)isB; // zero-gradient copy is component-agnostic

  std::string nameFunc = "Pic::apply_inflow_wall";
  timing_func(nameFunc);

  const BoxArray ba =
      get_boundary_active_ba(activeRegion, mf, Geom(iLev), nDim, iz_);
  const BoundaryBounds bnd(Geom(iLev), mf.boxArray().ixType(), &bc);

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    if (ba.contains(bxFab))
      continue;

    Array4<Real> const& arr = mf[mfi].array();
    const Array4<const int>& statusArr = status[mfi].array();

    ParallelFor(bxFab, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      if (!bit::is_lev_boundary(statusArr(i, j, k, 0)))
        return;

      const int ijk[3] = { i, j, k };
      int m[3] = { i, j, k };
      bool touched = false;

      for (int d = 0; d < nDim; ++d) {
        if (bnd.bcLo[d] == FieldBC::inflow && ijk[d] < bnd.loBnd[d]) {
          m[d] = bnd.loBnd[d];
          touched = true;
        } else if (bnd.bcHi[d] == FieldBC::inflow && ijk[d] > bnd.hiBnd[d]) {
          m[d] = bnd.hiBnd[d];
          touched = true;
        }
      }

      if (!touched)
        return;

      for (int iVar = 0; iVar < nComp; ++iVar) {
        arr(i, j, k, iStart + iVar) = arr(m[0], m[1], m[2], iStart + iVar);
      }
    });
  }
}

//==========================================================
void Pic::apply_wave_field(const iMultiFab& status, MultiFab& mf,
                           const int iStart, const int nComp, const int iLev,
                           const BoxBC<FieldBC::Type>& bc, int iField, Real t,
                           GETVALUE func) {
  // Pre-filter components active at time t for iField.
  // Hoisting envelope and time-phase evaluation out of the cell loops.
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

  std::vector<ActiveComp> activeComps;
  for (const auto& f : waveBC.faces) {
    if (f.direction >= nDim)
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
  if (activeComps.empty())
    return;

  std::string nameFunc = "Pic::apply_wave_field";
  timing_func(nameFunc);

  const BoxArray ba =
      get_boundary_active_ba(activeRegion, mf, Geom(iLev), nDim, iz_);
  const BoundaryBounds bnd(Geom(iLev), mf.boxArray().ixType(), &bc);

  const Real* plo = Geom(iLev).ProbLo();
  const Real* dx = Geom(iLev).CellSize();
  const Real offset[3] = { bnd.isNode[0] ? 0.0 : 0.5, bnd.isNode[1] ? 0.0 : 0.5,
                           bnd.isNode[2] ? 0.0 : 0.5 };

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    const Box& bxValid = mfi.validbox();
    if (ba.contains(bxFab))
      continue;

    // Box culling: check if this FAB touches any face active for wave
    // injection.
    bool boxTouchesWaveFace = false;
    for (const auto& ac : activeComps) {
      const int d = ac.dir;
      const int side = ac.side;
      if (side == 0) {
        if (bxFab.smallEnd(d) < bxValid.smallEnd(d) ||
            (bnd.isNode[d] && bxFab.smallEnd(d) <= bnd.loBnd[d])) {
          boxTouchesWaveFace = true;
          break;
        }
      } else {
        if (bxFab.bigEnd(d) > bxValid.bigEnd(d) ||
            (bnd.isNode[d] && bxFab.bigEnd(d) >= bnd.hiBnd[d])) {
          boxTouchesWaveFace = true;
          break;
        }
      }
    }
    if (!boxTouchesWaveFace)
      continue;

    Array4<Real> const& arr = mf[mfi].array();
    const Array4<const int>& statusArr = status[mfi].array();
    const Dim3 vLo = bxValid.smallEnd().dim3();
    const Dim3 vHi = bxValid.bigEnd().dim3();
    const int vLoArr[3] = { vLo.x, vLo.y, vLo.z };
    const int vHiArr[3] = { vHi.x, vHi.y, vHi.z };

    ParallelFor(bxFab, [&](int i, int j, int k) {
      const int ijk[3] = { i, j, k };
      const bool isGhost = bit::is_lev_boundary(statusArr(i, j, k, 0));

      Real pos[3] = { plo[0] + dx[0] * (i + offset[0]),
                      plo[1] + dx[1] * (j + offset[1]),
                      (nDim > 2) ? plo[2] + dx[2] * (k + offset[2]) : 0.0 };

      Real waveVal[3] = { 0.0, 0.0, 0.0 };
      bool hasWave = false;
      bool isBndNode = false;

      for (const auto& ac : activeComps) {
        const int d = ac.dir;
        const int side = ac.side;
        const int idx = ijk[d];
        const bool onGhostFace =
            isGhost && ((side == 0 && idx < bxValid.smallEnd(d)) ||
                        (side == 1 && idx > bxValid.bigEnd(d)));

        // Boundary node on physical wall (node-centred only).
        bool onNodeWall = false;
        if (bnd.isNode[d]) {
          const bool onLoNode = (side == 0 && idx == bnd.loBnd[d]);
          const bool onHiNode = (side == 1 && idx == bnd.hiBnd[d]);
          if (onLoNode || onHiNode) {
            bool inTangential = true;
            for (int od = 0; od < nDim; ++od) {
              if (od != d && (ijk[od] < vLoArr[od] || ijk[od] > vHiArr[od])) {
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
        if (isBndNode && func) {
          // Reset nodal physical boundary from base state before superimposing
          // wave.
          for (int iVar = 0; iVar < nComp; ++iVar) {
            const Real baseVal = (this->*func)(
                mfi, IntVect{ AMREX_D_DECL(i, j, k) }, iVar, iLev);
            arr(i, j, k, iStart + iVar) = baseVal + waveVal[iVar % 3];
          }
        } else {
          // Ghost cells already contain base state from apply_BC(); add wave
          // perturbation.
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
                             Real& dvz) {
  dvx = dvy = dvz = 0.0;
  if (!waveBC.active)
    return;
  for (const auto& f : waveBC.faces) {
    for (const auto& c : f.comps) {
      if (c.iField != 2 || c.amplitude == 0.0) // velocity kick
        continue;
      const Real val = waveBC.value(c, t, pos);
      if (val == 0.0)
        continue;
      dvx += val * c.pol[0];
      dvy += val * c.pol[1];
      dvz += val * c.pol[2];
    }
  }
}

//==========================================================

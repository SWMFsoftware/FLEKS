#include <algorithm>
#include <cmath>
#include <vector>

#include <AMReX_Loop.H>
#include <AMReX_MultiFabUtil.H>

#include "GridUtility.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

namespace {

struct BoundaryBounds {
  Dim3 domLo;
  Dim3 domHi;
  bool isNode[3] = { false, false, false };
  int loBnd[3] = { 0, 0, 0 };
  int hiBnd[3] = { 0, 0, 0 };
  int bcLo[3] = { 0, 0, 0 };
  int bcHi[3] = { 0, 0, 0 };

  BoundaryBounds() = default;
  BoundaryBounds(const Geometry& geom, IndexType ixType,
                 const BoxBC<FieldBC::Type>* bc = nullptr) {
    domLo = geom.Domain().smallEnd().dim3();
    domHi = geom.Domain().bigEnd().dim3();
    const int* dLo = geom.Domain().smallEnd().getVect();
    const int* dHi = geom.Domain().bigEnd().getVect();
    for (int d = 0; d < 3; ++d) {
      if (d < nDim) {
        isNode[d] = (ixType[d] == IndexType::NODE);
        loBnd[d] = dLo[d];
        hiBnd[d] = isNode[d] ? (dHi[d] + 1) : dHi[d];
        if (bc) {
          bcLo[d] = bc->face(d, 0);
          bcHi[d] = bc->face(d, 1);
        }
      }
    }
  }
};

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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool eval_use_float(const int i, const int j, const int k,
                    int &ip, int &jp, int &kp,
                    const BoxBC<FieldBC::Type> &bc,
                    const Box &bxValid, const int nDimLocal) {
  bool useFloat = false;
  ip = i;
  jp = j;
  kp = k;

  if (i < bxValid.smallEnd(0) &&
      (bc.lo[0] == FieldBC::outflow || bc.lo[0] == FieldBC::inflow)) {
    useFloat = true;
    ip = bxValid.smallEnd(0);
  }
  if (i > bxValid.bigEnd(0) &&
      (bc.hi[0] == FieldBC::outflow || bc.hi[0] == FieldBC::inflow)) {
    useFloat = true;
    ip = bxValid.bigEnd(0);
  }

  if (j < bxValid.smallEnd(1) &&
      (bc.lo[1] == FieldBC::outflow || bc.lo[1] == FieldBC::inflow)) {
    useFloat = true;
    jp = bxValid.smallEnd(1);
  }
  if (j > bxValid.bigEnd(1) &&
      (bc.hi[1] == FieldBC::outflow || bc.hi[1] == FieldBC::inflow)) {
    useFloat = true;
    jp = bxValid.bigEnd(1);
  }

  if (nDimLocal > 2) {
    if (k < bxValid.smallEnd(2) &&
        (bc.lo[2] == FieldBC::outflow || bc.lo[2] == FieldBC::inflow)) {
      useFloat = true;
      kp = bxValid.smallEnd(2);
    }
    if (k > bxValid.bigEnd(2) &&
        (bc.hi[2] == FieldBC::outflow || bc.hi[2] == FieldBC::inflow)) {
      useFloat = true;
      kp = bxValid.bigEnd(2);
    }
  }
  return useFloat;
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

  // Base fill: float on open faces, or evaluate state from func elsewhere.
  apply_BC(status, mf, iStart, nComp, func, iLev, &bcField);

  // Dedicated wall operators applied per configured face type.
  if (hasConductingBC_)
    apply_conducting_wall(status, mf, iStart, nComp, iLev, bcField, isB);

  if (hasAbsorbBC_)
    apply_absorbing_wall(status, mf, iStart, nComp, iLev, bcField, isB);

  if (hasInflowBC_ && fi->get_inflow_defined())
    apply_inflow_wall(status, mf, iStart, nComp, iLev, bcField, isB);

  // Wave boundary condition overwrites faces where active.
  if (waveBC.active) {
    const Real t = tc ? tc->get_time() : 0.0;
    apply_wave_field(status, mf, iStart, nComp, iLev, bcField, isB ? 0 : 1, t);
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

  const BoxArray ba =
      get_boundary_active_ba(activeRegion, mf, Geom(iLev), nDim, iz_);
  const int nDimLocal = nDim;

  constexpr int fModeNone = 0;
  constexpr int fModeZero = 1;
  constexpr int fModeCenterB = 2;
  constexpr int fModeNodeB = 3;
  constexpr int fModeFluidE = 4;

  int fMode = fModeNone;
  GETVALUE fnZero = &Pic::get_zero;
  GETVALUE fnCenterB = static_cast<amrex::Real (Pic::*)(amrex::MFIter&, amrex::IntVect, int, const int)>(&Pic::get_center_B);
  GETVALUE fnNodeB = static_cast<amrex::Real (Pic::*)(amrex::MFIter&, amrex::IntVect, int, const int)>(&Pic::get_node_B);
  GETVALUE fnNodeE = static_cast<amrex::Real (Pic::*)(amrex::MFIter&, amrex::IntVect, int, const int)>(&Pic::get_node_E);
  GETVALUE fnCenterE = static_cast<amrex::Real (Pic::*)(amrex::MFIter&, amrex::IntVect, int, const int)>(&Pic::get_center_E);

  if (func == fnZero) {
    fMode = fModeZero;
  } else if (func == fnCenterB) {
    fMode = fModeCenterB;
  } else if (func == fnNodeB) {
    fMode = fModeNodeB;
  } else if (func == fnNodeE || func == fnCenterE) {
    fMode = fModeFluidE;
  }

  const bool hasCenterB = (fMode == fModeCenterB) && (fi != nullptr) &&
                          (!fi->get_center_b(iLev).empty());
  const bool hasNodeFluid = (fMode == fModeNodeB || fMode == fModeFluidE) &&
                            (fi != nullptr) && (!fi->get_node_fluid(iLev).empty());
  const int cOffset = (fMode == fModeCenterB) ? 0 :
                      (fMode == fModeNodeB) ? fi->get_iBx() :
                      (fMode == fModeFluidE) ? fi->get_iEx() : 0;

  if (bc != nullptr) {
    const BoxBC<FieldBC::Type> bcVal = *bc;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bxFab = mfi.fabbox();
      const Box& bxValid = mfi.validbox();

      if (!ba.contains(bxFab)) {
        Array4<Real> const& arr = mf[mfi].array();
        const Array4<const int>& statusArr = status[mfi].array();
        Array4<const Real> srcArr = hasCenterB ? fi->get_center_b(iLev)[mfi].array()
                                  : hasNodeFluid ? fi->get_node_fluid(iLev)[mfi].array()
                                  : Array4<const Real>{};

        ParallelFor(bxFab, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          if (bit::is_lev_boundary(statusArr(i, j, k, 0))) {
            int ip = i, jp = j, kp = k;
            bool useFloat = eval_use_float(i, j, k, ip, jp, kp, bcVal, bxValid, nDimLocal);

            if (useFloat) {
              for (int iVar = iStart; iVar < iStart + nComp; ++iVar) {
                arr(i, j, k, iVar) = arr(ip, jp, kp, iVar);
              }
            } else if (fMode == fModeZero) {
              for (int iVar = iStart; iVar < iStart + nComp; ++iVar) {
                arr(i, j, k, iVar) = 0.0;
              }
            } else if (hasCenterB || hasNodeFluid) {
              for (int iVar = iStart; iVar < iStart + nComp; ++iVar) {
                arr(i, j, k, iVar) = srcArr(i, j, k, cOffset + (iVar - iStart));
              }
            }
          }
        });
      }
    }
    return;
  }

  if (func == nullptr) {
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bxFab = mfi.fabbox();
      const Box& bxValid = mfi.validbox();

      if (!ba.contains(bxFab)) {
        const Array4<Real> arr = mf[mfi].array();
        const Array4<const int> statusArr = status[mfi].array();

        Box box = bxValid;
        box.grow(1);

        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          if (bit::is_lev_boundary(statusArr(i, j, k, 0))) {
            bool isNeiFound = false;
            const int kmin = (nDimLocal > 2) ? -1 : 0;
            const int kmax = (nDimLocal > 2) ? 1 : 0;
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
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.fabbox();

      if (!ba.contains(bx)) {
        Array4<Real> const& arr = mf[mfi].array();
        const Array4<const int>& statusArr = status[mfi].array();
        Array4<const Real> srcArr = hasCenterB ? fi->get_center_b(iLev)[mfi].array()
                                  : hasNodeFluid ? fi->get_node_fluid(iLev)[mfi].array()
                                  : Array4<const Real>{};

        auto lo = IntVect(bx.loVect());
        auto hi = IntVect(bx.hiVect());
        if (nDim > 2 && Geom(iLev).Domain().bigEnd(iz_) ==
                            Geom(iLev).Domain().smallEnd(iz_)) {
          lo[iz_]++;
          hi[iz_]--;
        }

        Box box0(lo, hi);

        ParallelFor(box0, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          if (bit::is_lev_boundary(statusArr(i, j, k, 0))) {
            if (fMode == fModeZero) {
              for (int iVar = 0; iVar < nComp; ++iVar) {
                arr(i, j, k, iStart + iVar) = 0.0;
              }
            } else if (hasCenterB || hasNodeFluid) {
              for (int iVar = 0; iVar < nComp; ++iVar) {
                arr(i, j, k, iStart + iVar) = srcArr(i, j, k, cOffset + iVar);
              }
            }
          }
        });
      }
    }
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
            arr(i, j, k, comp) = 0.0;
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
            arr(i, j, k, comp) = 0.0;
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
// Mirror ion moments into physical-wall ghost cells for smooth Ohm/Hall
// stencils.
void Pic::apply_centerPlasma_BC(const iMultiFab& status, MultiFab& mf,
                                const int iLev) {
  if (Geom(iLev).isAllPeriodic() || mf.nGrow() == 0)
    return;

  std::string nameFunc = "Pic::apply_centerPlasma_BC";
  timing_func(nameFunc);

  const BoxArray ba =
      get_boundary_active_ba(activeRegion, mf, Geom(iLev), nDim, iz_);
  const Dim3 domLo = Geom(iLev).Domain().smallEnd().dim3();
  const Dim3 domHi = Geom(iLev).Domain().bigEnd().dim3();
  const int nComp = mf.nComp();

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
      const int dLo[3] = { domLo.x, domLo.y, domLo.z };
      const int dHi[3] = { domHi.x, domHi.y, domHi.z };
      int m[3] = { i, j, k };
      bool touched = false;

      for (int d = 0; d < nDim; ++d) {
        if (ijk[d] < dLo[d]) {
          m[d] = 2 * dLo[d] - 1 - ijk[d];
          touched = true;
        } else if (ijk[d] > dHi[d]) {
          m[d] = 2 * dHi[d] + 1 - ijk[d];
          touched = true;
        }
      }

      if (!touched)
        return;
      for (int comp = 0; comp < nComp; ++comp) {
        arr(i, j, k, comp) = arr(m[0], m[1], m[2], comp);
      }
    });
  }
}

//==========================================================
void Pic::apply_wave_field(const iMultiFab& status, MultiFab& mf,
                           const int iStart, const int nComp, const int iLev,
                           const BoxBC<FieldBC::Type>& bc, int iField, Real t) {
  (void)bc;
  bool hasField = false;
  for (const auto& f : waveBC.faces) {
    for (const auto& c : f.comps) {
      if (c.iField == iField) {
        hasField = true;
        break;
      }
    }
    if (hasField)
      break;
  }
  if (!hasField)
    return;

  std::string nameFunc = "Pic::apply_wave_field";
  timing_func(nameFunc);

  const BoxArray ba =
      get_boundary_active_ba(activeRegion, mf, Geom(iLev), nDim, iz_);
  const BoundaryBounds bnd(Geom(iLev), mf.boxArray().ixType());

  const Real* plo = Geom(iLev).ProbLo();
  const Real* dx = Geom(iLev).CellSize();
  const Real offset[3] = { bnd.isNode[0] ? 0.0 : 0.5, bnd.isNode[1] ? 0.0 : 0.5,
                           bnd.isNode[2] ? 0.0 : 0.5 };

#if defined(AMREX_USE_GPU)
  MultiFab h_mf(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrowVect(),
                MFInfo().SetArena(The_Pinned_Arena()));
  h_mf.ParallelCopy(mf);
  const iMultiFab& h_status = mf.ixType().cellCentered() ? host_cell_status(iLev)
                                                         : host_node_status(iLev);
  Gpu::streamSynchronize();
#else
  MultiFab& h_mf = mf;
  const iMultiFab& h_status = status;
#endif

  for (MFIter mfi(h_mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    const Box& bxValid = mfi.validbox();
    if (ba.contains(bxFab))
      continue;

    Array4<Real> const& arr = h_mf[mfi].array();
    const Array4<const int>& statusArr = h_status[mfi].array();

    // Host-only kernel: evaluates waveBC host structures and std::vector faces
    amrex::LoopOnCpu(bxFab, [&](int i, int j, int k) {
      if (!bit::is_lev_boundary(statusArr(i, j, k, 0)))
        return;

      Real pos[3] = { plo[0] + dx[0] * (i + offset[0]),
                      plo[1] + dx[1] * (j + offset[1]),
                      (nDim > 2) ? plo[2] + dx[2] * (k + offset[2]) : 0.0 };

      Real waveVal[3] = { 0.0, 0.0, 0.0 };
      bool hasWave = false;

      for (const auto& f : waveBC.faces) {
        const int d = f.direction;
        const int side = f.side;
        if (d >= nDim)
          continue;
        const int idx = (d == 0) ? i : (d == 1) ? j : k;
        bool onFace = (side == 0 && idx < bxValid.smallEnd(d)) ||
                      (side == 1 && idx > bxValid.bigEnd(d));
        if (!onFace)
          continue;

        for (const auto& c : f.comps) {
          if (c.iField != iField)
            continue;
          hasWave = true;
          const Real val = waveBC.value(c, t, pos);
          if (iField == 0 || iField == 1) {
            for (int iVar = 0; iVar < std::min(nComp, 3); ++iVar) {
              waveVal[iVar] += val * c.pol[iVar];
            }
          } else {
            waveVal[0] += val;
          }
        }
      }

      if (hasWave) {
        if (iField == 0 || iField == 1) {
          for (int iVar = 0; iVar < nComp; ++iVar) {
            arr(i, j, k, iStart + iVar) = waveVal[iVar % 3];
          }
        } else if (nComp > 0) {
          arr(i, j, k, iStart) = waveVal[0];
        }
      }
    });
  }

#if defined(AMREX_USE_GPU)
  mf.ParallelCopy(h_mf);
  Gpu::streamSynchronize();
#endif
}

//==========================================================
void Pic::wave_velocity_kick(const Real* pos, Real t, Real& dvx, Real& dvy,
                             Real& dvz) {
  dvx = dvy = dvz = 0.0;
  if (!waveBC.active)
    return;
  for (const auto& f : waveBC.faces) {
    for (const auto& c : f.comps) {
      if (c.iField != 2) // velocity kick
        continue;
      const Real val = waveBC.value(c, t, pos);
      dvx += val * c.pol[0];
      dvy += val * c.pol[1];
      dvz += val * c.pol[2];
    }
  }
}

//==========================================================

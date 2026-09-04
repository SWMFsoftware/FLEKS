#include <algorithm>
#include <cmath>
#include <vector>

#include <AMReX_MultiFabUtil.H>

#include "GridUtility.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

//==========================================================
void Pic::apply_field_bc(const iMultiFab& status, MultiFab& mf,
                         const int iStart, const int nComp, GETVALUE func,
                         const int iLev, const bool isB) {
  std::string nameFunc = "Pic::apply_field_bc";
  timing_func(nameFunc);

  if (Geom(iLev).isAllPeriodic())
    return;
  if (mf.nGrow() == 0)
    return;

  // Base fill: float on open faces, or evaluate state from func elsewhere.
  apply_BC(status, mf, iStart, nComp, func, iLev, &bcField);

  // Dedicated wall operators applied per configured face type.
  apply_conducting_wall(status, mf, iStart, nComp, iLev, bcField, isB);
  apply_absorbing_wall(status, mf, iStart, nComp, iLev, bcField, isB);
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
  std::string nameFunc = "Pic::apply_BC";
  timing_func(nameFunc);

  if (Geom(iLev).isAllPeriodic())
    return;
  if (mf.nGrow() == 0)
    return;

  bool useFloatBC = (func == nullptr);
  BoxArray ba = convert(activeRegion, mf.boxArray().ixType());

  const IntVect& ngrow = mf.nGrowVect();
  if (nDim > 2 &&
      Geom(iLev).Domain().bigEnd(iz_) == Geom(iLev).Domain().smallEnd(iz_)) {
    ba.grow(iz_, ngrow[iz_]);
  }

  if (bc != nullptr) {
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bxFab = mfi.fabbox();
      const Box& bxValid = mfi.validbox();

      if (!ba.contains(bxFab)) {
        Array4<Real> const& arr = mf[mfi].array();
        const Array4<const int>& statusArr = status[mfi].array();

        ParallelFor(bxFab, [&](int i, int j, int k) {
          if (bit::is_lev_boundary(statusArr(i, j, k, 0))) {
            int ip, jp, kp;
            bool useFloat = use_float(i, j, k, ip, jp, kp, *bc, bxValid);

            if (useFloat) {
              for (int iVar = iStart; iVar < iStart + nComp; iVar++) {
                arr(i, j, k, iVar) = arr(ip, jp, kp, iVar);
              }
            } else if (func) {
              for (int iVar = iStart; iVar < iStart + nComp; iVar++) {
                arr(i, j, k, iVar) = (this->*func)(
                    mfi, IntVect{ AMREX_D_DECL(i, j, k) }, iVar - iStart, iLev);
              }
            }
          }
        });
      }
    }
    return;
  }

  if (useFloatBC) {
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
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.fabbox();

      if (!ba.contains(bx)) {
        Array4<Real> const& arr = mf[mfi].array();
        const Array4<const int>& statusArr = status[mfi].array();

        auto lo = IntVect(bx.loVect());
        auto hi = IntVect(bx.hiVect());
        if (nDim > 2 && Geom(iLev).Domain().bigEnd(iz_) ==
                            Geom(iLev).Domain().smallEnd(iz_)) {
          lo[iz_]++;
          hi[iz_]--;
        }

        Box box0(lo, hi);

        ParallelFor(box0, nComp, [&](int i, int j, int k, int iVar) {
          if (bit::is_lev_boundary(statusArr(i, j, k, 0))) {
            arr(i, j, k, iStart + iVar) = (this->*func)(
                mfi, IntVect{ AMREX_D_DECL(i, j, k) }, iVar, iLev);
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

  if (Geom(iLev).isAllPeriodic())
    return;
  if (mf.nGrow() == 0)
    return;

  bool hasConducting = false;
  for (int d = 0; d < nDim; ++d) {
    if (bc.lo[d] == FieldBC::conducting || bc.hi[d] == FieldBC::conducting) {
      hasConducting = true;
      break;
    }
  }
  if (!hasConducting)
    return;

  BoxArray ba = convert(activeRegion, mf.boxArray().ixType());
  const IntVect& ngrow = mf.nGrowVect();
  if (nDim > 2 &&
      Geom(iLev).Domain().bigEnd(iz_) == Geom(iLev).Domain().smallEnd(iz_)) {
    ba.grow(iz_, ngrow[iz_]);
  }

  const IntVect domLo = Geom(iLev).Domain().smallEnd();
  const IntVect domHi = Geom(iLev).Domain().bigEnd();
  const IndexType ixType = mf.boxArray().ixType();

  bool isNode[3] = { false, false, false };
  int loBnd[3] = { 0, 0, 0 };
  int hiBnd[3] = { 0, 0, 0 };
  for (int d = 0; d < nDim; ++d) {
    isNode[d] = (ixType[d] == IndexType::NODE);
    loBnd[d] = domLo[d];
    hiBnd[d] = isNode[d] ? (domHi[d] + 1) : domHi[d];
  }

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    if (ba.contains(bxFab))
      continue;

    Array4<Real> const& arr = mf[mfi].array();
    const Array4<const int>& statusArr = status[mfi].array();
    const Box& bxValid = mfi.validbox();

    ParallelFor(bxFab, [&](int i, int j, int k) {
      IntVect ijk{ AMREX_D_DECL(i, j, k) };

      // 1. Boundary nodes on physical wall (node-centred only).
      for (int d = 0; d < nDim; ++d) {
        if (!isNode[d])
          continue;

        const bool onLoWall =
            (bc.lo[d] == FieldBC::conducting) && (ijk[d] == loBnd[d]);
        const bool onHiWall =
            (bc.hi[d] == FieldBC::conducting) && (ijk[d] == hiBnd[d]);
        if (onLoWall || onHiWall) {
          bool inValid = true;
          for (int od = 0; od < nDim; ++od) {
            if (od != d && (ijk[od] < bxValid.smallEnd(od) ||
                            ijk[od] > bxValid.bigEnd(od))) {
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

      for (int d = 0; d < nDim; ++d) {
        bool isLow = (bc.lo[d] == FieldBC::conducting) && (ijk[d] < loBnd[d]);
        bool isHigh = (bc.hi[d] == FieldBC::conducting) && (ijk[d] > hiBnd[d]);
        if (!isLow && !isHigh)
          continue;

        IntVect m = ijk;
        if (isNode[d])
          m[d] = isLow ? (2 * loBnd[d] - ijk[d]) : (2 * hiBnd[d] - ijk[d]);
        else
          m[d] =
              isLow ? (2 * loBnd[d] - 1 - ijk[d]) : (2 * hiBnd[d] + 1 - ijk[d]);

        for (int iVar = 0; iVar < nComp; ++iVar) {
          const int comp = iStart + iVar;
          if (isB) {
            if (iVar == d)
              arr(i, j, k, comp) = 0.0;
            else
              arr(i, j, k, comp) = arr(m, comp);
          } else {
            if (iVar != d)
              arr(i, j, k, comp) = 0.0;
            else
              arr(i, j, k, comp) = arr(m, comp);
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
  std::string nameFunc = "Pic::apply_absorbing_wall";
  timing_func(nameFunc);

  (void)isB;
  if (Geom(iLev).isAllPeriodic())
    return;
  if (mf.nGrow() == 0)
    return;

  bool hasAbsorb = false;
  for (int d = 0; d < nDim; ++d) {
    if (bc.lo[d] == FieldBC::absorb || bc.hi[d] == FieldBC::absorb) {
      hasAbsorb = true;
      break;
    }
  }
  if (!hasAbsorb)
    return;

  const Real dt = tc ? tc->get_dt() : 0.0;
  if (dt <= 0.0)
    return;

  // Characteristic speed; default c=1, override via #ABSORB.
  const Real cs = (absorbCharSpeed > 0.0) ? absorbCharSpeed : 1.0;

  BoxArray ba = convert(activeRegion, mf.boxArray().ixType());
  const IntVect& ngrow = mf.nGrowVect();
  if (nDim > 2 &&
      Geom(iLev).Domain().bigEnd(iz_) == Geom(iLev).Domain().smallEnd(iz_)) {
    ba.grow(iz_, ngrow[iz_]);
  }

  const IntVect domLo = Geom(iLev).Domain().smallEnd();
  const IntVect domHi = Geom(iLev).Domain().bigEnd();
  const IndexType ixType = mf.boxArray().ixType();
  const Real* dx = Geom(iLev).CellSize();

  bool isNode[3] = { false, false, false };
  int loBnd[3] = { 0, 0, 0 };
  int hiBnd[3] = { 0, 0, 0 };
  Real decay[3] = { 0.0, 0.0, 0.0 };
  Real drive[3] = { 0.0, 0.0, 0.0 };

  for (int d = 0; d < nDim; ++d) {
    isNode[d] = (ixType[d] == IndexType::NODE);
    loBnd[d] = domLo[d];
    hiBnd[d] = isNode[d] ? (domHi[d] + 1) : domHi[d];

    const Real drive0 = cs * dt / dx[d];
    decay[d] = (1.0 - drive0) / (1.0 + drive0);
    drive[d] = 2.0 * drive0 / (1.0 + drive0);
  }

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    if (ba.contains(bxFab))
      continue;

    Array4<Real> const& arr = mf[mfi].array();
    const Array4<const int>& statusArr = status[mfi].array();

    ParallelFor(bxFab, [&](int i, int j, int k) {
      if (!bit::is_lev_boundary(statusArr(i, j, k, 0)))
        return;

      IntVect ijk{ AMREX_D_DECL(i, j, k) };

      for (int d = 0; d < nDim; ++d) {
        bool isLow = (bc.lo[d] == FieldBC::absorb) && (ijk[d] < loBnd[d]);
        bool isHigh = (bc.hi[d] == FieldBC::absorb) && (ijk[d] > hiBnd[d]);
        if (!isLow && !isHigh)
          continue;

        IntVect m = ijk;
        if (isNode[d])
          m[d] = isLow ? (2 * loBnd[d] - ijk[d]) : (2 * hiBnd[d] - ijk[d]);
        else
          m[d] =
              isLow ? (2 * loBnd[d] - 1 - ijk[d]) : (2 * hiBnd[d] + 1 - ijk[d]);

        for (int iVar = 0; iVar < nComp; ++iVar) {
          const int comp = iStart + iVar;
          arr(i, j, k, comp) =
              decay[d] * arr(i, j, k, comp) + drive[d] * arr(m, comp);
        }
      }
    });
  }
}

//==========================================================
void Pic::apply_inflow_wall(const iMultiFab& status, MultiFab& mf,
                            const int iStart, const int nComp, const int iLev,
                            const BoxBC<FieldBC::Type>& bc, bool isB) {
  std::string nameFunc = "Pic::apply_inflow_wall";
  timing_func(nameFunc);

  (void)isB; // zero-gradient copy is component-agnostic

  if (!fi->get_inflow_defined())
    return;
  if (Geom(iLev).isAllPeriodic())
    return;
  if (mf.nGrow() == 0)
    return;

  bool hasInflow = false;
  for (int d = 0; d < nDim; ++d) {
    if (bc.lo[d] == FieldBC::inflow || bc.hi[d] == FieldBC::inflow) {
      hasInflow = true;
      break;
    }
  }
  if (!hasInflow)
    return;

  BoxArray ba = convert(activeRegion, mf.boxArray().ixType());
  const IntVect& ngrow = mf.nGrowVect();
  if (nDim > 2 &&
      Geom(iLev).Domain().bigEnd(iz_) == Geom(iLev).Domain().smallEnd(iz_)) {
    ba.grow(iz_, ngrow[iz_]);
  }

  const IntVect domLo = Geom(iLev).Domain().smallEnd();
  const IntVect domHi = Geom(iLev).Domain().bigEnd();
  const IndexType ixType = mf.boxArray().ixType();

  bool isNode[3] = { false, false, false };
  int loBnd[3] = { 0, 0, 0 };
  int hiBnd[3] = { 0, 0, 0 };
  for (int d = 0; d < nDim; ++d) {
    isNode[d] = (ixType[d] == IndexType::NODE);
    loBnd[d] = domLo[d];
    hiBnd[d] = isNode[d] ? (domHi[d] + 1) : domHi[d];
  }

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    if (ba.contains(bxFab))
      continue;

    Array4<Real> const& arr = mf[mfi].array();
    const Array4<const int>& statusArr = status[mfi].array();

    ParallelFor(bxFab, [&](int i, int j, int k) {
      if (!bit::is_lev_boundary(statusArr(i, j, k, 0)))
        return;

      IntVect ijk{ AMREX_D_DECL(i, j, k) };

      for (int d = 0; d < nDim; ++d) {
        bool isLow = (bc.lo[d] == FieldBC::inflow) && (ijk[d] < loBnd[d]);
        bool isHigh = (bc.hi[d] == FieldBC::inflow) && (ijk[d] > hiBnd[d]);
        if (!isLow && !isHigh)
          continue;

        IntVect m = ijk;
        m[d] = isLow ? loBnd[d] : hiBnd[d];

        for (int iVar = 0; iVar < nComp; ++iVar) {
          arr(i, j, k, iStart + iVar) = arr(m, iStart + iVar);
        }
      }
    });
  }
}

//==========================================================
// Mirror ion moments into physical-wall ghost cells for smooth Ohm/Hall
// stencils.
void Pic::apply_centerPlasma_BC(const iMultiFab& status, MultiFab& mf,
                                const int iLev) {
  std::string nameFunc = "Pic::apply_centerPlasma_BC";
  timing_func(nameFunc);

  if (Geom(iLev).isAllPeriodic())
    return;
  if (mf.nGrow() == 0)
    return;

  const IntVect& ngrow = mf.nGrowVect();
  BoxArray ba = convert(activeRegion, mf.boxArray().ixType());
  if (nDim > 2 &&
      Geom(iLev).Domain().bigEnd(iz_) == Geom(iLev).Domain().smallEnd(iz_)) {
    ba.grow(iz_, ngrow[iz_]);
  }

  const IntVect domLo = Geom(iLev).Domain().smallEnd();
  const IntVect domHi = Geom(iLev).Domain().bigEnd();
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

      IntVect ijk{ AMREX_D_DECL(i, j, k) };
      IntVect m = ijk;
      bool touched = false;

      for (int d = 0; d < nDim; ++d) {
        if (ijk[d] < domLo[d]) {
          m[d] = 2 * domLo[d] - 1 - ijk[d];
          touched = true;
        } else if (ijk[d] > domHi[d]) {
          m[d] = 2 * domHi[d] + 1 - ijk[d];
          touched = true;
        }
      }

      if (!touched)
        return;
      for (int comp = 0; comp < nComp; ++comp) {
        arr(i, j, k, comp) = arr(m, comp);
      }
    });
  }
}

//==========================================================
void Pic::apply_wave_field(const iMultiFab& status, MultiFab& mf,
                           const int iStart, const int nComp, const int iLev,
                           const BoxBC<FieldBC::Type>& bc, int iField, Real t) {
  std::string nameFunc = "Pic::apply_wave_field";
  timing_func(nameFunc);

  (void)bc;
  if (!waveBC.active)
    return;
  if (Geom(iLev).isAllPeriodic())
    return;
  if (mf.nGrow() == 0)
    return;

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

  BoxArray ba = convert(activeRegion, mf.boxArray().ixType());
  const IntVect& ngrow = mf.nGrowVect();
  if (nDim > 2 &&
      Geom(iLev).Domain().bigEnd(iz_) == Geom(iLev).Domain().smallEnd(iz_)) {
    ba.grow(iz_, ngrow[iz_]);
  }

  const Real* plo = Geom(iLev).ProbLo();
  const Real* dx = Geom(iLev).CellSize();

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const Box& bxFab = mfi.fabbox();
    const Box& bxValid = mfi.validbox();
    if (ba.contains(bxFab))
      continue;

    Array4<Real> const& arr = mf[mfi].array();
    const Array4<const int>& statusArr = status[mfi].array();

    ParallelFor(bxFab, [&](int i, int j, int k) {
      if (!bit::is_lev_boundary(statusArr(i, j, k, 0)))
        return;

      Real pos[3] = { plo[0] + dx[0] * i, plo[1] + dx[1] * j,
                      plo[2] + dx[2] * k };

      Real waveVal[3] = { 0.0, 0.0, 0.0 };
      bool hasWave = false;

      for (const auto& f : waveBC.faces) {
        const int d = f.direction;
        const int side = f.side;
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

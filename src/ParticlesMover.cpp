#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>

#include <AMReX_ParReduce.H>

#include "InitialCondition.h"
#include "Morton.h"
#include "Particles.h"
#include "SWMFInterface.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

namespace {

// Charge-conserving (Esirkepov) current of one particle of charge q moving
// from x0 to x1 (at most one cell per dimension) with the node-CIC shape S,
// deposited onto the nodes of jArr (components x, y, z):
//   face current  F_d = -q dx_d / (dt V) cumsum_d W_d   (faces between nodes)
//   node current  J_d = (F_d(i - 1/2) + F_d(i + 1/2)) / 2 along d.
// With the cell-centered charge rho_c = (node-CIC charge averaged over the
// 2^D nodes of each cell), rho_c^{new} - rho_c^{old} = -dt
// div_node_to_center(J) holds to round-off. In 2D the out-of-plane component is
// q vOut / V times the time-averaged shape. Returns false if the particle moved
// more than one cell.
inline bool deposit_esirkepov(const Array4<Real>& jArr, const Real* x0,
                              const Real* x1, Real q, Real vOut,
                              const Real* plo, const Real* invDx,
                              const Real* dx, Real invDt, Real invVol) {
  // A particle that stays in its cell touches 2 nodes per dimension, one that
  // crosses a cell face 3.
  constexpr int nS = 3;
  int base[3] = { 0, 0, 0 };
  int n[3] = { 1, 1, 1 };
  Real S0[3][nS] = {}, DS[3][nS] = {};
  for (int d = 0; d < 3; ++d) {
    if (d >= nDim) {
      S0[d][0] = 1.0; // a single node with no motion
      continue;
    }
    const Real s0 = (x0[d] - plo[d]) * invDx[d];
    const Real s1 = (x1[d] - plo[d]) * invDx[d];
    const int i0 = static_cast<int>(std::floor(s0));
    const int i1 = static_cast<int>(std::floor(s1));
    if (std::abs(i1 - i0) > 1)
      return false;
    base[d] = std::min(i0, i1);
    n[d] = 2 + std::abs(i1 - i0);
    const Real f0 = s0 - i0, f1 = s1 - i1;
    Real S1[nS] = {};
    S0[d][i0 - base[d]] += 1.0 - f0;
    S0[d][i0 - base[d] + 1] += f0;
    S1[i1 - base[d]] += 1.0 - f1;
    S1[i1 - base[d] + 1] += f1;
    for (int a = 0; a < n[d]; ++a)
      DS[d][a] = S1[a] - S0[d][a];
  }

  // Shape weight transverse to direction d: product of the time-averaged
  // shapes of the other two directions (Esirkepov's 1/2 and 1/3 terms).
  auto transverse = [&](int e, int f, int ie, int jf) {
    return S0[e][ie] * S0[f][jf] + 0.5 * DS[e][ie] * S0[f][jf] +
           0.5 * S0[e][ie] * DS[f][jf] + DS[e][ie] * DS[f][jf] / 3.0;
  };

  for (int d = 0; d < nDim; ++d) {
    const int e = (d + 1) % 3, f = (d + 2) % 3;
    const Real coef = -q * dx[d] * invDt * invVol;
    for (int ie = 0; ie < n[e]; ++ie)
      for (int jf = 0; jf < n[f]; ++jf) {
        const Real wt = transverse(e, f, ie, jf);
        Real fPrev = 0.0, fCum = 0.0;
        for (int a = 0; a < n[d]; ++a) {
          fCum += coef * DS[d][a] * wt;
          const Real jNode = 0.5 * (fCum + fPrev);
          fPrev = fCum;
          if (jNode == 0.0)
            continue;
          int idx[3];
          idx[d] = base[d] + a;
          idx[e] = base[e] + ie;
          idx[f] = base[f] + jf;
          jArr(IntVect(AMREX_D_DECL(idx[0], idx[1], idx[2])), d) += jNode;
        }
      }
  }

  if (nDim == 2) {
    const Real coef = q * vOut * invVol;
    for (int a = 0; a < n[0]; ++a)
      for (int b = 0; b < n[1]; ++b) {
        const Real w = transverse(0, 1, a, b);
        if (w != 0.0)
          jArr(IntVect(AMREX_D_DECL(base[0] + a, base[1] + b, 0)), iz_) +=
              coef * w;
      }
  }
  return true;
}

// E at a particle with the shape of the charge-conserving current: along
// component d, the average of the two nodes of the particle's cell; across d
// (and for the out-of-plane component in 2D), linear. loIdx and dShift come
// from find_node_interpolation (dShift = weight of the upper node).
inline void gather_E_matched_to_current(const Array4<Real const>& EArr,
                                        const IntVect& loIdx,
                                        const RealVect& dShift, Real* ep) {
  for (int c = 0; c < 3; ++c) {
    Real w[3][2];
    for (int d = 0; d < 3; ++d) {
      if (d >= nDim) {
        w[d][0] = 1.0;
        w[d][1] = 0.0;
      } else if (d == c) {
        w[d][0] = 0.5;
        w[d][1] = 0.5;
      } else {
        w[d][0] = 1.0 - dShift[d];
        w[d][1] = dShift[d];
      }
    }
    Real sum = 0.0;
    for (int k = 0; k < (nDim > 2 ? 2 : 1); ++k)
      for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 2; ++i) {
          const IntVect idx(AMREX_D_DECL(loIdx[0] + i, loIdx[1] + j,
                                         nDim > 2 ? loIdx[nDim - 1] + k : 0));
          sum += w[0][i] * w[1][j] * w[2][k] * EArr(idx, c);
        }
    ep[c] = sum;
  }
}

// Collect boundary/body absorption tallies deterministically across OpenMP
// threads without allocating an O(nTiles) array: uses a single inline tally on
// the serial fast path, and records only active boundary tiles sorted by tile
// ordinal when running multithreaded.
template <class TallyT> class OrderedTallyCollector {
public:
  explicit OrderedTallyCollector(bool forceSingle = false)
      : nSlots_(forceSingle ? 1 : omp_slot_count()),
        perSlot_(nSlots_ > 1 ? nSlots_ : 0) {}

  bool is_single() const noexcept { return nSlots_ <= 1; }

  TallyT* single_ptr() noexcept { return &single_; }

  void record(int orderKey, const TallyT& t) {
    if (t.active)
      perSlot_[omp_slot_id()].emplace_back(orderKey, t);
  }

  template <class ApplyFn> void flush(ApplyFn&& apply) {
    if (nSlots_ <= 1) {
      apply(single_);
      return;
    }
    std::size_t total = 0;
    for (const auto& v : perSlot_)
      total += v.size();
    if (total == 0)
      return;
    std::vector<std::pair<int, TallyT> > merged;
    merged.reserve(total);
    for (auto& v : perSlot_)
      for (auto& kv : v)
        merged.push_back(std::move(kv));
    std::sort(merged.begin(), merged.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    for (const auto& kv : merged)
      apply(kv.second);
  }

private:
  int nSlots_ = 1;
  TallyT single_{};
  std::vector<std::vector<std::pair<int, TallyT> > > perSlot_;
};

} // namespace

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::update_position_to_half_stage(
    const MultiFab& /*nodeEMF*/, const MultiFab& /*nodeBMF*/, Real dt) {
  timing_func("Pts::update_position_to_half_stage");

  Real dtLoc = 0.5 * dt;

  const int iLev = 0;
  const LevelGeomBox lb = make_level_geom_box(iLev);
  const BoxArray& ba = this->ParticleBoxArray(iLev);
  const bool hasBodyQNet = grid && grid->has_body_net_charge(iLev);
  const bool hasBodyQNetInc = grid && grid->has_body_net_charge_inc(iLev);

  struct DriftFabCtx {
    ActiveRegionBox ab;
    Array4<int const> status;
    Array4<Real const> bodyQNet;
    Array4<Real> bodyQNetInc;
  };
  const auto fabTable =
      build_local_fab_table(cell_status(iLev), [&](const MFIter& mfi) {
        const Box& bx = cell_status(iLev)[mfi].box();
        return DriftFabCtx{ make_active_region_box(lb, ba[mfi.index()], bx),
                            cell_status(iLev)[mfi].const_array(),
                            hasBodyQNet
                                ? grid->body_net_charge(iLev)[mfi].const_array()
                                : Array4<Real const>(),
                            hasBodyQNetInc
                                ? grid->body_net_charge_inc(iLev)[mfi].array()
                                : Array4<Real>() };
      });

  OrderedTallyCollector<MoveTallyDelta> tallies(hasBodyQNetInc);

  parallel_for_cell_tiles(
      *this, iLev,
      [&](const auto& item, int ord) {
        AoS& particles = item.ptile->GetArrayOfStructs();
        const auto& ctx = fabTable[item.grid];
        MoveTallyDelta localTally;
        MoveTallyDelta* tally =
            tallies.is_single() ? tallies.single_ptr() : &localTally;

        for (auto& p : particles) {
          if (p.id() < 0)
            continue;

          for (int iDim = 0; iDim < nDim; iDim++) {
            p.pos(iDim) += p.rdata(iup_ + iDim) * dtLoc;
          }

          // Mark for deletion
          if (reflect_or_delete_particle(p, ctx.status, ctx.ab, ctx.bodyQNet,
                                         ctx.bodyQNetInc, tally)) {
            p.id() = -1;
          }
        } // for p

        if (!tallies.is_single())
          tallies.record(ord, localTally);
      },
      hasBodyQNetInc);

  tallies.flush([&](const MoveTallyDelta& d) { apply_move_tally_delta(d); });

  redistribute_particles();
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::mover(
    const Vector<MultiFab>& nodeE, const Vector<MultiFab>& nodeB,
    const Vector<MultiFab>& nodeB0, const Vector<MultiFab>& eBg,
    const Vector<MultiFab>& uBg, Real dt, Real dtNext, int pushStage,
    Vector<MultiFab>* currentMF, bool trialOnly) {
  if (is_neutral()) {
    if (currentMF) {
      for (auto& mf : *currentMF)
        mf.setVal(0.0);
    }
    if (!trialOnly)
      neutral_mover(dt);
  } else {
    charged_particle_mover(nodeE, nodeB, nodeB0, eBg, uBg, dt, dtNext,
                           pushStage, currentMF, trialOnly);
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::charged_particle_mover(
    const Vector<MultiFab>& nodeE, const Vector<MultiFab>& nodeB,
    const Vector<MultiFab>& nodeB0, const Vector<MultiFab>& /*eBg*/,
    const Vector<MultiFab>& /*uBg*/, Real dt, Real dtNext, int pushStage,
    Vector<MultiFab>* currentMF, bool trialOnly) {
  timing_func("Pts::charged_particle_mover");

  // nodeB0 is empty when no intrinsic field is configured.
  const bool useB0 = nodeB0.size() == nodeB.size() && !nodeB0[0].empty();

  const Real qdto2mc = charge / mass * 0.5 * dt;
  Real dtLoc = 0.5 * (dt + dtNext);

  // Current of the step (Esirkepov), deposited as the particles move.
  const Real qFactor = charge * qomSign;
  const Real invDt = 1.0 / dt;

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    const LevelGeomBox lb = make_level_geom_box(iLev);
    const BoxArray& ba = this->ParticleBoxArray(iLev);
    const bool hasBodyQNet = grid && grid->has_body_net_charge(iLev);
    const bool hasBodyQNetInc = grid && grid->has_body_net_charge_inc(iLev);

    if (currentMF)
      (*currentMF)[iLev].setVal(0.0);

    struct MoverFabCtx {
      ActiveRegionBox ab;
      Array4<Real> jArr;
      Array4<Real const> bodyQNet;
      Array4<Real> bodyQNetInc;
      Array4<Real const> EArr;
      Array4<Real const> BArr;
      Array4<Real const> B0Arr;
      Array4<int const> status;
    };
    const auto fabTable =
        build_local_fab_table(cell_status(iLev), [&](const MFIter& mfi) {
          const Box& bx = cell_status(iLev)[mfi].box();
          return MoverFabCtx{
            make_active_region_box(lb, ba[mfi.index()], bx),
            currentMF ? (*currentMF)[iLev][mfi].array() : Array4<Real>(),
            hasBodyQNet ? grid->body_net_charge(iLev)[mfi].const_array()
                        : Array4<Real const>(),
            hasBodyQNetInc ? grid->body_net_charge_inc(iLev)[mfi].array()
                           : Array4<Real>(),
            nodeE[iLev][mfi].const_array(),
            nodeB[iLev][mfi].const_array(),
            useB0 ? nodeB0[iLev][mfi].const_array() : Array4<Real const>(),
            cell_status(iLev)[mfi].const_array()
          };
        });

    OrderedTallyCollector<MoveTallyDelta> tallies(hasBodyQNetInc);
    std::atomic<bool> movedTooFar{ false };

    auto processTile = [&](const auto& item, int ord, auto checkHomeTag,
                           bool strayPass) -> bool {
      AoS& particles = item.ptile->GetArrayOfStructs();
      const auto& ctx = fabTable[item.grid];
      const Real ploTile[3] = { lb.plo[0], lb.plo[1], lb.plo[2] };
      const Real invDxTile[3] = { lb.invDx[0], lb.invDx[1], lb.invDx[2] };
      const Real dxTile[3] = { lb.dx[0], lb.dx[1], lb.dx[2] };
      const Real qdto2mcTile = qdto2mc;
      const Real dtTile = dt;
      const Real dtLocTile = dtLoc;
      const Real invDtTile = invDt;
      const Real qFactorTile = qFactor;
      const Real invVolTile = invVol[iLev];
      const int stageTile = pushStage;
      const bool useB0Tile = useB0;
      const bool hasCurrentTile = (currentMF != nullptr);
      const bool trialOnlyTile = trialOnly;
      const int iup = iup_;
      const int ivp = ivp_;
      const int iwp = iwp_;
      const int iqp = iqp_;
      const int imu = imu_;
      const bool saveMu = (pMode == PartMode::PIC && imu < NStructReal);
      const bool singleTally = tallies.is_single();
      const Array4<Real const> fields2[2] = { ctx.BArr, ctx.EArr };
      const Array4<Real const> fields3[3] = { ctx.BArr, ctx.EArr, ctx.B0Arr };
      const Dim3 lo = init_dim3(0);
      const Dim3 hi = init_dim3(1);
      bool foundStray = false;
      MoveTallyDelta localTally;
      MoveTallyDelta* tally = singleTally ? tallies.single_ptr() : &localTally;

      for (auto& p : particles) {
        if (p.id() < 0)
          continue;

        if constexpr (decltype(checkHomeTag)::value) {
          bool isHome = true;
          for (int d = 0; d < nDim; ++d) {
            const int i0 = static_cast<int>(
                std::floor((p.pos(d) - ploTile[d]) * invDxTile[d]));
            if (i0 != item.cell[d]) {
              isHome = false;
              break;
            }
          }
          if (!strayPass) {
            if (!isHome) {
              foundStray = true;
              continue;
            }
          } else if (isHome) {
            continue;
          }
        }

        const Real u_start = p.rdata(iup);
        const Real v_start = p.rdata(ivp);
        const Real w_start = p.rdata(iwp);
        Real up = u_start;
        Real vp = v_start;
        Real wp = w_start;
        const Real xp = p.pos(ix_);
        const Real yp = p.pos(iy_);
        const Real zp = nDim > 2 ? p.pos(iz_) : 0;

        //-----calculate interpolate coef begin-------------
        // The stencil is centered on the node that contains the particle: plain
        // trilinear weights over the surrounding 2x2x2 nodes.
        IntVect loIdx;
        RealVect dShift;

        Real coef[2][2][2];
        RealVect posEval = p.pos();
        if (stageTile == 2) {
          posEval[ix_] += 0.5 * u_start * dtTile;
          posEval[iy_] += 0.5 * v_start * dtTile;
          if (nDim > 2)
            posEval[iz_] += 0.5 * w_start * dtTile;
        }
        find_node_interpolation(posEval, ploTile, invDxTile, loIdx, dShift,
                                coef);
        //-----calculate interpolate coef end-------------

        Real bp[3] = { 0, 0, 0 };
        Real ep[3] = { 0, 0, 0 };
        Real u0p[3] = { 0, 0, 0 };
        Real b0p[3] = { 0, 0, 0 };
        if (useB0Tile) {
          Real* values[3] = { bp, ep, b0p };
          interpolate_vector_fields(fields3, loIdx, coef, lo, hi, values);
          bp[ix_] += b0p[ix_];
          bp[iy_] += b0p[iy_];
          bp[iz_] += b0p[iz_];
        } else {
          Real* values[2] = { bp, ep };
          interpolate_vector_fields(fields2, loIdx, coef, lo, hi, values);
        }

        // The charge-conserving current of E_d has a zeroth-order shape along
        // d (the average of the two nodes of the cell the particle is in).
        // Gather E with the same shape, the transpose of the deposit, so the
        // work of the field on the particles matches the field energy loss.
        if (hasCurrentTile)
          gather_E_matched_to_current(ctx.EArr, loIdx, dShift, ep);

        up = up - u0p[ix_];
        vp = vp - u0p[iy_];
        wp = wp - u0p[iz_];

        const Real omx = qdto2mcTile * bp[ix_];
        const Real omy = qdto2mcTile * bp[iy_];
        const Real omz = qdto2mcTile * bp[iz_];

        // end interpolation
        const Real velocity[3] = { up, vp, wp };
        const Real ut = up + qdto2mcTile * ep[ix_];
        const Real vt = vp + qdto2mcTile * ep[iy_];
        const Real wt = wp + qdto2mcTile * ep[iz_];
        const Real electricVelocity[3] = { ut, vt, wt };
        const Real omega[3] = { omx, omy, omz };
        Real updatedVelocity[3];
        boris_push_nonrelativistic(velocity, electricVelocity, omega,
                                   updatedVelocity);

        const Real unp1 = updatedVelocity[ix_] + u0p[ix_];
        const Real vnp1 = updatedVelocity[iy_] + u0p[iy_];
        const Real wnp1 = updatedVelocity[iz_] + u0p[iz_];

        // Drift velocity: v^{n+1/2} (leapfrog) or the average of the start
        // and end velocities (predictor-corrector stages).
        const Real uDrift = stageTile > 0 ? 0.5 * (u_start + unp1) : unp1;
        const Real vDrift = stageTile > 0 ? 0.5 * (v_start + vnp1) : vnp1;
        const Real wDrift = stageTile > 0 ? 0.5 * (w_start + wnp1) : wnp1;
        const Real xOld[3] = { xp, yp, zp };
        const Real xNew[3] = { xp + uDrift * dtLocTile, yp + vDrift * dtLocTile,
                               zp + wDrift * dtLocTile };

        // J = q dx/dt for every component; in 2D the out-of-plane dx is
        // wDrift * dtLoc (dtLoc != dt when the time step changes).
        if (hasCurrentTile &&
            !deposit_esirkepov(ctx.jArr, xOld, xNew, qFactorTile * p.rdata(iqp),
                               wDrift * dtLocTile * invDtTile, ploTile,
                               invDxTile, dxTile, invDtTile, invVolTile))
          movedTooFar.store(true, std::memory_order_relaxed);

        // A trial push only deposits the current of the trial state.
        if (trialOnlyTile)
          continue;

        p.rdata(iup) = unp1;
        p.rdata(ivp) = vnp1;
        p.rdata(iwp) = wnp1;

        if (saveMu) {
          // Note: bp should be calculated at the new position. Now, bp at the
          // old position is used to save the calculation.
          p.rdata(imu) = cosine(p, bp);
        }

        p.pos(ix_) = xNew[ix_];
        p.pos(iy_) = xNew[iy_];
        if (nDim > 2)
          p.pos(iz_) = xNew[iz_];

        // Apply boundary condition (absorb: delete; reflect: mirror).
        if (reflect_or_delete_particle(p, ctx.status, ctx.ab, ctx.bodyQNet,
                                       ctx.bodyQNetInc, tally)) {
          p.id() = -1;
        }
      } // for p

      if (!singleTally)
        tallies.record(ord * 2 + (strayPass ? 1 : 0), localTally);

      return foundStray;
    };

    deposit_cell_tiles_colored(*this, iLev, currentMF ? 4 : 0, processTile,
                               hasBodyQNetInc);

    tallies.flush([&](const MoveTallyDelta& d) { apply_move_tally_delta(d); });

    if (movedTooFar.load(std::memory_order_relaxed))
      Abort("charged_particle_mover: a particle moved more than one cell in a "
            "step; the charge-conserving current needs a smaller time step.");

    if (currentMF) {
      fold_wall_ghost_nodes((*currentMF)[iLev], iLev, true);
      (*currentMF)[iLev].SumBoundary(Geom(iLev).periodicity());
    }
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::neutral_mover(Real dt) {
  timing_func("Pts::neutral_mover");

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    const LevelGeomBox lb = make_level_geom_box(iLev);
    const BoxArray& ba = this->ParticleBoxArray(iLev);
    const bool hasBodyQNet = grid && grid->has_body_net_charge(iLev);
    const bool hasBodyQNetInc = grid && grid->has_body_net_charge_inc(iLev);

    struct NeutralFabCtx {
      ActiveRegionBox ab;
      Array4<int const> status;
      Array4<Real const> bodyQNet;
      Array4<Real> bodyQNetInc;
    };
    const auto fabTable =
        build_local_fab_table(cell_status(iLev), [&](const MFIter& mfi) {
          const Box& bx = cell_status(iLev)[mfi].box();
          return NeutralFabCtx{
            make_active_region_box(lb, ba[mfi.index()], bx),
            cell_status(iLev)[mfi].const_array(),
            hasBodyQNet ? grid->body_net_charge(iLev)[mfi].const_array()
                        : Array4<Real const>(),
            hasBodyQNetInc ? grid->body_net_charge_inc(iLev)[mfi].array()
                           : Array4<Real>()
          };
        });

    OrderedTallyCollector<MoveTallyDelta> tallies(hasBodyQNetInc);

    parallel_for_cell_tiles(
        *this, iLev,
        [&](const auto& item, int ord) {
          AoS& particles = item.ptile->GetArrayOfStructs();
          const auto& ctx = fabTable[item.grid];
          MoveTallyDelta localTally;
          MoveTallyDelta* tally =
              tallies.is_single() ? tallies.single_ptr() : &localTally;

          for (auto& p : particles) {
            if (p.id() < 0)
              continue;

            const Real up = p.rdata(iup_);
            const Real vp = p.rdata(ivp_);
            const Real wp = p.rdata(iwp_);
            const Real xp = p.pos(ix_);
            const Real yp = p.pos(iy_);
            const Real zp = p.pos(iz_);

            p.pos(ix_) = xp + up * dt;
            p.pos(iy_) = yp + vp * dt;
            p.pos(iz_) = zp + wp * dt;

            // Apply boundary condition (absorb: delete; reflect: mirror).
            if (reflect_or_delete_particle(p, ctx.status, ctx.ab, ctx.bodyQNet,
                                           ctx.bodyQNetInc, tally)) {
              p.id() = -1;
            }
          } // for p

          if (!tallies.is_single())
            tallies.record(ord, localTally);
        },
        hasBodyQNetInc);

    tallies.flush([&](const MoveTallyDelta& d) { apply_move_tally_delta(d); });
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::divE_correct_position(
    const amrex::Vector<MultiFab>& phiMF, int iLev) {
  timing_func("Pts:divE_correct_position");

  const Real sign = charge / fabs(charge);
  const Real epsLimit = 0.1;

  const Real* const ploLoc = plo[iLev].begin();
  const Real* const invDxLoc = invDx[iLev].begin();
  const Real dx0 = dx[iLev][ix_];
  const Real dx1 = dx[iLev][iy_];
  const Real dx2 = (nDim > 2) ? dx[iLev][iz_] : 0.0;
  const Real invV = invVol[iLev];
#if AMREX_SPACEDIM > 2
  const int nD = get_dim();
#endif
  const bool useBody = (grid != nullptr && grid->use_body());
  const LevelGeomBox lb = make_level_geom_box(iLev);
  const BoxArray& ba = this->ParticleBoxArray(iLev);

  struct DivEFabCtx {
    ActiveRegionBox ab;
    Array4<Real const> phiArr;
    Array4<int const> status;
  };
  const auto fabTable =
      build_local_fab_table(cell_status(iLev), [&](const MFIter& mfi) {
        const Box& bx = cell_status(iLev)[mfi].box();
        return DivEFabCtx{ make_active_region_box(lb, ba[mfi.index()], bx),
                           phiMF[iLev][mfi].const_array(),
                           cell_status(iLev)[mfi].const_array() };
      });

  parallel_for_cell_tiles(*this, iLev, [&](const auto& item, int /*ord*/) {
    const auto& ctx = fabTable[item.grid];
    const Array4<Real const>& phiArr = ctx.phiArr;
    const Array4<int const>& status = ctx.status;
    const ActiveRegionBox& ab = ctx.ab;

    AoS& particles = item.ptile->GetArrayOfStructs();

    for (auto& p : particles) {
      if (p.id() == -1 || is_outside_active_region(p, status, ab)) {
        p.id() = -1;
        continue;
      }

      if (skip_particle_for_dive_cleaning(p.pos(), Geom(iLev), iLev, status) &&
          n_lev() > 1) {
        continue;
      }

      IntVect loIdx;
      RealVect dShift;
      find_cell_index(p.pos(), ploLoc, invDxLoc, loIdx, dShift);

      // Since the boundary condition for solving phi is not perfect,
      // correcting particles that are close to the boundaries may produce
      // artificial oscillations, which are seen in Earth's magnetotail
      // simulations. So, it is better to skip the boundary physical cells.
      if (iLev == 0) {
        bool isBoundaryPhysicalCell = false;
        for (int iz = 0; iz <= 1; iz++)
          for (int iy = 0; iy <= 1; iy++)
            for (int ix = 0; ix <= 1; ix++) {
              const int st = status(loIdx[ix_] + ix, loIdx[iy_] + iy,
#if (AMREX_SPACEDIM > 2)
                                    loIdx[iz_] + iz
#else
                                    0
#endif
              );
              if (bit::is_lev_boundary(st) || (useBody && bit::is_body(st)))
                isBoundaryPhysicalCell = true;
            }
        if (isBoundaryPhysicalCell)
          continue;
      }

      {
        //----- Mass matrix calculation begin--------------
        const Real xi0 = dShift[ix_] * dx0;
        const Real eta0 = dShift[iy_] * dx1;
        const Real zeta0 = nDim > 2 ? dShift[iz_] * dx2 : 0;
        const Real xi1 = dx0 - xi0;
        const Real eta1 = dx1 - eta0;
        const Real zeta1 = nDim > 2 ? dx2 - zeta0 : 1;

        const Real zeta12Vol = zeta1 * invV;

        RealVect eps_D = { AMREX_D_DECL(0, 0, 0) };
#if AMREX_SPACEDIM > 2
        const Real zeta02Vol = zeta0 * invV;
        const Real eta02Vol = eta0 * invV;
        const Real eta12Vol = eta1 * invV;
        alignas(32) const Real wx[8] = { -eta1 * zeta12Vol, eta1 * zeta12Vol,
                                         -eta0 * zeta12Vol, eta0 * zeta12Vol,
                                         -eta1 * zeta02Vol, eta1 * zeta02Vol,
                                         -eta0 * zeta02Vol, eta0 * zeta02Vol };
        alignas(32) const Real wy[8] = { -xi1 * zeta12Vol, -xi0 * zeta12Vol,
                                         xi1 * zeta12Vol,  xi0 * zeta12Vol,
                                         -xi1 * zeta02Vol, -xi0 * zeta02Vol,
                                         xi1 * zeta02Vol,  xi0 * zeta02Vol };
        alignas(32) const Real wz[8] = { -xi1 * eta12Vol, -xi0 * eta12Vol,
                                         -xi1 * eta02Vol, -xi0 * eta02Vol,
                                         xi1 * eta12Vol,  xi0 * eta12Vol,
                                         xi1 * eta02Vol,  xi0 * eta02Vol };
        const int i0 = loIdx[ix_], j0 = loIdx[iy_], k0 = loIdx[iz_];
        const Real* AMREX_RESTRICT p00 = phiArr.ptr(i0, j0, k0);
        const Real* AMREX_RESTRICT p10 = phiArr.ptr(i0, j0 + 1, k0);
        const Real* AMREX_RESTRICT p01 = phiArr.ptr(i0, j0, k0 + 1);
        const Real* AMREX_RESTRICT p11 = phiArr.ptr(i0, j0 + 1, k0 + 1);
        alignas(32) const Real phi8[8] = { p00[0], p00[1], p10[0], p10[1],
                                           p01[0], p01[1], p11[0], p11[1] };
        Real ex = 0.0, ey = 0.0, ez = 0.0;
#pragma omp simd reduction(+ : ex, ey, ez)
        for (int n = 0; n < 8; ++n) {
          ex += phi8[n] * wx[n];
          ey += phi8[n] * wy[n];
          ez += phi8[n] * wz[n];
        }
        eps_D[ix_] = ex;
        eps_D[iy_] = ey;
        eps_D[iz_] = (nD > 2) ? ez : 0.0;
#else
        alignas(32) const Real wx[4] = { -eta1 * zeta12Vol, eta1 * zeta12Vol,
                                         -eta0 * zeta12Vol, eta0 * zeta12Vol };
        alignas(32) const Real wy[4] = { -xi1 * zeta12Vol, -xi0 * zeta12Vol,
                                         xi1 * zeta12Vol, xi0 * zeta12Vol };
        const int i0 = loIdx[ix_], j0 = loIdx[iy_];
        const Real* AMREX_RESTRICT p00 = phiArr.ptr(i0, j0, 0);
        const Real* AMREX_RESTRICT p10 = phiArr.ptr(i0, j0 + 1, 0);
        alignas(32) const Real phi4[4] = { p00[0], p00[1], p10[0], p10[1] };
        Real ex = 0.0, ey = 0.0;
#pragma omp simd reduction(+ : ex, ey)
        for (int n = 0; n < 4; ++n) {
          ex += phi4[n] * wx[n];
          ey += phi4[n] * wy[n];
        }
        eps_D[ix_] = ex;
        eps_D[iy_] = ey;
#endif

        for (int iDim = 0; iDim < nDim; iDim++)
          eps_D[iDim] *= sign * fourPI;

        Real eps_D_dot_invDx_Max = 0.0;
        for (int iDim = 0; iDim < nDim; iDim++) {
          if (fabs(eps_D[iDim] * invDx[iLev][iDim]) > eps_D_dot_invDx_Max) {
            eps_D_dot_invDx_Max = fabs(eps_D[iDim] * invDx[iLev][iDim]);
          }
        }

        if (eps_D_dot_invDx_Max > epsLimit) {
          // If eps_D is too large, the underlying assumption of the particle
          // correction method will be not valid. Comparing each exp_D
          // component instead of the length dl saves the computational time.
          const Real dl = eps_D.vectorLength();
          const Real ratio = epsLimit * dx[iLev][ix_] / dl;
          for (int iDim = 0; iDim < nDim; iDim++)
            eps_D[iDim] *= ratio;
        }

        for (int iDim = 0; iDim < nDim; iDim++) {
          p.pos(iDim) += eps_D[iDim];
        }

        if (is_outside_active_region(p, status, ab)) {
          // Do not allow moving particles from physical cells to ghost cells
          // during divE correction.
          for (int iDim = 0; iDim < nDim; iDim++) {
            p.pos(iDim) -= eps_D[iDim];
          }

          // p.id() = -1;
        }
      }
    } // for p
  });
}

//==========================================================

#define INSTANTIATE_PARTICLES_MOVER(T)                                         \
  template void T::update_position_to_half_stage(const MultiFab&,              \
                                                 const MultiFab&, Real);       \
  template void T::mover(const Vector<MultiFab>&, const Vector<MultiFab>&,     \
                         const Vector<MultiFab>&, const Vector<MultiFab>&,     \
                         const Vector<MultiFab>&, Real, Real, int,             \
                         Vector<MultiFab>*, bool);                             \
  template void T::charged_particle_mover(                                     \
      const Vector<MultiFab>&, const Vector<MultiFab>&,                        \
      const Vector<MultiFab>&, const Vector<MultiFab>&,                        \
      const Vector<MultiFab>&, Real, Real, int, Vector<MultiFab>*, bool);      \
  template void T::neutral_mover(Real);                                        \
  template void T::divE_correct_position(const Vector<MultiFab>&, int);

INSTANTIATE_PARTICLES_MOVER(PicParticles)
INSTANTIATE_PARTICLES_MOVER(PTParticles)
#undef INSTANTIATE_PARTICLES_MOVER

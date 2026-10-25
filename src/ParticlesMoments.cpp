#include <cstdlib>

#include <AMReX_ParReduce.H>
#include <AMReX_ParticleInterpolators.H>
#include <AMReX_ParticleReduce.H>

#include "GpuBackend.h"
#include "InitialCondition.h"
#include "Morton.h"
#include "Particles.h"
#include "SWMFInterface.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void deposit_charge(
    Array4<Real> const& chargeArr, const IntVect& loIdx,
    const Real coef[2][2][2], Real charge, const Dim3& lo, const Dim3& hi) {
  for (int kk = lo.z; kk <= hi.z; ++kk)
    for (int jj = lo.y; jj <= hi.y; ++jj)
      for (int ii = lo.x; ii <= hi.x; ++ii) {
        const IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + ii, loIdx[iy_] + jj,
                                           loIdx[iz_] + kk) };
        deposit_add(chargeArr.ptr(ijk, 0), coef[ii][jj][kk] * charge);
      }
}

// Deterministic GPU deposit of the cell-centred mass matrix of the divE
// correction (see accumulate_mass_matrix_terms). The entry of the pair
// (n1, n2) of corners of a dual cell is
//   |q| invVol^3 sum_d sgn(n1_d) sgn(n2_d) prod_{e != d} h_e^2 P_e(s_e),
// where s is the position of the particle in the dual cell (0..1), h the cell
// size, sgn(1) = 1, sgn(0) = -1, and P_e is (1-s)^2, s(1-s) or s^2 for the
// corner pair (0,0), (0,1)/(1,0), (1,1). So a dual cell only needs the 27
// moments sum_p |q| sx^l sy^m sz^n (l, m, n = 0..2): the particles of a cell
// are summed in a fixed order (CellBins), then every cell centre gathers its
// entries from the dual cells around it. This is done for all the tiles of
// lb at once (see LevelBins): stateOf(t, ip, dShift, q) returns the shift and
// the charge of particle ip of tile t, and mmOf(t) the Array4 of tile t.
template <class StateOf, class MMOf>
void deposit_level_center_mm(const LevelBins& lb, const Real* dx, Real invVol,
                             StateOf stateOf, MMOf mmOf) {
  if (lb.n_particles() == 0)
    return;
  const int nBins = lb.n_bins(), nTile = lb.n_tiles();
  const DepositTile* const tl = lb.tiles();
  auto tile_of_bin = [=] AMREX_GPU_DEVICE(int b) {
    return find_deposit_tile(tl, nTile, b,
                             [](const DepositTile& d) { return d.bStart; });
  };

  // 1. The 27 moments of every dual cell.
  constexpr int nM = 27;
  Gpu::DeviceVector<Real> scratch(static_cast<std::size_t>(nBins) * nM);
  Real* const sc = scratch.data();
  const int* const perm = lb.perm();
  const int* const off = lb.offsets();
  // Three parts per cell, one per power n of sz, each summed by a warp (see
  // warp_sums).
  warp_sums<9>(
      nBins * 3,
      [=] AMREX_GPU_DEVICE(int t, int& begin, int& end) {
        const int b = t / 3;
        begin = off[b];
        end = off[b + 1];
      },
      [=] AMREX_GPU_DEVICE(int t, int i, Real* acc) {
        const int n = t % 3;
        const int iTile = tile_of_bin(t / 3);
        RealVect dShift;
        Real q;
        stateOf(iTile, perm[i] - tl[iTile].pStart, dShift, q);
        const Real w = amrex::Math::abs(q);
        const Real sz = nDim > 2 ? dShift[nDim > 2 ? 2 : 0] : 0;
        const Real px[3] = { 1, dShift[0], dShift[0] * dShift[0] };
        const Real py[3] = { 1, dShift[1], dShift[1] * dShift[1] };
        const Real pz[3] = { 1, sz, sz * sz };
        for (int l = 0; l < 3; ++l) {
          const Real wln = w * px[l] * pz[n];
          for (int k = 0; k < 3; ++k)
            acc[k * 3 + l] += wln * py[k];
        }
      },
      [=] AMREX_GPU_DEVICE(int t, int kl, Real v) {
        sc[static_cast<std::size_t>(t / 3) * nM + (t % 3) * 9 + kl] = v;
      });

  // 2. Every cell centre gathers the entries from the dual cells around it.
  const Real h2[3] = { dx[0] * dx[0], dx[1] * dx[1],
                       nDim > 2 ? dx[nDim > 2 ? 2 : 0] * dx[nDim > 2 ? 2 : 0]
                                : 1 };
  const Real coef = invVol * invVol * invVol;
  const int nK = nDim > 2 ? 2 : 1;
  // The node boxes of lb are the cell boxes grown by one on the high side.
  ParallelFor(lb.n_nodes(), [=] AMREX_GPU_DEVICE(int gNode) {
    const int iTile = find_deposit_tile(
        tl, nTile, gNode, [](const DepositTile& d) { return d.nodeStart; });
    const DepositTile& dt = tl[iTile];
    const int mNode = gNode - dt.nodeStart;
    const int k = dt.nodeLo.z + mNode / (dt.nodeLen.x * dt.nodeLen.y);
    const int j = dt.nodeLo.y + (mNode / dt.nodeLen.x) % dt.nodeLen.y;
    const int i = dt.nodeLo.x + mNode % dt.nodeLen.x;
    const int node[3] = { i, j, k };
    const int loC[3] = { dt.lo.x, dt.lo.y, dt.lo.z };
    const int lenC[3] = { dt.len.x, dt.len.y, dt.len.z };
    // Polynomial of the corner pair (a, b) in s (coefficients of 1, s, s^2).
    auto poly = [](int a, int b, int m) -> Real {
      constexpr Real p00[3] = { 1, -2, 1 }, p01[3] = { 0, 1, -1 },
                     p11[3] = { 0, 0, 1 };
      return (a == 0 && b == 0) ? p00[m] : (a != b ? p01[m] : p11[m]);
    };
    Real sum[nCMMComponents];
    for (int g = 0; g < nCMMComponents; ++g)
      sum[g] = 0;
    // The dual cells that have this centre as corner a.
    for (int cz = 0; cz < nK; ++cz)
      for (int cy = 0; cy < 2; ++cy)
        for (int cx = 0; cx < 2; ++cx) {
          const int a[3] = { cx, cy, cz };
          int c[3] = { 0, 0, 0 };
          bool inside = true;
          for (int d = 0; d < nDim; ++d) {
            c[d] = node[d] - a[d] - loC[d];
            inside = inside && c[d] >= 0 && c[d] < lenC[d];
          }
          if (!inside)
            continue;
          const Real* const m =
              sc + static_cast<std::size_t>(dt.bStart + c[0] +
                                            lenC[0] * (c[1] + lenC[1] * c[2])) *
                       nM;
          // The other corner b of the pair.
          for (int bz = 0; bz < nK; ++bz)
            for (int by = 0; by < 2; ++by)
              for (int bx = 0; bx < 2; ++bx) {
                const int bb[3] = { bx, by, bz };
                Real pe[3][3];
                for (int e = 0; e < 3; ++e)
                  for (int o = 0; o < 3; ++o)
                    pe[e][o] = poly(a[e], bb[e], o);
                Real value = 0;
                for (int d = 0; d < nDim; ++d) {
                  const Real sgn = (a[d] == bb[d]) ? 1 : -1;
                  Real hProd = 1;
                  for (int e = 0; e < 3; ++e)
                    if (e != d)
                      hProd *= h2[e];
                  Real t = 0;
                  for (int n = 0; n < 3; ++n) {
                    const Real fz = (d == 2) ? (n == 0) : pe[2][n];
                    if (fz == 0)
                      continue;
                    for (int kk = 0; kk < 3; ++kk) {
                      const Real fy = (d == 1) ? (kk == 0) : pe[1][kk];
                      if (fy == 0)
                        continue;
                      for (int l = 0; l < 3; ++l) {
                        const Real fx = (d == 0) ? (l == 0) : pe[0][l];
                        if (fx == 0)
                          continue;
                        t += fx * fy * fz * m[(n * 3 + kk) * 3 + l];
                      }
                    }
                  }
                  value += sgn * hProd * t;
                }
                const int ip = bx - cx + 1, jp = by - cy + 1;
                const int gp =
                    nDim > 2 ? (ip * 3 + jp) * 3 + (bz - cz + 1) : ip * 3 + jp;
                sum[gp] += coef * value;
              }
        }
    RealCMM& data = mmOf(iTile)(i, j, k);
    for (int g = 0; g < nCMMComponents; ++g)
      data[g] += sum[g];
  });
  // The kernels read the bins and the scratch, which are freed on return.
  sync_before_free();
}

// deposit_level_center_mm() for one tile: the cells of cellBox, cellOf and
// stateOf(ip, dShift, q) for its particles, as in deposit_to_nodes_with().
template <class CellOf, class StateOf>
void deposit_center_mm_by_cell(int np, const Box& cellBox,
                               const Array4<RealCMM>& mmArr, const Real* dx,
                               Real invVol, CellOf cellOf, StateOf stateOf) {
  if (np == 0)
    return;
  LevelBins lb;
  lb.build({ np }, { cellBox },
           [=] AMREX_GPU_DEVICE(int, int ip, IntVect& cell) {
             return cellOf(ip, cell);
           });
  deposit_level_center_mm(
      lb, dx, invVol,
      [=] AMREX_GPU_DEVICE(int, int ip, RealVect& dShift, Real& q) {
        stateOf(ip, dShift, q);
      },
      [=] AMREX_GPU_DEVICE(int) { return mmArr; });
}

// The data of a tile for the level deposit of sum_moments().
template <class P> struct MomentTile {
  const P* p;
  Array4<Real> arr;
};

// The data of a tile for the level deposits of sum_to_center().
template <class P> struct CenterTile {
  const P* p;
  Array4<Real> charge;
  Array4<RealCMM> mm;
};

} // namespace

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::sum_to_center(MultiFab& netChargeMF,
                                                       CenterMMFab& centerMM,
                                                       bool doNetChargeOnly,
                                                       int iLev) {
  timing_func("Pts::sum_to_center");

  const auto plo = Geom(iLev).ProbLoArray();
  const auto dxi = Geom(iLev).InvCellSizeArray();
  const Real invV = invVol[iLev];
  const int iqp = iqp_;
  Real dxLoc[3] = { 0, 0, 0 };
  for (int d = 0; d < nDim; ++d)
    dxLoc[d] = dx[iLev][d];

  if constexpr (useCellBinDeposits) {
    // GPU: the deposits of all the tiles at once (see LevelBins). The net
    // charge and the mass matrix have the same boxes, so they share the bins.
    std::vector<int> nps;
    std::vector<Box> cellBoxes;
    Vector<CenterTile<ParticleType> > tiles;
    for (PIter pti(*this, iLev, stream_ordered_pti()); pti.isValid(); ++pti) {
      // The cells whose 2^D-cell stencil lies in the FABs.
      Box cellBox = netChargeMF[pti].box();
      cellBox.growHi(-1);
      AMREX_ALWAYS_ASSERT(doNetChargeOnly ||
                          centerMM[pti].box() == netChargeMF[pti].box());
      nps.push_back(pti.numParticles());
      cellBoxes.push_back(cellBox);
      tiles.push_back(
          { pti.GetArrayOfStructs()().data(), netChargeMF[pti].array(),
            doNetChargeOnly ? Array4<RealCMM>() : centerMM[pti].array() });
    }
    if (tiles.empty())
      return;
    Gpu::AsyncArray<CenterTile<ParticleType> > tileStore(tiles.data(),
                                                         tiles.size());
    const CenterTile<ParticleType>* const tile = tileStore.data();

    // IDs of ghost cell particles are set to -1 inside divE_correct_position()
    // but these particles count here (see below).
    LevelBins lb;
    lb.build(nps, cellBoxes,
             [=] AMREX_GPU_DEVICE(int t, int ip, IntVect& cell) {
               ParticleInterpolator::Linear interp(tile[t].p[ip], plo, dxi);
               cell = IntVect(AMREX_D_DECL(interp.index[0], interp.index[1],
                                           interp.index[2]));
               return true;
             });
    deposit_level_to_nodes<1>(
        lb,
        [=] AMREX_GPU_DEVICE(int t, int ip, auto&& acc) {
          const ParticleType& p = tile[t].p[ip];
          ParticleInterpolator::Linear interp(p, plo, dxi);
          const Real pval = p.rdata(iqp) * invV;
          for (int kk = 0; kk <= ParticleInterpolator::Linear::nz; ++kk)
            for (int jj = 0; jj <= ParticleInterpolator::Linear::ny; ++jj)
              for (int ii = 0; ii <= ParticleInterpolator::Linear::nx; ++ii)
                acc[ii][jj][kk][0] +=
                    interp.w[ii] * interp.w[2 + jj] * interp.w[4 + kk] * pval;
        },
        [=] AMREX_GPU_DEVICE(int t, int i, int j, int k, int n, Real v) {
          tile[t].charge(i, j, k, n) += v;
        });
    if (!doNetChargeOnly)
      deposit_level_center_mm(
          lb, dxLoc, invV,
          [=] AMREX_GPU_DEVICE(int t, int ip, RealVect& dShift, Real& q) {
            const ParticleType& p = tile[t].p[ip];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            dShift =
                RealVect(AMREX_D_DECL(interp.w[1], interp.w[3], interp.w[5]));
            q = p.rdata(iqp);
          },
          [=] AMREX_GPU_DEVICE(int t) { return tile[t].mm; });
    return;
  }

  for (PIter pti(*this, iLev, stream_ordered_pti()); pti.isValid(); ++pti) {
    Array4<Real> const chargeArr = netChargeMF[pti].array();
    Array4<RealCMM> const mmArr = centerMM[pti].array();
    const ParticleType* const pData = pti.GetArrayOfStructs()().data();
    const int np = pti.numParticles();

    /*
    Q: Why do not check p.id() < 0?
    A: IDs of ghost cell particles are set to -1 inside
    divE_correct_position(), but these particles should be take into account
    here.
    */
    auto cellOf = [=] AMREX_GPU_DEVICE(int ip, IntVect& cell) {
      ParticleInterpolator::Linear interp(pData[ip], plo, dxi);
      cell = IntVect(
          AMREX_D_DECL(interp.index[0], interp.index[1], interp.index[2]));
      return true;
    };
    // The cells whose 2^D-cell stencil lies in the FABs.
    Box cellBox = netChargeMF[pti].box();
    cellBox.growHi(-1);

    // Net charge: the product of the weights in the order of
    // ParticleInterpolator::ParticleToMesh.
    auto chargeContrib = [=] AMREX_GPU_DEVICE(int ip, auto&& acc) {
      const ParticleType& p = pData[ip];
      ParticleInterpolator::Linear interp(p, plo, dxi);
      const Real pval = p.rdata(iqp) * invV;
      for (int kk = 0; kk <= ParticleInterpolator::Linear::nz; ++kk)
        for (int jj = 0; jj <= ParticleInterpolator::Linear::ny; ++jj)
          for (int ii = 0; ii <= ParticleInterpolator::Linear::nx; ++ii)
            acc[ii][jj][kk][0] +=
                interp.w[ii] * interp.w[2 + jj] * interp.w[4 + kk] * pval;
    };

    if constexpr (!useCellBinDeposits) {
      if (!doNetChargeOnly) {
        // CPU: the charge and the mass-matrix terms of a particle in one pass
        // (one read and one interpolation per particle), in particle order.
        auto addCharge = [=](int i, int j, int k, int, Real v) {
          chargeArr(i, j, k, 0) += v;
        };
        for (int ip = 0; ip < np; ++ip) {
          IntVect cell;
          cellOf(ip, cell);
          chargeContrib(ip, NodeAdder<decltype(addCharge)>{ &addCharge, cell });
          const ParticleType& p = pData[ip];
          ParticleInterpolator::Linear interp(p, plo, dxi);
          const RealVect dShift(
              AMREX_D_DECL(interp.w[1], interp.w[3], interp.w[5]));
          accumulate_mass_matrix_contribution(dxLoc, invV, cell, dShift,
                                              p.rdata(iqp), mmArr);
        }
        continue;
      }
    }

    deposit_to_nodes<1>(np, cellBox, chargeArr, 0, cellOf, chargeContrib);

    if (doNetChargeOnly)
      continue;

    {
      // GPU: per-cell sums of the mass-matrix terms (see deposit_to_nodes).
      Box mmBox = centerMM[pti].box();
      mmBox.growHi(-1);
      deposit_center_mm_by_cell(
          np, mmBox, mmArr, dxLoc, invV, cellOf,
          [=] AMREX_GPU_DEVICE(int ip, RealVect& dShift, Real& q) {
            const ParticleType& p = pData[ip];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            dShift =
                RealVect(AMREX_D_DECL(interp.w[1], interp.w[3], interp.w[5]));
            q = p.rdata(iqp);
          });
    }
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::sum_to_center_amr(
    MultiFab& netChargeMF, MultiFab& jc, MultiFab& jf, CenterMMFab& centerMM,
    bool doNetChargeOnly, int iLev) {
  timing_func("Pts::sum_to_center");

  int finer_level = iLev + 1;
  int coarser_level = iLev - 1;
  if (iLev == 0) {
    coarser_level = iLev;
  }
  if (iLev == (n_lev() - 1)) {
    finer_level = iLev;
  }

  // Copies, so that the kernels capture values, not host pointers.
  Real ploLev[3] = { 0, 0, 0 }, invDxLev[3] = { 0, 0, 0 };
  Real dxLev[3] = { 0, 0, 0 };
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    ploLev[d] = Geom(iLev).ProbLo(d);
    invDxLev[d] = Geom(iLev).InvCellSize(d);
  }
  for (int d = 0; d < nDim; ++d)
    dxLev[d] = dx[iLev][d];
  const Real invVolLev = invVol[iLev];
  const bool multiLevel = n_lev() > 1;
  const int iqp = iqp_;

  for (int nLev = finer_level; nLev >= coarser_level; nLev--) {
    Real ploN[3] = { 0, 0, 0 }, invDxN[3] = { 0, 0, 0 };
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      ploN[d] = Geom(nLev).ProbLo(d);
      invDxN[d] = Geom(nLev).InvCellSize(d);
    }

    if (iLev == nLev) {
      for (PIter pti(*this, nLev, stream_ordered_pti()); pti.isValid(); ++pti) {
        Array4<Real> const chargeArr = netChargeMF[pti].array();

        Array4<RealCMM> const mmArr = centerMM[pti].array();
        const Array4<int const> status = cell_status(nLev)[pti].const_array();
        const ParticleType* const pData = pti.GetArrayOfStructs()().data();
        const Dim3 lo = init_dim3(0);
        const Dim3 hi = init_dim3(1);
        const int np = pti.numParticles();
        if (useCellBinDeposits) {
          // GPU: per-cell sums in particle-key order (see deposit_to_nodes).
          // On this level every particle deposits its charge.
          auto cellOf = [=] AMREX_GPU_DEVICE(int ip, IntVect& cell) {
            IntVect loIdx;
            RealVect dShift;
            Real coef[2][2][2];
            find_cell_interpolation(pData[ip].pos(), ploLev, invDxLev, loIdx,
                                    dShift, coef);
            cell = loIdx;
            return true;
          };
          Box cellBox = netChargeMF[pti].box();
          cellBox.growHi(-1);
          deposit_to_nodes<1>(
              np, cellBox, chargeArr, 0, cellOf,
              [=] AMREX_GPU_DEVICE(int ip, auto&& acc) {
                IntVect loIdx;
                RealVect dShift;
                Real coef[2][2][2];
                find_cell_interpolation(pData[ip].pos(), ploLev, invDxLev,
                                        loIdx, dShift, coef);
                const Real cTmp = pData[ip].rdata(iqp) * invVolLev;
                for (int kk = lo.z; kk <= hi.z; ++kk)
                  for (int jj = lo.y; jj <= hi.y; ++jj)
                    for (int ii = lo.x; ii <= hi.x; ++ii)
                      acc[ii][jj][kk][0] += coef[ii][jj][kk] * cTmp;
              });
          if (doNetChargeOnly)
            continue;
          Box mmBox = centerMM[pti].box();
          mmBox.growHi(-1);
          deposit_center_mm_by_cell(
              np, mmBox, mmArr, dxLev, invVolLev,
              [=] AMREX_GPU_DEVICE(int ip, IntVect& cell) {
                if (multiLevel &&
                    skip_particle_for_dive_cleaning(pData[ip].pos(), ploLev,
                                                    invDxLev, iLev, status))
                  return false;
                return cellOf(ip, cell);
              },
              [=] AMREX_GPU_DEVICE(int ip, RealVect& dShift, Real& q) {
                IntVect loIdx;
                Real coef[2][2][2];
                find_cell_interpolation(pData[ip].pos(), ploLev, invDxLev,
                                        loIdx, dShift, coef);
                q = pData[ip].rdata(iqp);
              });
          continue;
        }
        ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) {
          const ParticleType& p = pData[ip];
          const Real qp = p.rdata(iqp);
          IntVect loIdx;
          RealVect dShift;
          IntVect realIdx;
          RealVect tmprv;
          find_cell_index_exp(p.pos(), ploN, invDxN, realIdx, tmprv);
          if (nLev == iLev || (bit::is_refined_neighbour(status(realIdx)) ||
                               bit::is_lev_edge(status(realIdx)))) {
            Real coef[2][2][2];
            find_cell_interpolation(p.pos(), ploLev, invDxLev, loIdx, dShift,
                                    coef);
            const Real cTmp = qp * invVolLev;
            deposit_charge(chargeArr, loIdx, coef, cTmp, lo, hi);
          }

          bool skipParticle = false;
          if (multiLevel) {
            skipParticle = skip_particle_for_dive_cleaning(
                p.pos(), ploLev, invDxLev, iLev, status);
          }
          if (!doNetChargeOnly && !skipParticle) {
            accumulate_mass_matrix_contribution(dxLev, invVolLev, loIdx, dShift,
                                                qp, mmArr);
          }
        });
      }
    }

    // A particle of level nLev deposits onto the charge of level iLev in jf
    // (finer level) or jc (coarser level).
    if (nLev != iLev) {
      MultiFab& depositMF = (nLev > iLev) ? jf : jc;
      for (PIter pti(*this, nLev, stream_ordered_pti()); pti.isValid(); ++pti) {
        Array4<Real> const chargeArr = depositMF[pti].array();
        const Array4<int const> status = cell_status(nLev)[pti].const_array();
        const ParticleType* const pData = pti.GetArrayOfStructs()().data();
        const Dim3 lo = init_dim3(0);
        const Dim3 hi = init_dim3(1);
        const int np = pti.numParticles();
        // A particle deposits if its cell (on its own level) is next to the
        // refinement or on the level edge.
        auto deposits = [=] AMREX_GPU_DEVICE(int ip) {
          IntVect realIdx;
          RealVect tmprv;
          find_cell_index_exp(pData[ip].pos(), ploN, invDxN, realIdx, tmprv);
          return bit::is_refined_neighbour(status(realIdx)) ||
                 bit::is_lev_edge(status(realIdx));
        };
        if (useCellBinDeposits) {
          // GPU: per-cell sums in particle-key order (see deposit_to_nodes).
          Box cellBox = depositMF[pti].box();
          cellBox.growHi(-1);
          deposit_to_nodes<1>(
              np, cellBox, chargeArr, 0,
              [=] AMREX_GPU_DEVICE(int ip, IntVect& cell) {
                if (!deposits(ip))
                  return false;
                IntVect loIdx;
                RealVect dShift;
                Real coef[2][2][2];
                find_cell_interpolation(pData[ip].pos(), ploLev, invDxLev,
                                        loIdx, dShift, coef);
                cell = loIdx;
                return true;
              },
              [=] AMREX_GPU_DEVICE(int ip, auto&& acc) {
                IntVect loIdx;
                RealVect dShift;
                Real coef[2][2][2];
                find_cell_interpolation(pData[ip].pos(), ploLev, invDxLev,
                                        loIdx, dShift, coef);
                const Real cTmp = pData[ip].rdata(iqp) * invVolLev;
                for (int kk = lo.z; kk <= hi.z; ++kk)
                  for (int jj = lo.y; jj <= hi.y; ++jj)
                    for (int ii = lo.x; ii <= hi.x; ++ii)
                      acc[ii][jj][kk][0] += coef[ii][jj][kk] * cTmp;
              });
          continue;
        }
        ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) {
          const ParticleType& p = pData[ip];
          const Real qp = p.rdata(iqp);
          IntVect loIdx;
          RealVect dShift;
          if (deposits(ip)) {
            Real coef[2][2][2];
            find_cell_interpolation(p.pos(), ploLev, invDxLev, loIdx, dShift,
                                    coef);
            const Real cTmp = qp * invVolLev;
            deposit_charge(chargeArr, loIdx, coef, cTmp, lo, hi);
          }
        });
      }
    }
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
std::array<Real, 5> Particles<NStructReal, NStructInt>::total_moments(
    bool localOnly) {
  timing_func("Pts::total_moments");

  const int iLev = 0;
  const int iup = iup_;
  const int ivp = ivp_;
  const int iwp = iwp_;
  const int iqp = iqp_;

  // Charge, momentum and energy, one partial sum per tile added in tile
  // order (the order of amrex::ParticleReduce in a CPU build), every partial
  // a deterministic_sum.
  DeferredSums partials;
  for (PIter pti(*this, iLev, stream_ordered_pti()); pti.isValid(); ++pti) {
    const ParticleType* const pData = pti.GetArrayOfStructs()().data();
    for (int c = 0; c < 5; ++c) {
      partials.add(pti.numParticles(), [=] AMREX_GPU_DEVICE(int ip) -> Real {
        const ParticleType& p = pData[ip];
        if (p.id() < 0)
          return 0.0;
        const Real up = p.rdata(iup);
        const Real vp = p.rdata(ivp);
        const Real wp = p.rdata(iwp);
        const Real qp = p.rdata(iqp);
        switch (c) {
          case 0:
            return qp;
          case 1:
            return qp * up;
          case 2:
            return qp * vp;
          case 3:
            return qp * wp;
          default:
            return 0.5 * qp * (up * up + vp * vp + wp * wp);
        }
      });
    }
  }
  Real r[5] = { 0, 0, 0, 0, 0 };
  const std::vector<Real> partial = partials.values();
  for (std::size_t t = 0; t < partial.size(); t += 5)
    for (int c = 0; c < 5; ++c)
      r[c] += partial[t + c];

  const Real factor = qomSign * get_mass();
  std::array<Real, 5> sum = { r[0] * factor, r[1] * factor, r[2] * factor,
                              r[3] * factor, r[4] * factor };

  if (!localOnly) {
    ParallelDescriptor::ReduceRealSum(sum.data(), sum.size(),
                                      ParallelDescriptor::IOProcessorNumber());
  }

  return sum;
}

//==========================================================

//==========================================================
// Mirror boundary condition for a node-centered deposit at non-periodic
// uncoupled walls (see the comment in sum_moments): double the edge node (or
// clear it for an odd component at a reflecting wall) and clear the first
// ghost-node layer. At a reflecting wall the ghost layer is first folded onto
// the first interior node with the same parity: it holds the part of the
// charge-conserving current of a particle that crossed the wall before it was
// mirrored, and the folded current is that of the reflected path (image
// method). isVector selects the parity rule: components (x, y, z) of a vector
// (odd: the normal one) instead of the moment layout (odd: normal momentum
// and normal-shear stresses).
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::fold_wall_ghost_nodes(MultiFab& mf,
                                                               int iLev,
                                                               bool isVector) {
  if (Geom(iLev).isAllPeriodic() || mf.nGrow() == 0)
    return;

  const Box& dom = Geom(iLev).Domain();
  const int nCompMF = mf.nComp();
  for (MFIter mfi(mf, stream_ordered_mfi()); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.validbox();
    Array4<Real> const& arr = mf[mfi].array();

    for (int iDim = 0; iDim < nDim; ++iDim) {
      if (Geom(iLev).isPeriodic(iDim))
        continue;

      for (int side = 0; side < 2; ++side) {
        const bool isLo = (side == 0);
        const int faceBc = isLo ? bc.lo[iDim] : bc.hi[iDim];
        // A particle face can no longer be `conducting`: that is a
        // field-only type, and on the particle side it is mapped to
        // `reflect` at parse time. The fold set covers all uncoupled,
        // non-periodic faces (reflect, inflow, outflow, vacuum, absorb)
        // where no exterior ghost particles exist to deposit into the
        // boundary node.
        const bool doFold =
            (faceBc == ParticleBC::reflect || faceBc == ParticleBC::inflow ||
             faceBc == ParticleBC::outflow || faceBc == ParticleBC::vacuum ||
             faceBc == ParticleBC::absorb);
        if (!doFold)
          continue;

        // NOTE: momentsMF is NODE-centered (nGrids = convert(cGrids,
        // nodeVector)), whose domain box extends one node beyond the
        // cell-centered Geometry box on the hi side.  The edge NODE is
        // dom.smallEnd (lo) or dom.bigEnd+1 (hi); using dom.bigEnd would
        // hit the second-last node and skip the real edge tile.
        const int domEdge = isLo ? dom.smallEnd(iDim) : dom.bigEnd(iDim) + 1;
        const int tileEdge = isLo ? bx.smallEnd(iDim) : bx.bigEnd(iDim);
        if (tileEdge != domEdge)
          continue;

        // Strip of the first ghost-node layer just outside the domain,
        // spanning the full transverse extent of this FAB (valid +
        // transverse ghosts).  Mirror the edge node into it, fold it into
        // the edge node (doubling), then clear the ghost.
        const Box& fbx = mfi.fabbox();
        IntVect gs = fbx.smallEnd();
        IntVect ge = fbx.bigEnd();
        gs[iDim] = domEdge + (isLo ? -1 : 1);
        ge[iDim] = gs[iDim];
        const Box strip(gs, ge);

        const bool isReflect = (faceBc == ParticleBC::reflect);

        ParallelFor(
            strip, nCompMF, [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) {
              // Odd parity components relative to iDim under
              // specular reflection: normal momentum and
              // off-diagonal normal-shear stresses.
              bool isOdd = false;
              if (isVector) {
                isOdd = (c == iDim);
              } else if (iDim == 0) {
                isOdd = (c == iMx_ || c == iPxy_ || c == iPxz_);
              } else if (iDim == 1) {
                isOdd = (c == iMy_ || c == iPxy_ || c == iPyz_);
              } else if (iDim == 2) {
                isOdd = (c == iMz_ || c == iPxz_ || c == iPyz_);
              }

#if (AMREX_SPACEDIM == 2)
              int ei = i, ej = j;
              if (iDim == 0)
                ei = domEdge;
              if (iDim == 1)
                ej = domEdge;
              const Real val = arr(ei, ej, 0, c);
              if (isReflect && isOdd) {
                arr(ei, ej, 0, c) = 0.0;
              } else {
                arr(ei, ej, 0, c) += val; // fold mirror back -> 2x
              }
              if (isReflect) {
                const int ii = iDim == 0 ? 2 * domEdge - i : i;
                const int jj = iDim == 1 ? 2 * domEdge - j : j;
                arr(ii, jj, 0, c) += isOdd ? -arr(i, j, 0, c) : arr(i, j, 0, c);
              }
              arr(i, j, 0, c) = 0.0; // zero the ghost layer
#elif (AMREX_SPACEDIM == 3)
                int ei = i, ej = j, ek = k;
                if (iDim == 0) ei = domEdge;
                if (iDim == 1) ej = domEdge;
                if (iDim == 2) ek = domEdge;
                const Real val = arr(ei, ej, ek, c);
                if (isReflect && isOdd) {
                  arr(ei, ej, ek, c) = 0.0;
                } else {
                  arr(ei, ej, ek, c) += val; // fold mirror back -> 2x
                }
                if (isReflect) {
                  const int ii = iDim == 0 ? 2 * domEdge - i : i;
                  const int jj = iDim == 1 ? 2 * domEdge - j : j;
                  const int kk = iDim == 2 ? 2 * domEdge - k : k;
                  arr(ii, jj, kk, c) +=
                      isOdd ? -arr(i, j, k, c) : arr(i, j, k, c);
                }
                arr(i, j, k, c) = 0.0; // zero the ghost layer
#endif
            });
      }
    }
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
Real Particles<NStructReal, NStructInt>::sum_moments(
    Vector<MultiFab>& momentsMF, Vector<MultiFab>& nodeBMF, Real dt) {
  timing_func("Pts::sum_moments");

  // The kinetic energy is summed on the fly (see KernelSum).
  KernelSum energySum;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    timing_func("Pts::sum_moments_node_deposit");
    momentsMF[iLev].setVal(0.0);
    const auto ploLev = Geom(iLev).ProbLoArray();
    const auto invDxLev = Geom(iLev).InvCellSizeArray();
    // The cell of a particle is the cell whose lower node is the lower node
    // of its interpolation stencil.
    auto cellOfParticle = [=] AMREX_GPU_DEVICE(const ParticleType& p,
                                               IntVect& cell) {
      if (p.id() < 0)
        return false;
      RealVect dShift;
      Real coef[2][2][2];
      find_node_interpolation(p.pos(), ploLev.data(), invDxLev.data(), cell,
                              dShift, coef);
      return true;
    };
    auto momentsOfParticle = [=] AMREX_GPU_DEVICE(const ParticleType& p,
                                                  auto&& acc) {
      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      IntVect loIdx;
      RealVect dShift;
      Real coef[2][2][2];
      find_node_interpolation(p.pos(), ploLev.data(), invDxLev.data(), loIdx,
                              dShift, coef);
      //-----calculate interpolate coef end-------------

      //-------nodePlasma begin---------
      Real pMoments[nMoments];

      pMoments[iNum_] = 1;
      pMoments[iRho_] = qp;

      {
        const Real mx = qp * up;
        const Real my = qp * vp;
        const Real mz = qp * wp;
        pMoments[iMx_] = mx;
        pMoments[iMy_] = my;
        pMoments[iMz_] = mz;

        pMoments[iPxx_] = mx * up;
        pMoments[iPyy_] = my * vp;
        pMoments[iPzz_] = mz * wp;

        pMoments[iPxy_] = mx * vp;
        pMoments[iPxz_] = mx * wp;
        pMoments[iPyz_] = my * wp;
      }

      for (int iVar = 0; iVar < nMoments; iVar++)
        for (int kk = 0; kk < 2; ++kk)
          for (int jj = 0; jj < 2; ++jj)
            for (int ii = 0; ii < 2; ++ii)
              acc[ii][jj][kk][iVar] += coef[ii][jj][kk] * pMoments[iVar];
      //-------nodePlasma end---------
    };

    // A GPU build deposits all the tiles at once (see LevelBins).
    std::vector<int> nps;
    std::vector<Box> cellBoxes;
    Vector<MomentTile<ParticleType> > tiles;
    for (PIter pti(*this, iLev, stream_ordered_pti()); pti.isValid(); ++pti) {
      Array4<Real> const momentsArr = momentsMF[iLev][pti].array();

      const ParticleType* const pData = pti.GetArrayOfStructs()().data();

      const int np = pti.numParticles();

      const Box cellBox = enclosedCells(momentsMF[iLev][pti].box());
      if constexpr (useCellBinDeposits) {
        nps.push_back(np);
        cellBoxes.push_back(cellBox);
        tiles.push_back({ pData, momentsArr });
      } else {
        deposit_to_nodes<nMoments>(
            np, cellBox, momentsArr, 0,
            [=] AMREX_GPU_DEVICE(int ip, IntVect& cell) {
              return cellOfParticle(pData[ip], cell);
            },
            [=] AMREX_GPU_DEVICE(int ip, auto&& acc) {
              momentsOfParticle(pData[ip], acc);
            });
      }

      energySum.add(np, [=] AMREX_GPU_DEVICE(int ip) -> GpuTuple<Real> {
        const ParticleType& p = pData[ip];
        if (p.id() < 0)
          return { 0.0 };
        const Real up = p.rdata(iup_);
        const Real vp = p.rdata(ivp_);
        const Real wp = p.rdata(iwp_);
        return { p.rdata(iqp_) * (up * up + vp * vp + wp * wp) };
      });
    }

    if (!tiles.empty()) {
      Gpu::AsyncArray<MomentTile<ParticleType> > tileStore(tiles.data(),
                                                           tiles.size());
      const MomentTile<ParticleType>* const tile = tileStore.data();
      LevelBins lb;
      lb.build(nps, cellBoxes,
               [=] AMREX_GPU_DEVICE(int t, int ip, IntVect& cell) {
                 return cellOfParticle(tile[t].p[ip], cell);
               });
      deposit_level_to_nodes<nMoments>(
          lb,
          [=] AMREX_GPU_DEVICE(int t, int ip, auto&& acc) {
            momentsOfParticle(tile[t].p[ip], acc);
          },
          [=] AMREX_GPU_DEVICE(int t, int i, int j, int k, int n, Real v) {
            tile[t].arr(i, j, k, n) += v;
          });
    }

    // Exclude the number density.
    momentsMF[iLev].mult(invVol[iLev], 0, nMoments - 1,
                         momentsMF[iLev].nGrow());

    //----- Mirror boundary condition for the node-centered deposit ------------
    // A node-centered CIC deposit at a non-periodic wall leaves the boundary
    // node with only ~half the charge of an interior node: with the node at the
    // domain face, a particle in the edge cell deposits its weight between that
    // edge node and the next interior node, and NO particle lies on the
    // exterior side to contribute the complementary weight.  This is why the
    // reflecting-wall node of the 1D shock test showed rhoS0 ~ 6.25 instead of
    // ~12.5 before this fold was introduced.
    //
    // The physically correct boundary condition for a reflecting wall (and an
    // open-inflow face, whose injected particles live just inside the edge
    // cell) is MIRROR SYMMETRY across the boundary node: the ghost node
    // must carry the same moment value as the edge node (Neumann mirror).  We
    // therefore copy the edge node into the first ghost node and fold the ghost
    // back into the edge node (equivalent to doubling the edge node), then zero
    // the ghost layer so SumBoundary does not re-use it.  We apply this ONLY at
    // reflect / inflow faces (`conducting` is a field-only type; on the
    // particle side it is mapped to `reflect` at parse time).  At coupled
    // (MHD-AEPIC), outflow, vacuum, absorb, and periodic faces the
    // ghost layer may hold real injected/leaving particles, so we leave it
    // untouched to avoid double counting.  Periodic directions are skipped
    // (SumBoundary wraps them).
    // NOTE: the fold runs BEFORE SumBoundary on purpose so a transverse ghost
    // target that lies in a tile-neighbour's valid region is combined exactly
    // once.
    fold_wall_ghost_nodes(momentsMF[iLev], iLev, false);

    sum_boundary(momentsMF[iLev], Geom(iLev).periodicity());
  }

  for (int iLev = n_lev() - 2; iLev >= 0; iLev--) {
    timing_func("Pts::sum_moments_coarse_fine_interface");
    sum_two_lev_interface_node(
        momentsMF[iLev], momentsMF[iLev + 1], 0, momentsMF[iLev].nComp(),
        get_ref_ratio(iLev), Geom(iLev), Geom(iLev + 1), node_status(iLev + 1));
  }

  // Correct domain edge nodes
  for (int iLev = 0; iLev < n_lev() - 1; iLev++) {
    timing_func("Pts::sum_moments_domain_edge_correction");
    interp_from_coarse_to_fine_for_domain_edge(
        momentsMF[iLev], momentsMF[iLev + 1], 0, momentsMF[iLev].nComp(),
        get_ref_ratio(iLev), Geom(iLev), Geom(iLev + 1), node_status(iLev + 1));
  }

  Real energy = energySum.value();
  energy *= 0.5 * qomSign * get_mass();

  return energy;
}

//==========================================================

template <int NStructReal, int NStructInt>
Real Particles<NStructReal, NStructInt>::calc_max_thermal_velocity(
    MultiFab& momentsMF) {

  constexpr Real c1over3 = 1. / 3;
  auto const& ma = momentsMF.const_arrays();

  Real uthMax = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, momentsMF,
                          IntVect(0), // zero ghost cells
                          [=] AMREX_GPU_DEVICE(int nb, int i, int j, int k)
                              noexcept -> GpuTuple<Real> {
                                Array4<Real const> const& arr = ma[nb];
                                Real rho = arr(i, j, k, iRho_);
                                if (rho == 0)
                                  return 0.0;
                                Real p =
                                    (arr(i, j, k, iPxx_) + arr(i, j, k, iPyy_) +
                                     arr(i, j, k, iPzz_)) *
                                    c1over3;
                                Real uth = sqrt(p / rho);
                                return uth;
                              });

  return uthMax;
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::convert_to_fluid_moments(
    Vector<MultiFab>& momentsMF) {

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    MultiFab tmpMF(momentsMF[iLev], make_alias, iRho_, iPyz_ - iRho_ + 1);
    tmpMF.mult(qomSign * get_mass(), tmpMF.nGrow());

    for (MFIter mfi(momentsMF[iLev], stream_ordered_mfi()); mfi.isValid();
         ++mfi) {
      const Array4<Real> arr = momentsMF[iLev][mfi].array();

      ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const Real rho = arr(i, j, k, iRho_);
        if (rho > 0) {
          const Real ux = arr(i, j, k, iUx_) / rho;
          const Real uy = arr(i, j, k, iUy_) / rho;
          const Real uz = arr(i, j, k, iUz_) / rho;
          arr(i, j, k, iPxx_) = arr(i, j, k, iPxx_) - rho * ux * ux;
          arr(i, j, k, iPyy_) = arr(i, j, k, iPyy_) - rho * uy * uy;
          arr(i, j, k, iPzz_) = arr(i, j, k, iPzz_) - rho * uz * uz;

          arr(i, j, k, iPxy_) = arr(i, j, k, iPxy_) - rho * ux * uy;
          arr(i, j, k, iPxz_) = arr(i, j, k, iPxz_) - rho * ux * uz;
          arr(i, j, k, iPyz_) = arr(i, j, k, iPyz_) - rho * uy * uz;
        }
      });
    }
  }
}

//==========================================================

template void PicParticles::sum_to_center(MultiFab&, CenterMMFab&, bool, int);
template void PTParticles::sum_to_center(MultiFab&, CenterMMFab&, bool, int);
template void PicParticles::sum_to_center_amr(MultiFab&, MultiFab&, MultiFab&,
                                              CenterMMFab&, bool, int);
template void PTParticles::sum_to_center_amr(MultiFab&, MultiFab&, MultiFab&,
                                             CenterMMFab&, bool, int);
template std::array<Real, 5> PicParticles::total_moments(bool);
template std::array<Real, 5> PTParticles::total_moments(bool);
template Real PicParticles::sum_moments(Vector<MultiFab>&, Vector<MultiFab>&,
                                        Real);
template void PicParticles::fold_wall_ghost_nodes(MultiFab&, int, bool);
template void PTParticles::fold_wall_ghost_nodes(MultiFab&, int, bool);
template Real PTParticles::sum_moments(Vector<MultiFab>&, Vector<MultiFab>&,
                                       Real);
template Real PicParticles::calc_max_thermal_velocity(MultiFab&);
template Real PTParticles::calc_max_thermal_velocity(MultiFab&);
template void PicParticles::convert_to_fluid_moments(Vector<MultiFab>&);
template void PTParticles::convert_to_fluid_moments(Vector<MultiFab>&);

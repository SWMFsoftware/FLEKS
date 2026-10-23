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

// Deterministic deposit of jHat and of the kp > 0 half of the node mass
// matrix (see calc_mass_matrix) for one tile: the particles of each cell
// are summed in a fixed order (CellBins), then every node gathers from the
// cells around it. No floating-point atomics, so the result is the same from
// run to run.
//
// With the trilinear weight w0 = 1 - s, w1 = s of the lower and upper node in
// each direction, the mass-matrix entry of a node pair in a cell is
//   sum_p alpha_p prod_d w_a(s_d) w_b(s_d),
// and w_a w_b is a quadratic in s_d: (1-s)^2, s(1-s) or s^2. So a cell only
// needs the 27 moments sum_p alpha_p sx^l sy^m sz^n (l, m, n = 0, 1, 2) of
// each of the 9 alpha components, and the entries are fixed combinations of
// them: 243 sums per cell instead of 432 scattered additions per particle.
//
// The entries of the pairs (node, node + o) with o_z + 1 >= kpMin are filled:
// kpMin = 1 for a single level (the kp = 0 half is mirrored afterwards),
// kpMin = 0 for the AMR mass matrix, which stores all of them.
template <class CellOf, class IdOf, class State>
void deposit_mass_matrix_by_cell(int np, const Box& cellBox,
                                 const Array4<Real>& jArr,
                                 const Array4<RealMM>& mmArr, CellOf cellOf,
                                 IdOf idOf, State state, int kpMin = 1) {
  if (np == 0)
    return;
  const Dim3 lo = lbound(cellBox);
  const Dim3 len = length(cellBox);
  const int nBins = static_cast<int>(cellBox.numPts());

  CellBins bins;
  bins.build(
      np, nBins,
      [=] AMREX_GPU_DEVICE(int ip) {
        IntVect c;
        if (!cellOf(ip, c))
          return nBins;
        const int i = c[0] - lo.x, j = c[1] - lo.y, k = c[2] - lo.z;
        if (i < 0 || i >= len.x || j < 0 || j >= len.y || k < 0 || k >= len.z)
          return nBins;
        return i + len.x * (j + len.y * k);
      },
      idOf);

  // 1. Per cell: the jHat contribution of its 8 nodes and the 27 x 9 moments.
  constexpr int nJ = 8 * 3, nS = 27 * 9, nAcc = nJ + nS;
  Gpu::DeviceVector<Real> scratch(static_cast<std::size_t>(nBins) * nAcc);
  Real* const sc = scratch.data();
  const int* const perm = bins.perm();
  const int* const off = bins.offsets();
  ParallelFor(nBins, [=] AMREX_GPU_DEVICE(int b) {
    Real acc[nAcc];
    for (int n = 0; n < nAcc; ++n)
      acc[n] = 0;
    for (int i = off[b]; i < off[b + 1]; ++i) {
      IntVect loIdx;
      RealVect dShift;
      Real coef[2][2][2], alpha[9], currents[3];
      state(perm[i], loIdx, dShift, coef, alpha, currents);
      for (int kk = 0; kk < 2; ++kk)
        for (int jj = 0; jj < 2; ++jj)
          for (int ii = 0; ii < 2; ++ii)
            for (int d = 0; d < 3; ++d)
              acc[((kk * 2 + jj) * 2 + ii) * 3 + d] +=
                  coef[ii][jj][kk] * currents[d];
      const Real sx[3] = { 1, dShift[0], dShift[0] * dShift[0] };
      const Real sy[3] = { 1, dShift[1], dShift[1] * dShift[1] };
      const Real sz[3] = { 1, dShift[2], dShift[2] * dShift[2] };
      for (int n = 0; n < 3; ++n)
        for (int m = 0; m < 3; ++m)
          for (int l = 0; l < 3; ++l) {
            const Real w = sx[l] * sy[m] * sz[n];
            Real* const a = acc + nJ + ((n * 3 + m) * 3 + l) * 9;
            for (int idx = 0; idx < 9; ++idx)
              a[idx] += alpha[idx] * w;
          }
    }
    Real* const out = sc + static_cast<std::size_t>(b) * nAcc;
    for (int n = 0; n < nAcc; ++n)
      out[n] = acc[n];
  });

  // 2. Per node: gather jHat and the mass-matrix entries from the cells.
  ParallelFor(surroundingNodes(cellBox), [=] AMREX_GPU_DEVICE(int i, int j,
                                                              int k) {
    const int node[3] = { i, j, k };
    const int loC[3] = { lo.x, lo.y, lo.z };
    const int lenC[3] = { len.x, len.y, len.z };
    // Polynomial of w_a w_b in s (coefficients of 1, s, s^2).
    auto poly = [](int a, int b, int m) -> Real {
      constexpr Real p00[3] = { 1, -2, 1 }, p01[3] = { 0, 1, -1 },
                     p11[3] = { 0, 0, 1 };
      return (a == 0 && b == 0) ? p00[m] : (a != b ? p01[m] : p11[m]);
    };

    // jHat: the cells that have this node as a corner.
    Real jSum[3] = { 0, 0, 0 };
    for (int kk = 0; kk < 2; ++kk)
      for (int jj = 0; jj < 2; ++jj)
        for (int ii = 0; ii < 2; ++ii) {
          const int corner[3] = { ii, jj, kk };
          int c[3];
          bool inside = true;
          for (int d = 0; d < 3; ++d) {
            c[d] = node[d] - corner[d] - loC[d];
            inside = inside && c[d] >= 0 && c[d] < lenC[d];
          }
          if (!inside)
            continue;
          const Real* const in = sc +
                                 static_cast<std::size_t>(
                                     c[0] + lenC[0] * (c[1] + lenC[1] * c[2])) *
                                     nAcc +
                                 ((kk * 2 + jj) * 2 + ii) * 3;
          for (int d = 0; d < 3; ++d)
            jSum[d] += in[d];
        }
    for (int d = 0; d < 3; ++d)
      jArr(i, j, k, d) += jSum[d];

    // Mass matrix: the pairs (node, node + o), kp = o_z + 1 >= kpMin.
    RealMM& data = mmArr(i, j, k);
    for (int kp = kpMin; kp < 3; ++kp)
      for (int jp = 0; jp < 3; ++jp)
        for (int ip = 0; ip < 3; ++ip) {
          const int o[3] = { ip - 1, jp - 1, kp - 1 };
          Real value[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
          // The cells that contain both nodes.
          for (int cz = 0; cz < 2; ++cz)
            for (int cy = 0; cy < 2; ++cy)
              for (int cx = 0; cx < 2; ++cx) {
                const int a[3] = { cx, cy, cz }; // corner of node in the cell
                int c[3], bb[3];
                bool ok = true;
                for (int d = 0; d < 3; ++d) {
                  c[d] = node[d] - a[d] - loC[d];
                  bb[d] = a[d] + o[d]; // corner of node + o
                  ok = ok && c[d] >= 0 && c[d] < lenC[d] && bb[d] >= 0 &&
                       bb[d] <= 1;
                }
                if (!ok)
                  continue;
                const Real* const in =
                    sc +
                    static_cast<std::size_t>(
                        c[0] + lenC[0] * (c[1] + lenC[1] * c[2])) *
                        nAcc +
                    nJ;
                for (int n = 0; n < 3; ++n) {
                  const Real pz = poly(a[2], bb[2], n);
                  if (pz == 0)
                    continue;
                  for (int m = 0; m < 3; ++m) {
                    const Real py = poly(a[1], bb[1], m);
                    if (py == 0)
                      continue;
                    for (int l = 0; l < 3; ++l) {
                      const Real px = poly(a[0], bb[0], l);
                      if (px == 0)
                        continue;
                      const Real w = px * py * pz;
                      const Real* const s = in + ((n * 3 + m) * 3 + l) * 9;
                      for (int idx = 0; idx < 9; ++idx)
                        value[idx] += w * s[idx];
                    }
                  }
                }
              }
          const int idx0 = kp * 81 + jp * 27 + ip * 9;
          for (int idx = 0; idx < 9; ++idx)
            data[idx0 + idx] += value[idx];
        }
  });
  // The kernels read the bins and the scratch, which are freed on return.
  Gpu::streamSynchronize();
}

} // namespace

// Whether the cell of particle p is next to a refined cell.
template <class P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
is_next_to_refined(const P& p, const Array4<int const>& status,
                   const Real* plo, const Real* invDx) {
  IntVect cell;
  RealVect dShift;
  find_cell_index_exp(p.pos(), plo, invDx, cell, dShift);
  return bit::is_refined_neighbour(status(cell));
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::calc_mass_matrix(
    NodeMMFab& nodeMM, MultiFab& jHat, MultiFab& nodeBMF,
    const MultiFab* nodeB0MF, MultiFab& u0MF, Real dt, int iLev,
    bool solveInCoMov) {
  timing_func("Pts::calc_mass_matrix");

  const Real qdto2mc = charge / mass * 0.5 * dt;
  const bool useB0 = (nodeB0MF != nullptr);
  const auto ploLev = Geom(iLev).ProbLoArray();
  const auto invDxLev = Geom(iLev).InvCellSizeArray();
  const Real invVolLev = invVol[iLev];

  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real const> const nodeBArr = nodeBMF[pti].const_array();
    const Array4<Real const> nodeB0Arr =
        useB0 ? (*nodeB0MF)[pti].const_array() : Array4<Real const>();
    Array4<Real> const jArr = jHat[pti].array();
    Array4<RealMM> const mmArr = nodeMM[pti].array();

    Array4<Real const> const u0Arr = u0MF[pti].const_array();

    const ParticleType* const pData = pti.GetArrayOfStructs()().data();

    const Dim3 lo = init_dim3(0);
    const Dim3 hi = init_dim3(1);

    const int np = pti.numParticles();

    // The state of particle ip that the mass matrix and jHat need: its
    // interpolation stencil and weights, alpha and its current.
    auto particle_state =
        [=] AMREX_GPU_DEVICE(int ip, IntVect& loIdx, RealVect& dShift,
                             Real(&coef)[2][2][2], Real(&alpha)[9],
                             Real(&currents)[3]) -> bool {
      const ParticleType& p = pData[ip];
      if (p.id() < 0)
        return false;

      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);

      //-----calculate interpolation coef begin-------------
      find_node_interpolation(p.pos(), ploLev.data(), invDxLev.data(), loIdx,
                              dShift, coef);
      //-----calculate interpolate coef end-------------

      //----- Mass matrix calculation begin--------------
      Real u0[3] = { 0, 0, 0 };
      Real bp[3] = { 0, 0, 0 };

      if (solveInCoMov && useB0) {
        // The gyro frequency that enters the mass matrix is the one of the
        // total field B1 + B0.
        const Array4<Real const> fields[3] = { nodeBArr, u0Arr, nodeB0Arr };
        Real b0[3] = { 0, 0, 0 };
        Real* values[3] = { bp, u0, b0 };
        interpolate_vector_fields(fields, loIdx, coef, lo, hi, values);
        bp[ix_] += b0[ix_];
        bp[iy_] += b0[iy_];
        bp[iz_] += b0[iz_];
      } else if (solveInCoMov) {
        const Array4<Real const> fields[2] = { nodeBArr, u0Arr };
        Real* values[2] = { bp, u0 };
        interpolate_vector_fields(fields, loIdx, coef, lo, hi, values);
      } else if (useB0) {
        const Array4<Real const> fields[2] = { nodeBArr, nodeB0Arr };
        Real b0[3] = { 0, 0, 0 };
        Real* values[2] = { bp, b0 };
        interpolate_vector_fields(fields, loIdx, coef, lo, hi, values);
        bp[ix_] += b0[ix_];
        bp[iy_] += b0[iy_];
        bp[iz_] += b0[iz_];
      } else {
        interpolate_vector_field(nodeBArr, loIdx, coef, lo, hi, bp);
      }

      const Real omx = qdto2mc * bp[ix_];
      const Real omy = qdto2mc * bp[iy_];
      const Real omz = qdto2mc * bp[iz_];

      // end interpolation
      const Real omx2 = omx * omx;
      const Real omy2 = omy * omy;
      const Real omz2 = omz * omz;
      const Real omxomy = omx * omy;
      const Real omxomz = omx * omz;
      const Real omyomz = omy * omz;
      const Real omsq = omx2 + omy2 + omz2;
      const Real denom = 1.0 / (1.0 + omsq);

      const Real c0 = denom * invVolLev * qp * qdto2mc;

      alpha[0] = (1 + omx2) * c0;
      alpha[1] = (omz + omxomy) * c0;
      alpha[2] = (-omy + omxomz) * c0;
      alpha[3] = (-omz + omxomy) * c0;
      alpha[4] = (1 + omy2) * c0;
      alpha[5] = (omx + omyomz) * c0;
      alpha[6] = (omy + omxomz) * c0;
      alpha[7] = (-omx + omyomz) * c0;
      alpha[8] = (1 + omz2) * c0;

      {
        // jHat

        const Real up1 = up - u0[0];
        const Real vp1 = vp - u0[1];
        const Real wp1 = wp - u0[2];

        const Real udotOm1 = up1 * omx + vp1 * omy + wp1 * omz;

        {
          const Real coef1 = denom * qp;
          currents[ix_] =
              (up1 + (vp1 * omz - wp1 * omy + udotOm1 * omx)) * coef1;
          currents[iy_] =
              (vp1 + (wp1 * omx - up1 * omz + udotOm1 * omy)) * coef1;
          currents[iz_] =
              (wp1 + (up1 * omy - vp1 * omx + udotOm1 * omz)) * coef1;
        }
      }
      return true;
    };

    if (useCellBinDeposits && AMREX_SPACEDIM == 3) {
      // GPU: deterministic per-cell sums (see deposit_mass_matrix_by_cell).
      auto cellOf = [=] AMREX_GPU_DEVICE(int ip, IntVect& cell) {
        const ParticleType& p = pData[ip];
        if (p.id() < 0)
          return false;
        RealVect dShift;
        Real coef[2][2][2];
        find_node_interpolation(p.pos(), ploLev.data(), invDxLev.data(), cell,
                                dShift, coef);
        return true;
      };
      deposit_mass_matrix_by_cell(
          np, enclosedCells(nodeMM[pti].box()), jArr, mmArr, cellOf,
          [=] AMREX_GPU_DEVICE(int ip) { return particle_key(pData[ip]); },
          particle_state);
      continue;
    }

    ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) {
      IntVect loIdx;
      RealVect dShift;
      Real coef[2][2][2];
      Real alpha[9];
      Real currents[3];
      if (!particle_state(ip, loIdx, dShift, coef, alpha, currents))
        return;

      deposit_vector_field(jArr, loIdx, coef, lo, hi, currents, 3);

      const int iMin = loIdx[ix_];
      const int jMin = loIdx[iy_];
      const int kMin = nDim > 2 ? loIdx[iz_] : 0;
      const int iMax = iMin + 1;
      const int jMax = jMin + 1;
      const int kMax = nDim > 2 ? kMin + 1 : 0;

      for (int k1 = kMin; k1 <= kMax; k1++)
        for (int j1 = jMin; j1 <= jMax; j1++)
          for (int i1 = iMin; i1 <= iMax; i1++) {
            const Real wg = coef[i1 - iMin][j1 - jMin][k1 - kMin];
            auto& data0 = mmArr(i1, j1, k1);
            for (int k2 = kMin; k2 <= kMax; k2++) {
              const int kp = k2 - k1 + 1;
              if (kp > 0) {
                for (int j2 = jMin; j2 <= jMax; j2++) {
                  const int jp = j2 - j1 + 1;
                  for (int i2 = iMin; i2 <= iMax; i2++) {
                    const Real weight =
                        wg * coef[i2 - iMin][j2 - jMin][k2 - kMin];
#if (AMREX_SPACEDIM == 2)
                    const int idx0 = jp * 27 + (i2 - i1 + 1) * 9;
#else
                    const int idx0 = kp * 81 + jp * 27 + (i2 - i1 + 1) * 9;
#endif

                    Real* const data = &(data0[idx0]);
                    for (int idx = 0; idx < 9; idx++) {
                      deposit_add(&data[idx], alpha[idx] * weight);
                    }
                  } // k2

                } // j2
              } // if (ip > 0)
            } // i2
          } // k1

      //----- Mass matrix calculation end--------------
    });
  }

#if (AMREX_SPACEDIM != 2)
  for (MFIter mfi(nodeMM); mfi.isValid(); ++mfi) {
    // Finalize the mass matrix calculation: copy the kp = 2 half that the
    // particles deposited into the mirrored kp = 0 entries of the receiving
    // nodes. Every receiving entry has exactly one sender, so the nodes can
    // be processed in parallel.
    const Box box = mfi.validbox();
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Array4<RealMM> const mmArr = nodeMM[mfi].array();

    // We only need the mass matrix on the physical nodes. But the first layer
    // of the ghost nodes may contributes to the physical nodes below (ghost
    // node constributes as a sender). So, we need the '-1' and '+1' staff.
    const int iMin = lo.x - 1, jMin = lo.y - 1, kMin = nDim > 2 ? lo.z - 1 : 0;
    const int iMax = hi.x + 1, jMax = hi.y + 1, kMax = nDim > 2 ? hi.z + 1 : 0;
    const Box senders(IntVect(AMREX_D_DECL(iMin, jMin, kMin)),
                      IntVect(AMREX_D_DECL(iMax, jMax, kMax)));

    ParallelFor(senders, [=] AMREX_GPU_DEVICE(int i1, int j1, int k1) {
      const int kp = 2;
      const int kr = nDim > 2 ? k1 + kp - 1 : 0;
      if (kr > kMax || kr < kMin)
        return;
      const auto& datas0 = mmArr(i1, j1, k1);
      for (int jp = 0; jp < 3; jp++) {
        const int jr = j1 + jp - 1;
        if (jr > jMax || jr < jMin)
          continue;
        const int jpr = 2 - jp;
        for (int ip = 0; ip < 3; ip++) {
          const int ir = i1 + ip - 1;
          if (ir > iMax || ir < iMin)
            continue;
          const int ipr = 2 - ip;
          const int gpr = jpr * 3 + ipr;    // gp_receive
          const int gps = 18 + jp * 3 + ip; // gp_send = kp*9+jp*3+ip

          Real* const datar = &(mmArr(ir, jr, kr)[gpr * 9]);
          const Real* const datas = &(datas0[gps * 9]);
          for (int idx = 0; idx < 9; idx++) {
            datar[idx] = datas[idx];
          } // idx
        } // ip
      } // jp
    });
  }
#endif
}
//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::calc_mass_matrix_amr(
    NodeMMFab& nodeMM, amrex::Vector<amrex::Vector<NodeMMFab> >& nmmc,
    amrex::Vector<NodeMMFab>& nmmf, MultiFab& jHat,
    amrex::Vector<amrex::Vector<amrex::MultiFab> >& jhc,
    amrex::Vector<amrex::MultiFab>& jhf, MultiFab& nodeBMF,
    const MultiFab* nodeB0MF, MultiFab& u0MF, Real dt, int iLev,
    bool solveInCoMov, const amrex::Vector<amrex::iMultiFab>& cellstatus) {
  timing_func("Pts::calc_mass_matrix");

  Real qdto2mc = charge / mass * 0.5 * dt;
  const bool useB0 = (nodeB0MF != nullptr);

  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real const> const nodeBArr = nodeBMF[pti].const_array();
    const Array4<Real const> nodeB0Arr =
        useB0 ? (*nodeB0MF)[pti].const_array() : Array4<Real const>();
    Array4<Real> const jArr = jHat[pti].array();
    Array4<RealMM> const mmArr = nodeMM[pti].array();
    Array4<Real const> const u0Arr = u0MF[pti].const_array();
    // The particles also deposit onto the finer level if their cell is next
    // to a refined cell (a tile is a single cell on the CPU, a box on the
    // GPU, so the status is read at the particle's cell).
    const Array4<int const> status = cellstatus[iLev][pti].const_array();
    // jhf/nmmf exist below the finest level.
    const bool hasFiner = iLev < n_lev() - 1;
    constexpr int maxCoef = 8;
    AMREX_ALWAYS_ASSERT(iLev + 2 <= maxCoef);

    // Coarser levels, this level, and (possibly) the finer level.
    Array4<Real> jArrt[maxCoef];
    Array4<RealMM> mmArrt[maxCoef];
    Box mmBoxt[maxCoef]; // the (nodal) boxes of the mass matrices
    int count = 0;
    if (iLev > 0) {
      for (int i = 0; i < iLev; i++) {
        jArrt[count] = jhc[iLev][i][pti].array();
        mmArrt[count] = nmmc[iLev][i][pti].array();
        mmBoxt[count] = nmmc[iLev][i][pti].box();
        count++;
      }
    }
    jArrt[count] = jArr;
    mmArrt[count] = mmArr;
    mmBoxt[count] = nodeMM[pti].box();
    count++;

    if (hasFiner) {
      jArrt[count] = jhf[iLev][pti].array();
      mmArrt[count] = nmmf[iLev][pti].array();
      mmBoxt[count] = nmmf[iLev][pti].box();
    }

    // Copies, so that the kernels capture values, not host pointers.
    const int nCoefMax = hasFiner ? iLev + 2 : iLev + 1;
    Real ploAll[maxCoef][AMREX_SPACEDIM], invDxAll[maxCoef][AMREX_SPACEDIM];
    for (int i = 0; i < nCoefMax; i++) {
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        ploAll[i][d] = Geom(i).ProbLo(d);
        invDxAll[i][d] = Geom(i).InvCellSize(d);
      }
    }
    const Real invVolLev = invVol[iLev];

    const ParticleType* const pData = pti.GetArrayOfStructs()().data();

    const Dim3 lo = init_dim3(0);
    const Dim3 hi = init_dim3(1);

    // The state of particle ip: its stencils on the nCoef levels it deposits
    // onto, alpha and its current. False for a deleted particle.
    auto particle_state =
        [=] AMREX_GPU_DEVICE(int ip, IntVect(&loIdx)[maxCoef],
                             RealVect(&dShift)[maxCoef],
                             Real(&coef)[maxCoef][2][2][2], Real(&alpha)[9],
                             Real(&currents)[3], int& nCoef) -> bool {
      const ParticleType& p = pData[ip];
      if (p.id() < 0)
        return false;

      const bool refinedneighbour =
          hasFiner && is_next_to_refined(p, status, ploAll[iLev],
                                         invDxAll[iLev]);
      nCoef = iLev + 1 + refinedneighbour;
      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      for (int i = 0; i < nCoef; i++) {
        find_node_interpolation(p.pos(), ploAll[i], invDxAll[i], loIdx[i],
                                dShift[i], coef[i]);
      }

      //-----calculate interpolate coef end-------------

      //----- Mass matrix calculation begin--------------
      Real u0[3] = { 0, 0, 0 };
      Real bp[3] = { 0, 0, 0 };

      if (solveInCoMov && useB0) {
        const Array4<Real const> fields[3] = { nodeBArr, u0Arr, nodeB0Arr };
        Real b0[3] = { 0, 0, 0 };
        Real* values[3] = { bp, u0, b0 };
        interpolate_vector_fields(fields, loIdx[iLev], coef[iLev], lo, hi,
                                  values);
        bp[ix_] += b0[ix_];
        bp[iy_] += b0[iy_];
        bp[iz_] += b0[iz_];
      } else if (solveInCoMov) {
        const Array4<Real const> fields[2] = { nodeBArr, u0Arr };
        Real* values[2] = { bp, u0 };
        interpolate_vector_fields(fields, loIdx[iLev], coef[iLev], lo, hi,
                                  values);
      } else if (useB0) {
        const Array4<Real const> fields[2] = { nodeBArr, nodeB0Arr };
        Real b0[3] = { 0, 0, 0 };
        Real* values[2] = { bp, b0 };
        interpolate_vector_fields(fields, loIdx[iLev], coef[iLev], lo, hi,
                                  values);
        bp[ix_] += b0[ix_];
        bp[iy_] += b0[iy_];
        bp[iz_] += b0[iz_];
      } else {
        interpolate_vector_field(nodeBArr, loIdx[iLev], coef[iLev], lo, hi, bp);
      }

      const Real omx = qdto2mc * bp[ix_];
      const Real omy = qdto2mc * bp[iy_];
      const Real omz = qdto2mc * bp[iz_];

      // end interpolation
      const Real omx2 = omx * omx;
      const Real omy2 = omy * omy;
      const Real omz2 = omz * omz;
      const Real omxomy = omx * omy;
      const Real omxomz = omx * omz;
      const Real omyomz = omy * omz;
      const Real omsq = omx2 + omy2 + omz2;
      const Real denom = 1.0 / (1.0 + omsq);

      const Real c0 = denom * invVolLev * qp * qdto2mc;

      alpha[0] = (1 + omx2) * c0;
      alpha[1] = (omz + omxomy) * c0;
      alpha[2] = (-omy + omxomz) * c0;
      alpha[3] = (-omz + omxomy) * c0;
      alpha[4] = (1 + omy2) * c0;
      alpha[5] = (omx + omyomz) * c0;
      alpha[6] = (omy + omxomz) * c0;
      alpha[7] = (-omx + omyomz) * c0;
      alpha[8] = (1 + omz2) * c0;

      // jHat

      const Real up1 = up - u0[0];
      const Real vp1 = vp - u0[1];
      const Real wp1 = wp - u0[2];

      const Real udotOm1 = up1 * omx + vp1 * omy + wp1 * omz;

      const Real coef1 = denom * qp;
      currents[ix_] = (up1 + (vp1 * omz - wp1 * omy + udotOm1 * omx)) * coef1;
      currents[iy_] = (vp1 + (wp1 * omx - up1 * omz + udotOm1 * omy)) * coef1;
      currents[iz_] = (wp1 + (up1 * omy - vp1 * omx + udotOm1 * omz)) * coef1;

      return true;
    };

    const int np = pti.numParticles();

    if (useCellBinDeposits && AMREX_SPACEDIM == 3) {
      // GPU: deterministic per-cell sums (see deposit_mass_matrix_by_cell),
      // one level at a time.
      for (int iT = 0; iT < nCoefMax; ++iT) {
        GpuArray<Real, 3> ploT = { 0, 0, 0 }, invDxT = { 0, 0, 0 };
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          ploT[d] = ploAll[iT][d];
          invDxT[d] = invDxAll[iT][d];
        }
        const bool finer = iT == iLev + 1;
        deposit_mass_matrix_by_cell(
            np, enclosedCells(mmBoxt[iT]), jArrt[iT], mmArrt[iT],
            [=] AMREX_GPU_DEVICE(int ip, IntVect& cell) {
              const ParticleType& p = pData[ip];
              if (p.id() < 0 ||
                  (finer && !is_next_to_refined(p, status, ploAll[iLev],
                                                invDxAll[iLev])))
                return false;
              RealVect dS;
              Real c[2][2][2];
              find_node_interpolation(p.pos(), ploT.data(), invDxT.data(), cell,
                                      dS, c);
              return true;
            },
            [=] AMREX_GPU_DEVICE(int ip) { return particle_key(pData[ip]); },
            [=] AMREX_GPU_DEVICE(int ip, IntVect& loIdx1, RealVect& dShift1,
                                 Real(&coef1)[2][2][2], Real(&alpha)[9],
                                 Real(&currents)[3]) {
              IntVect loIdx[maxCoef];
              RealVect dShift[maxCoef];
              Real coef[maxCoef][2][2][2];
              int nCoef = 0;
              if (!particle_state(ip, loIdx, dShift, coef, alpha, currents,
                                  nCoef) ||
                  iT >= nCoef)
                return false;
              loIdx1 = loIdx[iT];
              dShift1 = dShift[iT];
              for (int kk = 0; kk < 2; ++kk)
                for (int jj = 0; jj < 2; ++jj)
                  for (int ii = 0; ii < 2; ++ii)
                    coef1[ii][jj][kk] = coef[iT][ii][jj][kk];
              return true;
            },
            /*kpMin=*/0);
      }
      continue;
    }

    ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) {
      IntVect loIdx[maxCoef];
      RealVect dShift[maxCoef];
      Real coef[maxCoef][2][2][2];
      Real alpha[9];
      Real currents[3];
      int nCoef = 0;
      if (!particle_state(ip, loIdx, dShift, coef, alpha, currents, nCoef))
        return;

      for (int i = 0; i < nCoef; ++i)
        deposit_vector_field(jArrt[i], loIdx[i], coef[i], lo, hi, currents, 3);

      for (int i = 0; i < nCoef; i++) {
        const int iMin = loIdx[i][ix_];
        const int jMin = loIdx[i][iy_];
        const int kMin = nDim > 2 ? loIdx[i][iz_] : 0;
        const int iMax = iMin + 1;
        const int jMax = jMin + 1;
        const int kMax = nDim > 2 ? kMin + 1 : 0;

        for (int k1 = kMin; k1 <= kMax; k1++)
          for (int j1 = jMin; j1 <= jMax; j1++)
            for (int i1 = iMin; i1 <= iMax; i1++) {
              const Real wg = coef[i][i1 - iMin][j1 - jMin][k1 - kMin];
              auto& data0 = mmArrt[i](i1, j1, k1);
              for (int k2 = kMin; k2 <= kMax; k2++) {
#if (AMREX_SPACEDIM == 3)
                const int kp = k2 - k1 + 1;
#endif
                // if (kp > 0)
                if (true) {
                  for (int j2 = jMin; j2 <= jMax; j2++) {
                    const int jp = j2 - j1 + 1;
                    for (int i2 = iMin; i2 <= iMax; i2++) {
                      const Real weight =
                          wg * coef[i][i2 - iMin][j2 - jMin][k2 - kMin];
#if (AMREX_SPACEDIM == 2)
                      const int idx0 = jp * 27 + (i2 - i1 + 1) * 9;
#else
                      const int idx0 = kp * 81 + jp * 27 + (i2 - i1 + 1) * 9;
#endif

                      Real* const data = &(data0[idx0]);
                      for (int idx = 0; idx < 9; idx++) {
                        deposit_add(&data[idx], alpha[idx] * weight);
                      }
                    } // k2

                  } // j2
                } // if (ip > 0)
              } // i2
            } // k1
      }
    });
  }
}
//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::calc_jhat(MultiFab& jHat,
                                                   MultiFab& nodeBMF,
                                                   const MultiFab* nodeB0MF,
                                                   Real dt) {
  timing_func("Pts::calc_jhat");

  Real qdto2mc = charge / mass * 0.5 * dt;
  const bool useB0 = (nodeB0MF != nullptr);

  const int iLev = 0;
  const auto ploLev = Geom(iLev).ProbLoArray();
  const auto invDxLev = Geom(iLev).InvCellSizeArray();
  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real const> const nodeBArr = nodeBMF[pti].const_array();
    const Array4<Real const> nodeB0Arr =
        useB0 ? (*nodeB0MF)[pti].const_array() : Array4<Real const>();
    Array4<Real> const jArr = jHat[pti].array();

    const ParticleType* const pData = pti.GetArrayOfStructs()().data();

    ParallelFor(pti.numParticles(), [=] AMREX_GPU_DEVICE(int ip) {
      const ParticleType& p = pData[ip];
      if (p.id() < 0)
        return;

      // Print()<<"p = "<<p<<std::endl;
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

      Real Bxl = 0, Byl = 0, Bzl = 0; // should be bp[3];

      for (int kk = 0; kk < 2; ++kk)
        for (int jj = 0; jj < 2; ++jj)
          for (int ii = 0; ii < 2; ++ii) {
            Bxl += nodeBArr(loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk,
                            ix_) *
                   coef[ii][jj][kk];
            Byl += nodeBArr(loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk,
                            iy_) *
                   coef[ii][jj][kk];
            Bzl += nodeBArr(loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk,
                            iz_) *
                   coef[ii][jj][kk];
          }

      if (useB0) {
        // Same gather for the frozen intrinsic field.
        for (int kk = 0; kk < 2; ++kk)
          for (int jj = 0; jj < 2; ++jj)
            for (int ii = 0; ii < 2; ++ii) {
              Bxl += nodeB0Arr(loIdx[ix_] + ii, loIdx[iy_] + jj,
                               loIdx[iz_] + kk, ix_) *
                     coef[ii][jj][kk];
              Byl += nodeB0Arr(loIdx[ix_] + ii, loIdx[iy_] + jj,
                               loIdx[iz_] + kk, iy_) *
                     coef[ii][jj][kk];
              Bzl += nodeB0Arr(loIdx[ix_] + ii, loIdx[iy_] + jj,
                               loIdx[iz_] + kk, iz_) *
                     coef[ii][jj][kk];
            }
      }

      const Real omx = qdto2mc * Bxl;
      const Real omy = qdto2mc * Byl;
      const Real omz = qdto2mc * Bzl;

      // end interpolation
      const Real omsq = (omx * omx + omy * omy + omz * omz);
      const Real denom = 1.0 / (1.0 + omsq);
      const Real udotOm = up * omx + vp * omy + wp * omz;

      {
        // jHat
        Real currents[3];

        {
          const Real coef1 = denom * qp;
          currents[ix_] = (up + (vp * omz - wp * omy + udotOm * omx)) * coef1;
          currents[iy_] = (vp + (wp * omx - up * omz + udotOm * omy)) * coef1;
          currents[iz_] = (wp + (up * omy - vp * omx + udotOm * omz)) * coef1;
        }

        deposit_vector_field(jArr, loIdx, coef, init_dim3(0), init_dim3(1),
                             currents, nDim);
      }
    });
  }
}

//==========================================================
// Apply mirror BC for node-centered current at reflect and inflow domain faces.
// For reflect faces, specular symmetry zeroes normal current and doubles
// tangential current. For inflow faces, all components double to compensate
// for half-space node weighting.

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::apply_jhat_mirror(MultiFab& jHat,
                                                           int iLev) {
  if (!Geom(iLev).isAllPeriodic() && jHat.nGrow() > 0) {
    const Box& dom = Geom(iLev).Domain();
    const int nCompJ = jHat.nComp();
    for (MFIter mfi(jHat); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.validbox();
      Array4<Real> const& jArr = jHat[mfi].array();

      for (int iDim = 0; iDim < nDim; ++iDim) {
        if (Geom(iLev).isPeriodic(iDim))
          continue;

        for (int side = 0; side < 2; ++side) {
          const bool isLo = (side == 0);
          const int faceBc = isLo ? bc.lo[iDim] : bc.hi[iDim];
          if (faceBc != ParticleBC::reflect && faceBc != ParticleBC::inflow)
            continue;

          // Upper boundary node is dom.bigEnd + 1 due to node-centering.
          const int domEdge = isLo ? dom.smallEnd(iDim) : dom.bigEnd(iDim) + 1;
          if ((isLo ? bx.smallEnd(iDim) : bx.bigEnd(iDim)) != domEdge)
            continue;

          const Box& fbx = mfi.fabbox();
          IntVect gs = fbx.smallEnd();
          IntVect ge = fbx.bigEnd();
          gs[iDim] = domEdge + (isLo ? -1 : 1);
          ge[iDim] = gs[iDim];
          const Box strip(gs, ge);

          const bool isReflect = (faceBc == ParticleBC::reflect);
          const int di = (iDim == 0) ? (isLo ? 1 : -1) : 0;
          const int dj = (iDim == 1) ? (isLo ? 1 : -1) : 0;
          const int dk = (iDim == 2) ? (isLo ? 1 : -1) : 0;

          ParallelFor(strip, nCompJ,
                      [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) {
                        const int ei = i + di, ej = j + dj, ek = k + dk;
                        if (isReflect && c == iDim) {
                          jArr(ei, ej, ek, c) = 0.0;
                        } else {
                          jArr(ei, ej, ek, c) *= 2.0;
                        }
                        jArr(i, j, k, c) = 0.0;
                      });
        }
      }
    }
  }
}

//==========================================================

#define INSTANTIATE_PARTICLES_MASS_MATRIX(T)                                   \
  template void T::calc_mass_matrix(NodeMMFab&, MultiFab&, MultiFab&,          \
                                    MultiFab const*, MultiFab&, Real, int,     \
                                    bool);                                     \
  template void T::calc_mass_matrix_amr(                                       \
      NodeMMFab&, Vector<Vector<NodeMMFab> >&, Vector<NodeMMFab>&, MultiFab&,  \
      Vector<Vector<MultiFab> >&, Vector<MultiFab>&, MultiFab&,                \
      MultiFab const*, MultiFab&, Real, int, bool, const Vector<iMultiFab>&);  \
  template void T::calc_jhat(MultiFab&, MultiFab&, MultiFab const*, Real);     \
  template void T::apply_jhat_mirror(MultiFab&, int);

INSTANTIATE_PARTICLES_MASS_MATRIX(PicParticles)
INSTANTIATE_PARTICLES_MASS_MATRIX(PTParticles)
#undef INSTANTIATE_PARTICLES_MASS_MATRIX

#include <cstdlib>

#include <AMReX_ParReduce.H>

#include "InitialCondition.h"
#include "Morton.h"
#include "Particles.h"
#include "SWMFInterface.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::sum_to_center(MultiFab& netChargeMF,
                                                       CenterMMFab& centerMM,
                                                       bool doNetChargeOnly,
                                                       int iLev) {
  timing_func("Pts::sum_to_center");

  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real> const& chargeArr = netChargeMF[pti].array();
    Array4<RealCMM> const& mmArr = centerMM[pti].array();
    const AoS& particles = pti.GetArrayOfStructs();

    const Dim3 lo = init_dim3(0);
    const Dim3 hi = init_dim3(1);

    for (const auto& p : particles) {
      /*
      Q: Why do not check p.id() < 0?
      A: IDs of ghost cell particles are set to -1 inside
      divE_correct_position(), but these particles should be take into account
      here.
      */

      // Print() << "particle = " << p << std::endl;

      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      IntVect loIdx;
      RealVect dShift;
      find_cell_index(p.pos(), Geom(iLev).ProbLo(), Geom(iLev).InvCellSize(),
                      loIdx, dShift);
      Real coef[2][2][2];
      linear_interpolation_coef(dShift, coef);
      //-----calculate interpolate coef end-------------

      const Real cTmp = qp * invVol[iLev];
      for (int kk = lo.z; kk <= hi.z; ++kk)
        for (int jj = lo.y; jj <= hi.y; ++jj)
          for (int ii = lo.x; ii <= hi.x; ++ii) {
            const IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + ii, loIdx[iy_] + jj,
                                               loIdx[iz_] + kk) };
            chargeArr(ijk) += coef[ii][jj][kk] * cTmp;
          }

      if (!doNetChargeOnly) {
        accumulate_mass_matrix_contribution(iLev, loIdx, dShift, qp, mmArr);
      } // if doChargeOnly

    } // for p
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
  for (int nLev = finer_level; nLev >= coarser_level; nLev--) {
    if (iLev == nLev) {
      for (PIter pti(*this, nLev); pti.isValid(); ++pti) {
        Array4<Real> const& chargeArr = netChargeMF[pti].array();

        Array4<RealCMM> const& mmArr = centerMM[pti].array();
        const Array4<int const>& status = cell_status(nLev)[pti].array();
        const AoS& particles = pti.GetArrayOfStructs();
        const Dim3 lo = init_dim3(0);
        const Dim3 hi = init_dim3(1);
        for (const auto& p : particles) {
          const Real qp = p.rdata(iqp_);
          IntVect loIdx;
          RealVect dShift;
          IntVect realIdx;
          RealVect tmprv;
          find_cell_index_exp(p.pos(), Geom(nLev).ProbLo(),
                              Geom(nLev).InvCellSize(), realIdx, tmprv);
          if (nLev == iLev || (bit::is_refined_neighbour(status(realIdx)) ||
                               bit::is_lev_edge(status(realIdx)))) {
            find_cell_index(p.pos(), Geom(iLev).ProbLo(),
                            Geom(iLev).InvCellSize(), loIdx, dShift);
            Real coef[2][2][2];
            linear_interpolation_coef(dShift, coef);
            const Real cTmp = qp * invVol[iLev];
            for (int kk = lo.z; kk <= hi.z; ++kk)
              for (int jj = lo.y; jj <= hi.y; ++jj)
                for (int ii = lo.x; ii <= hi.x; ++ii) {
                  const IntVect ijk = { AMREX_D_DECL(
                      loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk) };
                  chargeArr(ijk) += coef[ii][jj][kk] * cTmp;
                }
          }

          bool skipParticle = false;
          if (n_lev() > 1) {
            skipParticle = skip_particle_for_dive_cleaning(p.pos(), Geom(iLev),
                                                           iLev, status);
          }
          if (!doNetChargeOnly && !skipParticle) {
            accumulate_mass_matrix_contribution(iLev, loIdx, dShift, qp, mmArr);
          }
        }
      }
    }
    if (nLev > iLev) {
      for (PIter pti(*this, nLev); pti.isValid(); ++pti) {
        Array4<Real> const& chargeArr = jf[pti].array();
        const Array4<int const>& status = cell_status(nLev)[pti].array();
        const AoS& particles = pti.GetArrayOfStructs();
        const Dim3 lo = init_dim3(0);
        const Dim3 hi = init_dim3(1);
        for (const auto& p : particles) {
          const Real qp = p.rdata(iqp_);
          IntVect loIdx;
          RealVect dShift;
          IntVect realIdx;
          RealVect tmprv;
          find_cell_index_exp(p.pos(), Geom(nLev).ProbLo(),
                              Geom(nLev).InvCellSize(), realIdx, tmprv);
          if (nLev == iLev || (bit::is_refined_neighbour(status(realIdx)) ||
                               bit::is_lev_edge(status(realIdx)))) {
            find_cell_index(p.pos(), Geom(iLev).ProbLo(),
                            Geom(iLev).InvCellSize(), loIdx, dShift);
            Real coef[2][2][2];
            linear_interpolation_coef(dShift, coef);
            const Real cTmp = qp * invVol[iLev];
            for (int kk = lo.z; kk <= hi.z; ++kk)
              for (int jj = lo.y; jj <= hi.y; ++jj)
                for (int ii = lo.x; ii <= hi.x; ++ii) {
                  const IntVect ijk = { AMREX_D_DECL(
                      loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk) };

                  chargeArr(ijk) += coef[ii][jj][kk] * cTmp;
                }
          }
        }
      }
    }
    if (nLev < iLev) {
      for (PIter pti(*this, nLev); pti.isValid(); ++pti) {
        Array4<Real> const& chargeArr = jc[pti].array();
        const Array4<int const>& status = cell_status(nLev)[pti].array();
        const AoS& particles = pti.GetArrayOfStructs();
        const Dim3 lo = init_dim3(0);
        const Dim3 hi = init_dim3(1);
        for (const auto& p : particles) {
          const Real qp = p.rdata(iqp_);
          IntVect loIdx;
          RealVect dShift;
          IntVect realIdx;
          RealVect tmprv;
          find_cell_index_exp(p.pos(), Geom(nLev).ProbLo(),
                              Geom(nLev).InvCellSize(), realIdx, tmprv);
          if (nLev == iLev || (bit::is_refined_neighbour(status(realIdx)) ||
                               bit::is_lev_edge(status(realIdx)))) {
            find_cell_index(p.pos(), Geom(iLev).ProbLo(),
                            Geom(iLev).InvCellSize(), loIdx, dShift);
            Real coef[2][2][2];
            linear_interpolation_coef(dShift, coef);
            const Real cTmp = qp * invVol[iLev];
            for (int kk = lo.z; kk <= hi.z; ++kk)
              for (int jj = lo.y; jj <= hi.y; ++jj)
                for (int ii = lo.x; ii <= hi.x; ++ii) {
                  const IntVect ijk = { AMREX_D_DECL(
                      loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk) };

                  chargeArr(ijk) += coef[ii][jj][kk] * cTmp;
                }
          }
        }
      }
    }
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
std::array<Real, 5> Particles<NStructReal, NStructInt>::total_moments(
    bool localOnly) {
  timing_func("Pts::total_moments");

  std::array<Real, 5> sum = { 0, 0, 0, 0, 0 };

  for (int i = 0; i < 5; ++i)
    sum[i] = 0;

  const int iLev = 0;
  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    const AoS& particles = pti.GetArrayOfStructs();
    for (const auto& p : particles) {
      if (p.id() < 0)
        continue;

      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);

      sum[0] += qp;
      sum[1] += qp * up;
      sum[2] += qp * vp;
      sum[3] += qp * wp;
      sum[4] += 0.5 * qp * (up * up + vp * vp + wp * wp);
    }
  }

  for (int i = 0; i < 5; ++i)
    sum[i] *= qomSign * get_mass();

  if (!localOnly) {
    ParallelDescriptor::ReduceRealSum(sum.data(), sum.size(),
                                      ParallelDescriptor::IOProcessorNumber());
  }

  return sum;
}

//==========================================================

template <int NStructReal, int NStructInt>
Real Particles<NStructReal, NStructInt>::sum_moments(
    Vector<MultiFab>& momentsMF, Vector<MultiFab>& nodeBMF, Real dt) {
  timing_func("Pts::sum_moments");

  Real energy = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    timing_func("Pts::sum_moments_node_deposit");
    momentsMF[iLev].setVal(0.0);
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      Array4<Real> const& momentsArr = momentsMF[iLev][pti].array();

      const AoS& particles = pti.GetArrayOfStructs();

      const Dim3 lo = init_dim3(0);
      const Dim3 hi = init_dim3(1);

      // Print() << "iLev = " << iLev << std::endl;
      for (const auto& p : particles) {
        if (p.id() < 0)
          continue;

        // Print() << "p = " << p << std::endl;
        const Real up = p.rdata(iup_);
        const Real vp = p.rdata(ivp_);
        const Real wp = p.rdata(iwp_);
        const Real qp = p.rdata(iqp_);

        //-----calculate interpolate coef begin-------------
        IntVect loIdx;
        RealVect dShift;
        find_node_index(p.pos(), Geom(iLev).ProbLo(), Geom(iLev).InvCellSize(),
                        loIdx, dShift);
        Real coef[2][2][2];
        linear_interpolation_coef(dShift, coef);
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
          for (int kk = lo.z; kk <= hi.z; ++kk)
            for (int jj = lo.y; jj <= hi.y; ++jj)
              for (int ii = lo.x; ii <= hi.x; ++ii) {
                const IntVect ijk = { AMREX_D_DECL(
                    loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk) };
                momentsArr(ijk, iVar) += coef[ii][jj][kk] * pMoments[iVar];
              }

        //-------nodePlasma end---------

        energy += qp * (up * up + vp * vp + wp * wp);
      } // for p
    }

    // Exclude the number density.
    momentsMF[iLev].mult(invVol[iLev], 0, nMoments - 1,
                         momentsMF[iLev].nGrow());

    //----- Mirror boundary condition for the node-centred deposit ------------
    // A node-centred CIC deposit at a non-periodic wall leaves the boundary
    // node with only ~half the charge of an interior node: with the node at the
    // domain face, a particle in the edge cell deposits its weight between that
    // edge node and the next interior node, and NO particle lies on the
    // exterior side to contribute the complementary weight.  This is why, in
    // full-PIC mode, the reflecting-wall node of the 1D shock test showed rhoS0
    // ~ 6.25 instead of ~12.5, while the hybrid path (cell-centred deposit +
    // sum_moments_cell_centered) did not.
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
    if (!Geom(iLev).isAllPeriodic() && momentsMF[iLev].nGrow() > 0) {
      const Box& dom = Geom(iLev).Domain();
      const int nCompMF = momentsMF[iLev].nComp();
      for (MFIter mfi(momentsMF[iLev]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> const& arr = momentsMF[iLev][mfi].array();

        for (int iDim = 0; iDim < nDim; ++iDim) {
          if (Geom(iLev).isPeriodic(iDim))
            continue;

          for (int side = 0; side < 2; ++side) {
            const bool isLo = (side == 0);
            const int faceBc = isLo ? bc.lo[iDim] : bc.hi[iDim];
            // A particle face can no longer be `conducting`: that is a
            // field-only type, and on the particle side it is mapped to
            // `reflect` at parse time.  The fold set is therefore
            // reflect + inflow.
            const bool doFold =
                (faceBc == ParticleBC::reflect || faceBc == ParticleBC::inflow);
            if (!doFold)
              continue;

            // NOTE: momentsMF is NODE-centred (nGrids = convert(cGrids,
            // nodeVector)), whose domain box extends one node beyond the
            // cell-centred Geometry box on the hi side.  The edge NODE is
            // dom.smallEnd (lo) or dom.bigEnd+1 (hi); using dom.bigEnd would
            // hit the second-last node and skip the real edge tile.
            const int domEdge =
                isLo ? dom.smallEnd(iDim) : dom.bigEnd(iDim) + 1;
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

            ParallelFor(strip, nCompMF,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) {
                          // Odd parity components relative to iDim under
                          // specular reflection: normal momentum and
                          // off-diagonal normal-shear stresses.
                          bool isOdd = false;
                          if (iDim == 0) {
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
                  arr(i, j, k, c) = 0.0; // zero the ghost layer
#endif
                        });
          }
        }
      }
    }

    momentsMF[iLev].SumBoundary(Geom(iLev).periodicity());
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

  energy *= 0.5 * qomSign * get_mass();

  return energy;
}

//==========================================================
// Cell-centred moment deposit. The raw rho / momentum / pressure-tensor moments
// are scattered to the cell-centred momentsMF (centerPlasma[iSpecies]) with a
// plain cell-centred trilinear scatter (find_cell_index + linear weights).

template <int NStructReal, int NStructInt>
Real Particles<NStructReal, NStructInt>::sum_moments_cell_centered(
    Vector<MultiFab>& momentsMF) {
  timing_func("Pts::sum_moments_cell_centered");

  Real energy = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    momentsMF[iLev].setVal(0.0);
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      Array4<Real> const& momentsArr = momentsMF[iLev][pti].array();

      const AoS& particles = pti.GetArrayOfStructs();

      const Dim3 lo = init_dim3(0);
      const Dim3 hi = init_dim3(1);

      for (const auto& p : particles) {
        if (p.id() < 0)
          continue;

        const Real up = p.rdata(iup_);
        const Real vp = p.rdata(ivp_);
        const Real wp = p.rdata(iwp_);
        const Real qp = p.rdata(iqp_);

        //-----calculate interpolate coef begin-------------
        IntVect loIdx;
        RealVect dShift;
        // Cell-centred deposit: find the containing cell (find_cell_index) and
        // interpolate between its centre and the next cell centre (trilinear).
        find_cell_index(p.pos(), Geom(iLev).ProbLo(), Geom(iLev).InvCellSize(),
                        loIdx, dShift);
        Real coef[2][2][2];
        linear_interpolation_coef(dShift, coef);
        //-----calculate interpolate coef end-------------

        //-------cell-centred moments begin---------
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
          for (int kk = lo.z; kk <= hi.z; ++kk)
            for (int jj = lo.y; jj <= hi.y; ++jj)
              for (int ii = lo.x; ii <= hi.x; ++ii) {
                const IntVect ijk = { AMREX_D_DECL(
                    loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk) };
                momentsArr(ijk, iVar) += coef[ii][jj][kk] * pMoments[iVar];
              }
        //-------cell-centred moments end---------

        energy += qp * (up * up + vp * vp + wp * wp);
      } // for p
    }

    // Exclude the number density.
    momentsMF[iLev].mult(invVol[iLev], 0, nMoments - 1,
                         momentsMF[iLev].nGrow());

    // Fold first ghost layer back into edge cells at non-periodic domain
    // boundaries to conserve charge from outer half-cell CIC deposits.
    // Executed before SumBoundary so transverse ghost tails are correctly
    // combined across tile boundaries.
    if (!Geom(iLev).isAllPeriodic() && momentsMF[iLev].nGrow() > 0) {
      const Box& dom = Geom(iLev).Domain();
      const int nCompMF = momentsMF[iLev].nComp();
      for (MFIter mfi(momentsMF[iLev]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> const& arr = momentsMF[iLev][mfi].array();

        for (int iDim = 0; iDim < nDim; ++iDim) {
          if (Geom(iLev).isPeriodic(iDim))
            continue;

          for (int side = 0; side < 2; ++side) {
            const bool isLo = (side == 0);
            const int domEdge = isLo ? dom.smallEnd(iDim) : dom.bigEnd(iDim);
            if ((isLo ? bx.smallEnd(iDim) : bx.bigEnd(iDim)) != domEdge)
              continue;

            const Box& fbx = mfi.fabbox();
            IntVect gs = fbx.smallEnd();
            IntVect ge = fbx.bigEnd();
            gs[iDim] = domEdge + (isLo ? -1 : 1);
            ge[iDim] = gs[iDim];
            const Box strip(gs, ge);

            const int di = (iDim == 0) ? (isLo ? 1 : -1) : 0;
            const int dj = (iDim == 1) ? (isLo ? 1 : -1) : 0;
            const int dk = (iDim == 2) ? (isLo ? 1 : -1) : 0;

            ParallelFor(strip, nCompMF,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) {
                          arr(i + di, j + dj, k + dk, c) += arr(i, j, k, c);
                          arr(i, j, k, c) = 0.0;
                        });
          }
        }
      }
    }

    momentsMF[iLev].SumBoundary(Geom(iLev).periodicity());
  }

  // Cell-centred coarse-fine interface for AMR. Unlike the node-centred path,
  // we only need the fine→coarse summation (average_down + Add) for
  // cell-centred moments. The coarse→fine overwrite of level-edge cells is NOT
  // used here because it discards valid fine-level particle deposits (the fine
  // level-edge cells contain real particles), causing worse energy growth.
  for (int iLev = n_lev() - 2; iLev >= 0; iLev--) {
    timing_func("Pts::sum_moments_coarse_fine_interface_cell");
    sum_fine_to_coarse_lev_bny_cell(momentsMF[iLev], momentsMF[iLev + 1], 0,
                                    momentsMF[iLev].nComp(),
                                    get_ref_ratio(iLev));
  }

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

    for (MFIter mfi(momentsMF[iLev]); mfi.isValid(); ++mfi) {
      FArrayBox& fab = momentsMF[iLev][mfi];
      const Box& box = mfi.fabbox();
      const Array4<Real>& arr = fab.array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
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
          }
    }
  }
}

//==========================================================

// Explicit template instantiations.
template class Particles<nPicPartReal, nPicPartInt>;
template class Particles<nPTPartReal, nPTPartInt>;

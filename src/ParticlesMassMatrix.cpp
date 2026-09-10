#include <cstdlib>

#include <AMReX_ParReduce.H>

#include "InitialCondition.h"
#include "Morton.h"
#include "Particles.h"
#include "SWMFInterface.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::accumulate_mass_matrix_contribution(
    int iLev, const IntVect& loIdx, const RealVect& dShift, Real qp,
    Array4<RealCMM> const& mmArr) {

  Real weights_IIID[2][2][2][nDim3];
  //----- Mass matrix calculation begin--------------
  const Real xi0 = dShift[ix_] * dx[iLev][ix_];
  const Real eta0 = dShift[iy_] * dx[iLev][iy_];
  const Real zeta0 = nDim > 2 ? dShift[iz_] * dx[iLev][iz_] : 0;
  const Real xi1 = dx[iLev][ix_] - xi0;
  const Real eta1 = dx[iLev][iy_] - eta0;
  const Real zeta1 = nDim > 2 ? dx[iLev][iz_] - zeta0 : 1;

  weights_IIID[1][1][1][ix_] = eta0 * zeta0 * invVol[iLev];
  weights_IIID[1][1][1][iy_] = xi0 * zeta0 * invVol[iLev];
  weights_IIID[1][1][1][iz_] = xi0 * eta0 * invVol[iLev];

  // xi0*eta0*zeta1*invVol[iLev];
  weights_IIID[1][1][0][ix_] = eta0 * zeta1 * invVol[iLev];
  weights_IIID[1][1][0][iy_] = xi0 * zeta1 * invVol[iLev];
  weights_IIID[1][1][0][iz_] = -xi0 * eta0 * invVol[iLev];

  // xi0*eta1*zeta0*invVol[iLev];
  weights_IIID[1][0][1][ix_] = eta1 * zeta0 * invVol[iLev];
  weights_IIID[1][0][1][iy_] = -xi0 * zeta0 * invVol[iLev];
  weights_IIID[1][0][1][iz_] = xi0 * eta1 * invVol[iLev];

  // xi0*eta1*zeta1*invVol[iLev];
  weights_IIID[1][0][0][ix_] = eta1 * zeta1 * invVol[iLev];
  weights_IIID[1][0][0][iy_] = -xi0 * zeta1 * invVol[iLev];
  weights_IIID[1][0][0][iz_] = -xi0 * eta1 * invVol[iLev];

  // xi1*eta0*zeta0*invVol[iLev];
  weights_IIID[0][1][1][ix_] = -eta0 * zeta0 * invVol[iLev];
  weights_IIID[0][1][1][iy_] = xi1 * zeta0 * invVol[iLev];
  weights_IIID[0][1][1][iz_] = xi1 * eta0 * invVol[iLev];

  // xi1*eta0*zeta1*invVol[iLev];
  weights_IIID[0][1][0][ix_] = -eta0 * zeta1 * invVol[iLev];
  weights_IIID[0][1][0][iy_] = xi1 * zeta1 * invVol[iLev];
  weights_IIID[0][1][0][iz_] = -xi1 * eta0 * invVol[iLev];

  // xi1*eta1*zeta0*invVol[iLev];
  weights_IIID[0][0][1][ix_] = -eta1 * zeta0 * invVol[iLev];
  weights_IIID[0][0][1][iy_] = -xi1 * zeta0 * invVol[iLev];
  weights_IIID[0][0][1][iz_] = xi1 * eta1 * invVol[iLev];

  // xi1*eta1*zeta1*invVol[iLev];
  weights_IIID[0][0][0][ix_] = -eta1 * zeta1 * invVol[iLev];
  weights_IIID[0][0][0][iy_] = -xi1 * zeta1 * invVol[iLev];
  weights_IIID[0][0][0][iz_] = -xi1 * eta1 * invVol[iLev];

  const int iMin = loIdx[ix_];
  const int jMin = loIdx[iy_];
  const int kMin = nDim > 2 ? loIdx[iz_] : 0;
  const int iMax = iMin + 1;
  const int jMax = jMin + 1;
  const int kMax = nDim > 2 ? kMin + 1 : 0;

  const Real coef = fabs(qp) * invVol[iLev];
  Real wg_D[nDim3];
  for (int k1 = kMin; k1 <= kMax; k1++)
    for (int j1 = jMin; j1 <= jMax; j1++)
      for (int i1 = iMin; i1 <= iMax; i1++) {

        for (int iDim = 0; iDim < nDim; iDim++) {
          wg_D[iDim] =
              coef * weights_IIID[i1 - iMin][j1 - jMin][k1 - kMin][iDim];
        }

        auto& data = mmArr(i1, j1, k1);
        // Real weights[27] = { 0 };
        for (int i2 = iMin; i2 <= iMax; i2++) {
          int ip = i2 - i1 + 1;
          const int gp0 = ip * 9;
          for (int j2 = jMin; j2 <= jMax; j2++) {
            int jp = j2 - j1 + 1;
            const int gp1 = gp0 + jp * nDim3;
            for (int k2 = kMin; k2 <= kMax; k2++) {
              const Real(&wg1_D)[nDim3] =
                  weights_IIID[i2 - iMin][j2 - jMin][k2 - kMin];

              // const int kp = k2 - k1 + 1;
              const int gp = gp1 + k2 - k1 + 1;
              for (int iDim = 0; iDim < nDim; iDim++) {
                data[gp] += wg_D[iDim] * wg1_D[iDim];
              }
            }
          }
        }
      }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::calc_mass_matrix(
    NodeMMFab& nodeMM, MultiFab& jHat, MultiFab& nodeBMF, MultiFab& u0MF,
    Real dt, int iLev, bool solveInCoMov) {
  timing_func("Pts::calc_mass_matrix");

  Real qdto2mc = charge / mass * 0.5 * dt;

  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real const> const& nodeBArr = nodeBMF[pti].array();
    Array4<Real> const& jArr = jHat[pti].array();
    Array4<RealMM> const& mmArr = nodeMM[pti].array();

    Array4<Real const> const& u0Arr = u0MF[pti].array();

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

      //-----calculate interpolation coef begin-------------
      IntVect loIdx;
      RealVect dShift;

      find_node_index(p.pos(), Geom(iLev).ProbLo(), Geom(iLev).InvCellSize(),
                      loIdx, dShift);

      Real coef[2][2][2];
      linear_interpolation_coef(dShift, coef);
      //-----calculate interpolate coef end-------------

      //----- Mass matrix calculation begin--------------
      Real u0[3] = { 0, 0, 0 };
      Real bp[3] = { 0, 0, 0 };

      for (int kk = lo.z; kk <= hi.z; ++kk)
        for (int jj = lo.y; jj <= hi.y; ++jj)
          for (int ii = lo.x; ii <= hi.x; ++ii) {
            const IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + ii, loIdx[iy_] + jj,
                                               loIdx[iz_] + kk) };
            for (int iDim = 0; iDim < nDim3; iDim++) {
              bp[iDim] += nodeBArr(ijk, iDim) * coef[ii][jj][kk];

              if (solveInCoMov)
                u0[iDim] += u0Arr(ijk, iDim) * coef[ii][jj][kk];
            }
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

      const Real c0 = denom * invVol[iLev] * qp * qdto2mc;

      Real alpha[9];
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
        Real currents[3];

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

        for (int iVar = 0; iVar < 3; iVar++)
          for (int kk = lo.z; kk <= hi.z; ++kk)
            for (int jj = lo.y; jj <= hi.y; ++jj)
              for (int ii = lo.x; ii <= hi.x; ++ii) {
                IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + ii, loIdx[iy_] + jj,
                                             loIdx[iz_] + kk) };
                jArr(ijk, iVar) += coef[ii][jj][kk] * currents[iVar];
              }
      }

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
                    const int idx0 = kp * 81 + jp * 27 + (i2 - i1 + 1) * 9;

                    Real* const data = &(data0[idx0]);
                    for (int idx = 0; idx < 9; idx++) {
                      data[idx] += alpha[idx] * weight;
                    }
                  } // k2

                } // j2
              } // if (ip > 0)
            } // i2
          } // k1

      //----- Mass matrix calculation end--------------

    } // for p
  }

  for (MFIter mfi(nodeMM); mfi.isValid(); ++mfi) {
    // Finalize the mass matrix calculation.
    const Box box = mfi.validbox();
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Array4<RealMM> const& mmArr = nodeMM[mfi].array();

    // We only need the mass matrix on the physical nodes. But the first layer
    // of the ghost nodes may contributes to the physical nodes below (ghost
    // node constributes as a sender). So, we need the '-1' and '+1' staff.
    const int iMin = lo.x - 1, jMin = lo.y - 1, kMin = nDim > 2 ? lo.z - 1 : 0;
    const int iMax = hi.x + 1, jMax = hi.y + 1, kMax = nDim > 2 ? hi.z + 1 : 0;

    int gps, gpr; // gp_send, gp_receive
    for (int k1 = kMin; k1 <= kMax; k1++)
      for (int j1 = jMin; j1 <= jMax; j1++)
        for (int i1 = iMin; i1 <= iMax; i1++) {
          const int kp = 2;
          const int kr = nDim > 2 ? k1 + kp - 1 : 0;
          if (kr > kMax || kr < kMin)
            continue;
          auto& datas0 = mmArr(i1, j1, k1);
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
              gpr = jpr * 3 + ipr;
              gps = 18 + jp * 3 + ip; // gps = kp*9+jp*3+kp

              Real* const datar = &(mmArr(ir, jr, kr)[gpr * 9]);
              const Real* const datas = &(datas0[gps * 9]);
              for (int idx = 0; idx < 9; idx++) {
                datar[idx] = datas[idx];
              } // idx
            } // kp
          } // jp
        } // k1
  }
}
//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::calc_mass_matrix_amr(
    NodeMMFab& nodeMM, amrex::Vector<amrex::Vector<NodeMMFab> >& nmmc,
    amrex::Vector<NodeMMFab>& nmmf, MultiFab& jHat,
    amrex::Vector<amrex::Vector<amrex::MultiFab> >& jhc,
    amrex::Vector<amrex::MultiFab>& jhf, MultiFab& nodeBMF, MultiFab& u0MF,
    Real dt, int iLev, bool solveInCoMov,
    amrex::Vector<amrex::iMultiFab>& cellstatus) {
  timing_func("Pts::calc_mass_matrix");

  Real qdto2mc = charge / mass * 0.5 * dt;

  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real const> const& nodeBArr = nodeBMF[pti].array();
    Array4<Real> const& jArr = jHat[pti].array();
    Array4<RealMM> const& mmArr = nodeMM[pti].array();
    Array4<Real const> const& u0Arr = u0MF[pti].array();
    const Array4<int const>& status = cellstatus[iLev][pti].array();
    Box bx = pti.tilebox();
    IntVect ibx = bx.smallEnd();
    bool refinedneighbour = false;
    if (bit::is_refined_neighbour(status(ibx))) {
      refinedneighbour = true;
    }
    amrex::Vector<Array4<Real> > jArrt;
    amrex::Vector<Array4<RealMM> > mmArrt;
    if (iLev > 0) {
      for (int i = 0; i < iLev; i++) {
        jArrt.push_back(jhc[iLev][i][pti].array());
        mmArrt.push_back(nmmc[iLev][i][pti].array());
      }
    }
    jArrt.push_back(jArr);
    mmArrt.push_back(mmArr);

    if (refinedneighbour) {
      jArrt.push_back(jhf[iLev][pti].array());
      mmArrt.push_back(nmmf[iLev][pti].array());
    }

    amrex::Vector<IntVect> loIdx;
    amrex::Vector<RealVect> dShift;
    const int nCoef = iLev + 1 + refinedneighbour;
    loIdx.resize(nCoef);
    dShift.resize(nCoef);
    using InterpolationCoef = Real[2][2][2];
    auto coef = std::make_unique<InterpolationCoef[]>(nCoef);

    const AoS& particles = pti.GetArrayOfStructs();

    const Dim3 lo = init_dim3(0);
    const Dim3 hi = init_dim3(1);

    for (const auto& p : particles) {
      if (p.id() < 0)
        continue;

      // Print()<<"p = "<<p<<std::endl;
      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      for (int i = 0; i < nCoef; i++) {
        find_node_index(p.pos(), Geom(i).ProbLo(), Geom(i).InvCellSize(),
                        loIdx[i], dShift[i]);

        linear_interpolation_coef(dShift[i], coef[i]);
      }

      //-----calculate interpolate coef end-------------

      //----- Mass matrix calculation begin--------------
      Real u0[3] = { 0, 0, 0 };
      Real bp[3] = { 0, 0, 0 };

      for (int kk = lo.z; kk <= hi.z; ++kk)
        for (int jj = lo.y; jj <= hi.y; ++jj)
          for (int ii = lo.x; ii <= hi.x; ++ii) {
            const IntVect ijk = { AMREX_D_DECL(loIdx[iLev][ix_] + ii,
                                               loIdx[iLev][iy_] + jj,
                                               loIdx[iLev][iz_] + kk) };
            for (int iDim = 0; iDim < nDim3; iDim++) {
              bp[iDim] += nodeBArr(ijk, iDim) * coef[iLev][ii][jj][kk];

              if (solveInCoMov)
                u0[iDim] += u0Arr(ijk, iDim) * coef[iLev][ii][jj][kk];
            }
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

      const Real c0 = denom * invVol[iLev] * qp * qdto2mc;

      Real alpha[9];
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
      Real currents[3];

      const Real up1 = up - u0[0];
      const Real vp1 = vp - u0[1];
      const Real wp1 = wp - u0[2];

      const Real udotOm1 = up1 * omx + vp1 * omy + wp1 * omz;

      const Real coef1 = denom * qp;
      currents[ix_] = (up1 + (vp1 * omz - wp1 * omy + udotOm1 * omx)) * coef1;
      currents[iy_] = (vp1 + (wp1 * omx - up1 * omz + udotOm1 * omy)) * coef1;
      currents[iz_] = (wp1 + (up1 * omy - vp1 * omx + udotOm1 * omz)) * coef1;

      for (int iVar = 0; iVar < 3; iVar++)
        for (int kk = lo.z; kk <= hi.z; ++kk)
          for (int jj = lo.y; jj <= hi.y; ++jj)
            for (int ii = lo.x; ii <= hi.x; ++ii) {
              for (int i = 0; i < nCoef; i++) {
                IntVect ijk = { AMREX_D_DECL(loIdx[i][ix_] + ii,
                                             loIdx[i][iy_] + jj,
                                             loIdx[i][iz_] + kk) };
                jArrt[i](ijk, iVar) += coef[i][ii][jj][kk] * currents[iVar];
              }
            }

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
                const int kp = k2 - k1 + 1;
                // if (kp > 0)
                if (true) {
                  for (int j2 = jMin; j2 <= jMax; j2++) {
                    const int jp = j2 - j1 + 1;
                    for (int i2 = iMin; i2 <= iMax; i2++) {
                      const Real weight =
                          wg * coef[i][i2 - iMin][j2 - jMin][k2 - kMin];
                      const int idx0 = kp * 81 + jp * 27 + (i2 - i1 + 1) * 9;

                      Real* const data = &(data0[idx0]);
                      for (int idx = 0; idx < 9; idx++) {
                        data[idx] += alpha[idx] * weight;
                      }
                    } // k2

                  } // j2
                } // if (ip > 0)
              } // i2
            } // k1
      }
    } // for p
  }

  // for (MFIter mfi(nodeMM); mfi.isValid(); ++mfi) {
  //   // Finalize the mass matrix calculation.
  //   const Box box = mfi.validbox();
  //   const auto lo = lbound(box);
  //   const auto hi = ubound(box);

  //   Array4<RealMM> const& mmArr = nodeMM[mfi].array();

  //   // We only need the mass matrix on the physical nodes. But the first
  //   // layer
  //   // of the ghost nodes may contributes to the physical nodes below (ghost
  //   // node constributes as a sender). So, we need the '-1' and '+1' staff.
  //   const int iMin = lo.x - 1, jMin = lo.y - 1, kMin = nDim > 2 ? lo.z - 1 :
  //   0; const int iMax = hi.x + 1, jMax = hi.y + 1, kMax = nDim > 2 ? hi.z + 1
  //   : 0;

  //   int gps, gpr; // gp_send, gp_receive
  //   for (int k1 = kMin; k1 <= kMax; k1++)
  //     for (int j1 = jMin; j1 <= jMax; j1++)
  //       for (int i1 = iMin; i1 <= iMax; i1++) {
  //         const int kp = 2;
  //         const int kr = nDim > 2 ? k1 + kp - 1 : 0;
  //         if (kr > kMax || kr < kMin)
  //           continue;
  //         auto& datas0 = mmArr(i1, j1, k1);
  //         for (int jp = 0; jp < 3; jp++) {
  //           const int jr = j1 + jp - 1;
  //           if (jr > jMax || jr < jMin)
  //             continue;
  //           const int jpr = 2 - jp;
  //           for (int ip = 0; ip < 3; ip++) {
  //             const int ir = i1 + ip - 1;
  //             if (ir > iMax || ir < iMin)
  //               continue;
  //             const int ipr = 2 - ip;
  //             gpr = jpr * 3 + ipr;
  //             gps = 18 + jp * 3 + ip; // gps = kp*9+jp*3+kp

  //             Real* const datar = &(mmArr(ir, jr, kr)[gpr * 9]);
  //             const Real* const datas = &(datas0[gps * 9]);
  //             for (int idx = 0; idx < 9; idx++) {
  //               datar[idx] = datas[idx];
  //             } // idx
  //           } // kp
  //         } // jp
  //       } // k1
  // }
}
//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::calc_jhat(MultiFab& jHat,
                                                   MultiFab& nodeBMF, Real dt) {
  timing_func("Pts::calc_jhat");

  Real qdto2mc = charge / mass * 0.5 * dt;

  const int iLev = 0;
  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real const> const& nodeBArr = nodeBMF[pti].array();
    Array4<Real> const& jArr = jHat[pti].array();

    const AoS& particles = pti.GetArrayOfStructs();

    for (const auto& p : particles) {
      if (p.id() < 0)
        continue;

      // Print()<<"p = "<<p<<std::endl;
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

        for (int iVar = 0; iVar < nDim; iVar++)
          for (int kk = 0; kk < 2; ++kk)
            for (int jj = 0; jj < 2; ++jj)
              for (int ii = 0; ii < 2; ++ii) {
                jArr(loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk, iVar) +=
                    coef[ii][jj][kk] * currents[iVar];
              }
      }

    } // for p
  }
}

//==========================================================
// Apply mirror BC for node-centred current at reflect and inflow domain faces.
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

// Explicit template instantiations.
template class Particles<nPicPartReal, nPicPartInt>;
template class Particles<nPTPartReal, nPTPartInt>;

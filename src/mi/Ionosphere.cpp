#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include "GpuBackend.h"
#include "Ionosphere.h"

using namespace amrex;

namespace {

constexpr Real cDegToRad = dPI / 180.0;

inline int idx2d(int i, int j, int nPhi) { return i * nPhi + j; }

Real dot_vec(const Real *AMREX_RESTRICT a, const Real *AMREX_RESTRICT b,
             size_t n) {
  Real s = 0.0;
#pragma omp simd reduction(+ : s)
  for (size_t i = 0; i < n; ++i)
    s += a[i] * b[i];
  return s;
}

Real norm2_vec(const Real *a, size_t n) { return std::sqrt(dot_vec(a, a, n)); }

// Right-preconditioned BiCGSTAB for A x = b using persistent workspace.
bool solve_bicgstab(HemisphereWorkspace &hw, Real tol, int maxIter,
                    int &nIterOut, Real &residOut) {
  const size_t n = hw.rhs.size();
  if (n == 0) {
    nIterOut = 0;
    residOut = 0.0;
    return true;
  }

  Real *AMREX_RESTRICT x = hw.x.data();
  const Real *AMREX_RESTRICT b = hw.rhs.data();
  Real *AMREX_RESTRICT r = hw.r.data();
  Real *AMREX_RESTRICT rHat = hw.rHat.data();
  Real *AMREX_RESTRICT p = hw.p.data();
  Real *AMREX_RESTRICT v = hw.v.data();
  Real *AMREX_RESTRICT s = hw.s.data();
  Real *AMREX_RESTRICT t = hw.t.data();
  Real *AMREX_RESTRICT y = hw.y.data();
  Real *AMREX_RESTRICT z = hw.z.data();

  const Real bNorm = norm2_vec(b, n);
  if (bNorm <= Real(1.0e-30)) {
    std::fill(x, x + n, Real(0.0));
    nIterOut = 0;
    residOut = 0.0;
    return true;
  }

  hw.matvec(hw.x, hw.v);
#pragma omp simd
  for (size_t i = 0; i < n; ++i) {
    const Real ri = b[i] - v[i];
    r[i] = ri;
    rHat[i] = ri;
    p[i] = ri;
  }

  Real rNorm = norm2_vec(r, n);
  residOut = rNorm / bNorm;
  if (residOut <= tol) {
    nIterOut = 0;
    return true;
  }

  const Real invBNorm = Real(1.0) / bNorm;
  Real rhoOld = 1.0;
  Real alpha = 1.0;
  Real omega = 1.0;

  for (int iter = 1; iter <= maxIter; ++iter) {
    const Real rhoNew = dot_vec(rHat, r, n);
    if (std::abs(rhoNew) < Real(1.0e-30)) {
      nIterOut = iter;
      residOut = rNorm * invBNorm;
      return residOut <= tol;
    }

    if (iter > 1) {
      const Real beta = (rhoNew / rhoOld) * (alpha / omega);
      const Real negBetaOmega = -beta * omega;
#pragma omp simd
      for (size_t i = 0; i < n; ++i) {
        p[i] = r[i] + beta * p[i] + negBetaOmega * v[i];
      }
    }

    hw.apply_precond(hw.p, hw.y);
    hw.matvec(hw.y, hw.v);

    const Real rHatV = dot_vec(rHat, v, n);
    if (std::abs(rHatV) < Real(1.0e-30)) {
      nIterOut = iter;
      residOut = rNorm * invBNorm;
      return residOut <= tol;
    }
    alpha = rhoNew / rHatV;

    Real s2 = 0.0;
#pragma omp simd reduction(+ : s2)
    for (size_t i = 0; i < n; ++i) {
      const Real si = r[i] - alpha * v[i];
      s[i] = si;
      s2 += si * si;
    }

    const Real sNorm = std::sqrt(s2);
    if (sNorm * invBNorm <= tol) {
#pragma omp simd
      for (size_t i = 0; i < n; ++i)
        x[i] += alpha * y[i];
      nIterOut = iter;
      residOut = sNorm * invBNorm;
      return true;
    }

    hw.apply_precond(hw.s, hw.z);
    hw.matvec(hw.z, hw.t);

    // Fuse tt and ts dot products in a single pass.
    Real tt = 0.0;
    Real ts = 0.0;
#pragma omp simd reduction(+ : tt, ts)
    for (size_t i = 0; i < n; ++i) {
      const Real ti = t[i];
      tt += ti * ti;
      ts += ti * s[i];
    }

    if (tt < Real(1.0e-30)) {
#pragma omp simd
      for (size_t i = 0; i < n; ++i)
        x[i] += alpha * y[i];
      nIterOut = iter;
      residOut = sNorm * invBNorm;
      return residOut <= tol;
    }
    omega = ts / tt;

    // Fuse x update, r update, and ||r||^2 in a single pass.
    Real r2 = 0.0;
#pragma omp simd reduction(+ : r2)
    for (size_t i = 0; i < n; ++i) {
      x[i] += alpha * y[i] + omega * z[i];
      const Real ri = s[i] - omega * t[i];
      r[i] = ri;
      r2 += ri * ri;
    }

    rNorm = std::sqrt(r2);
    residOut = rNorm * invBNorm;
    nIterOut = iter;
    if (residOut <= tol)
      return true;

    if (std::abs(omega) < Real(1.0e-30))
      return false;
    rhoOld = rhoNew;
  }

  return residOut <= tol;
}

} // namespace

//==========================================================
void HemisphereWorkspace::allocate(int iMinIn, int iMaxIn, int nPhiUsed,
                                   bool northIn, bool neumannIn) {
  iMin = iMinIn;
  iMax = iMaxIn;
  nTh = (iMax >= iMin) ? (iMax - iMin + 1) : 0;
  nPh = nPhiUsed;
  isNorth = northIn;
  isNeumann = neumannIn;

  const size_t nX = static_cast<size_t>(nTh) * nPh;
  const int nDiag = (nTh > 0 && nPh > 0) ? (nTh + nPh - 1) : 0;
  diags.resize(static_cast<size_t>(nDiag));

  int curOff = 0;
  for (int d = 0; d < nDiag; ++d) {
    const int iStart = std::max(0, d - (nPh - 1));
    const int iEnd = std::min(nTh - 1, d);
    const int len = iEnd - iStart + 1;
    diags[d].off = curOff;
    diags[d].len = len;
    diags[d].iStart = iStart;
    diags[d].prevBase = 0;
    diags[d].nextBase = 0;
    curOff += len;
  }
  for (int d = 0; d < nDiag; ++d) {
    if (d > 0) {
      diags[d].prevBase =
          diags[d - 1].off + (diags[d].iStart - diags[d - 1].iStart);
    }
    if (d + 1 < nDiag) {
      diags[d].nextBase =
          diags[d + 1].off + (diags[d].iStart - diags[d + 1].iStart);
    }
  }

  poleWfIdx.resize(static_cast<size_t>(std::max(0, nPh)));
  if (nTh > 0 && nPh > 0) {
    const int iPole = isNorth ? 0 : (nTh - 1);
    for (int j = 0; j < nPh; ++j) {
      poleWfIdx[j] = wf_index(iPole, j);
    }
    kPin = wf_index(isNorth ? (nTh - 1) : 0, 0);
  } else {
    kPin = 0;
  }

  cA.assign(nX, 0.0);
  cB.assign(nX, 0.0);
  cC.assign(nX, 0.0);
  cD.assign(nX, 0.0);
  cE.assign(nX, 0.0);

  dInv.assign(nX, 0.0);
  eL.assign(nX, 0.0);
  fU.assign(nX, 0.0);
  fUs.assign(nX, 0.0);
  e1L.assign(nX, 0.0);
  f1U.assign(nX, 0.0);
  f1Us.assign(nX, 0.0);

  rhs.assign(nX, 0.0);
  x.assign(nX, 0.0);
  r.assign(nX, 0.0);
  rHat.assign(nX, 0.0);
  p.assign(nX, 0.0);
  v.assign(nX, 0.0);
  s.assign(nX, 0.0);
  t.assign(nX, 0.0);
  y.assign(nX, 0.0);
  z.assign(nX, 0.0);
}

//==========================================================
void HemisphereWorkspace::build_ilu0() {
  const int n = nTh * nPh;
  if (n <= 0)
    return;

  const int nDiag = static_cast<int>(diags.size());

  const Real *AMREX_RESTRICT a = cA.data();
  const Real *AMREX_RESTRICT b = cB.data();
  const Real *AMREX_RESTRICT c = cC.data();
  const Real *AMREX_RESTRICT dC = cD.data();
  const Real *AMREX_RESTRICT eC = cE.data();

  Real *AMREX_RESTRICT dI = dInv.data();
  Real *AMREX_RESTRICT e0 = eL.data();
  Real *AMREX_RESTRICT f0 = fU.data();
  Real *AMREX_RESTRICT f0s = fUs.data();
  Real *AMREX_RESTRICT e1 = e1L.data();
  Real *AMREX_RESTRICT f1 = f1U.data();
  Real *AMREX_RESTRICT f1s = f1Us.data();

  // Initialize triangular stencil weights in wavefront order with endpoint
  // peeling so interior diagonal loops are branchless SIMD.
  for (int d = 0; d < nDiag; ++d) {
    const WavefrontDiag &dg = diags[d];
    const int len = dg.len;
    const int off = dg.off;

    auto init_endpoint = [&](int m) noexcept {
      const int i = dg.iStart + m;
      const int j = d - i;
      const int k = off + m;
      e0[k] = (i > 0) ? b[k] : Real(0.0);
      f0[k] = (i < nTh - 1) ? c[k] : Real(0.0);
      e1[k] = (j > 0) ? dC[k] : Real(0.0);
      f1[k] = (j < nPh - 1) ? eC[k] : Real(0.0);
    };

    init_endpoint(0);
    if (len >= 3) {
      const Real *AMREX_RESTRICT bCur = b + off;
      const Real *AMREX_RESTRICT cCur = c + off;
      const Real *AMREX_RESTRICT dCur = dC + off;
      const Real *AMREX_RESTRICT eCur = eC + off;
      Real *AMREX_RESTRICT e0Cur = e0 + off;
      Real *AMREX_RESTRICT f0Cur = f0 + off;
      Real *AMREX_RESTRICT e1Cur = e1 + off;
      Real *AMREX_RESTRICT f1Cur = f1 + off;
#pragma omp simd
      for (int m = 1; m < len - 1; ++m) {
        e0Cur[m] = bCur[m];
        f0Cur[m] = cCur[m];
        e1Cur[m] = dCur[m];
        f1Cur[m] = eCur[m];
      }
    }
    if (len > 1) {
      init_endpoint(len - 1);
    }
  }

  // Wavefront MILU(0.5) factorization: all entries along diagonal d depend only
  // on diagonal d - 1, so interior entries along d vectorize cleanly.
  for (int d = 0; d < nDiag; ++d) {
    const WavefrontDiag &dg = diags[d];
    const int len = dg.len;
    const int off = dg.off;

    auto factor_endpoint = [&](int m) noexcept {
      const int i = dg.iStart + m;
      const int j = d - i;
      const int k = off + m;
      Real pivot = a[k];
      if (i > 0) {
        const int km = dg.prevBase - 1 + m;
        pivot -= e0[k] * (f0[km] + Real(0.5) * f1[km]) * dI[km];
      }
      if (j > 0) {
        const int km = dg.prevBase + m;
        pivot -= e1[k] * (f1[km] + Real(0.5) * f0[km]) * dI[km];
      }
      if (std::abs(pivot) < Real(1.0e-30))
        pivot = (pivot >= Real(0.0) ? Real(1.0e-30) : Real(-1.0e-30));
      const Real invP = Real(1.0) / pivot;
      dI[k] = invP;
      f0s[k] = f0[k] * invP;
      f1s[k] = f1[k] * invP;
    };

    factor_endpoint(0);
    if (len >= 3) {
      const Real *aCur = a + off;
      const Real *e0Cur = e0 + off;
      const Real *e1Cur = e1 + off;
      const Real *f0Cur = f0 + off;
      const Real *f1Cur = f1 + off;
      const Real *dPrevI = dI + (dg.prevBase - 1);
      const Real *f0PrevI = f0 + (dg.prevBase - 1);
      const Real *f1PrevI = f1 + (dg.prevBase - 1);
      const Real *dPrevJ = dI + dg.prevBase;
      const Real *f0PrevJ = f0 + dg.prevBase;
      const Real *f1PrevJ = f1 + dg.prevBase;
      Real *dICur = dI + off;
      Real *f0sCur = f0s + off;
      Real *f1sCur = f1s + off;

#pragma omp simd
      for (int m = 1; m < len - 1; ++m) {
        Real pivot =
            aCur[m] -
            e0Cur[m] * (f0PrevI[m] + Real(0.5) * f1PrevI[m]) * dPrevI[m] -
            e1Cur[m] * (f1PrevJ[m] + Real(0.5) * f0PrevJ[m]) * dPrevJ[m];
        pivot = (std::abs(pivot) < Real(1.0e-30))
                    ? (pivot >= Real(0.0) ? Real(1.0e-30) : Real(-1.0e-30))
                    : pivot;
        const Real invP = Real(1.0) / pivot;
        dICur[m] = invP;
        f0sCur[m] = f0Cur[m] * invP;
        f1sCur[m] = f1Cur[m] * invP;
      }
    }
    if (len > 1) {
      factor_endpoint(len - 1);
    }
  }
}

//==========================================================
void HemisphereWorkspace::apply_precond(const std::vector<Real> &in,
                                        std::vector<Real> &out) const {
  const int n = nTh * nPh;
  if (n <= 0)
    return;

  const int nDiag = static_cast<int>(diags.size());
  const Real *AMREX_RESTRICT rIn = in.data();
  Real *zOut = out.data();
  const Real *AMREX_RESTRICT dI = dInv.data();
  const Real *AMREX_RESTRICT e0 = eL.data();
  const Real *AMREX_RESTRICT e1 = e1L.data();
  const Real *AMREX_RESTRICT fu = fUs.data();
  const Real *AMREX_RESTRICT f1u = f1Us.data();

  // Forward wavefront solve: L y = r (d = 0 .. nDiag - 1).
  for (int d = 0; d < nDiag; ++d) {
    const WavefrontDiag &dg = diags[d];
    const int len = dg.len;
    const int off = dg.off;

    // Endpoint m = 0
    {
      const int i = dg.iStart;
      const int j = d - i;
      Real val = rIn[off];
      if (i > 0)
        val -= e0[off] * zOut[dg.prevBase - 1];
      if (j > 0)
        val -= e1[off] * zOut[dg.prevBase];
      zOut[off] = val * dI[off];
    }

    // Interior m = 1 .. len - 2: both (i-1, j) and (i, j-1) live contiguously
    // on diagonal d - 1, so this loop is data-parallel and branchless.
    if (len >= 3) {
      const Real *AMREX_RESTRICT rCur = rIn + off;
      const Real *AMREX_RESTRICT e0Cur = e0 + off;
      const Real *AMREX_RESTRICT e1Cur = e1 + off;
      const Real *AMREX_RESTRICT dICur = dI + off;
      const Real *zPrevI = zOut + (dg.prevBase - 1);
      const Real *zPrevJ = zOut + dg.prevBase;
      Real *zCur = zOut + off;
#pragma omp simd
      for (int m = 1; m < len - 1; ++m) {
        zCur[m] =
            (rCur[m] - e0Cur[m] * zPrevI[m] - e1Cur[m] * zPrevJ[m]) * dICur[m];
      }
    }

    // Endpoint m = len - 1
    if (len > 1) {
      const int m = len - 1;
      const int i = dg.iStart + m;
      const int j = d - i;
      const int k = off + m;
      Real val = rIn[k];
      if (i > 0)
        val -= e0[k] * zOut[dg.prevBase - 1 + m];
      if (j > 0)
        val -= e1[k] * zOut[dg.prevBase + m];
      zOut[k] = val * dI[k];
    }
  }

  // Backward wavefront solve: U z = y using pre-scaled fUs, f1Us (d = nDiag - 2
  // down to 0).
  for (int d = nDiag - 2; d >= 0; --d) {
    const WavefrontDiag &dg = diags[d];
    const int len = dg.len;
    const int off = dg.off;

    // Endpoint m = 0
    {
      const int i = dg.iStart;
      const int j = d - i;
      Real sub = 0.0;
      if (i < nTh - 1)
        sub += fu[off] * zOut[dg.nextBase + 1];
      if (j < nPh - 1)
        sub += f1u[off] * zOut[dg.nextBase];
      zOut[off] -= sub;
    }

    // Interior m = 1 .. len - 2: both (i+1, j) and (i, j+1) live contiguously
    // on diagonal d + 1.
    if (len >= 3) {
      const Real *AMREX_RESTRICT fuCur = fu + off;
      const Real *AMREX_RESTRICT f1uCur = f1u + off;
      const Real *zNextI = zOut + (dg.nextBase + 1);
      const Real *zNextJ = zOut + dg.nextBase;
      Real *zCur = zOut + off;
#pragma omp simd
      for (int m = 1; m < len - 1; ++m) {
        zCur[m] -= fuCur[m] * zNextI[m] + f1uCur[m] * zNextJ[m];
      }
    }

    // Endpoint m = len - 1
    if (len > 1) {
      const int m = len - 1;
      const int i = dg.iStart + m;
      const int j = d - i;
      const int k = off + m;
      Real sub = 0.0;
      if (i < nTh - 1)
        sub += fu[k] * zOut[dg.nextBase + 1 + m];
      if (j < nPh - 1)
        sub += f1u[k] * zOut[dg.nextBase + m];
      zOut[k] -= sub;
    }
  }
}

//==========================================================
void HemisphereWorkspace::matvec(const std::vector<Real> &in,
                                 std::vector<Real> &out) const {
  if (nTh <= 0 || nPh <= 0)
    return;

  const int nDiag = static_cast<int>(diags.size());
  const Real *AMREX_RESTRICT xIn = in.data();
  Real *AMREX_RESTRICT yOut = out.data();
  const Real *AMREX_RESTRICT a = cA.data();
  const Real *AMREX_RESTRICT b = cB.data();
  const Real *AMREX_RESTRICT c = cC.data();
  const Real *AMREX_RESTRICT dC = cD.data();
  const Real *AMREX_RESTRICT eC = cE.data();

  // Compute pole average from the row adjacent to the pole.
  Real poleAvg = 0.0;
  for (int j = 0; j < nPh; ++j) {
    poleAvg += xIn[poleWfIdx[j]];
  }
  poleAvg /= static_cast<Real>(nPh);

  const Real neuFactor = isNeumann ? Real(1.0) : Real(0.0);

  for (int d = 0; d < nDiag; ++d) {
    const WavefrontDiag &dg = diags[d];
    const int len = dg.len;
    const int off = dg.off;

    auto eval_endpoint = [&](int m) noexcept {
      const int i = dg.iStart + m;
      const int j = d - i;
      const int k = off + m;
      const Real xIm1 = (i > 0) ? xIn[dg.prevBase - 1 + m]
                                : (isNorth ? poleAvg : neuFactor * xIn[k]);
      const Real xIp1 = (i < nTh - 1)
                            ? xIn[dg.nextBase + 1 + m]
                            : (isNorth ? neuFactor * xIn[k] : poleAvg);
      const Real xJm1 = xIn[wf_index(i, (j == 0) ? (nPh - 1) : (j - 1))];
      const Real xJp1 = xIn[wf_index(i, (j == nPh - 1) ? 0 : (j + 1))];
      yOut[k] = a[k] * xIn[k] + b[k] * xIm1 + c[k] * xIp1 + dC[k] * xJm1 +
                eC[k] * xJp1;
    };

    eval_endpoint(0);

    // Interior m = 1 .. len - 2: all 4 neighbors are strictly in-bounds and
    // unit-stride contiguous on diagonals d - 1 and d + 1.
    if (len >= 3) {
      const Real *xCur = xIn + off;
      const Real *xPrevI = xIn + (dg.prevBase - 1);
      const Real *xNextI = xIn + (dg.nextBase + 1);
      const Real *xPrevJ = xIn + dg.prevBase;
      const Real *xNextJ = xIn + dg.nextBase;
      const Real *AMREX_RESTRICT aCur = a + off;
      const Real *AMREX_RESTRICT bCur = b + off;
      const Real *AMREX_RESTRICT cCur = c + off;
      const Real *AMREX_RESTRICT dCur = dC + off;
      const Real *AMREX_RESTRICT eCur = eC + off;
      Real *AMREX_RESTRICT yCur = yOut + off;

#pragma omp simd
      for (int m = 1; m < len - 1; ++m) {
        yCur[m] = aCur[m] * xCur[m] + bCur[m] * xPrevI[m] +
                  cCur[m] * xNextI[m] + dCur[m] * xPrevJ[m] +
                  eCur[m] * xNextJ[m];
      }
    }

    if (len > 1) {
      eval_endpoint(len - 1);
    }
  }

  // For pure Neumann low-lat BC, pin one gauge point to remove the nullspace.
  if (isNeumann) {
    yOut[kPin] = xIn[kPin];
  }
}

//==========================================================
void Ionosphere::build_dipole_frame(const Real dipoleMoment[3], Real rot[3][3],
                                    Real &dipoleSign) {
  const Real mNorm = std::sqrt(dipoleMoment[0] * dipoleMoment[0] +
                               dipoleMoment[1] * dipoleMoment[1] +
                               dipoleMoment[2] * dipoleMoment[2]);
  if (mNorm <= Real(1.0e-30)) {
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        rot[a][b] = (a == b) ? Real(1.0) : Real(0.0);
    dipoleSign = -1.0;
    return;
  }

  Real ax[3] = {
    dipoleMoment[0] / mNorm,
    dipoleMoment[1] / mNorm,
    dipoleMoment[2] / mNorm,
  };

  dipoleSign = (ax[2] < Real(0.0)) ? Real(-1.0) : Real(1.0);
  Real ez[3] = { dipoleSign * ax[0], dipoleSign * ax[1], dipoleSign * ax[2] };

  Real ref[3] = { 1.0, 0.0, 0.0 };
  if (std::abs(ez[0]) > Real(0.99)) {
    ref[0] = 0.0;
    ref[1] = 1.0;
    ref[2] = 0.0;
  }
  const Real dotRef = ref[0] * ez[0] + ref[1] * ez[1] + ref[2] * ez[2];
  Real ex[3] = {
    ref[0] - dotRef * ez[0],
    ref[1] - dotRef * ez[1],
    ref[2] - dotRef * ez[2],
  };
  const Real exNorm = std::sqrt(ex[0] * ex[0] + ex[1] * ex[1] + ex[2] * ex[2]);
  ex[0] /= exNorm;
  ex[1] /= exNorm;
  ex[2] /= exNorm;

  Real ey[3] = {
    ez[1] * ex[2] - ez[2] * ex[1],
    ez[2] * ex[0] - ez[0] * ex[2],
    ez[0] * ex[1] - ez[1] * ex[0],
  };

  for (int b = 0; b < 3; ++b) {
    rot[0][b] = ex[b];
    rot[1][b] = ey[b];
    rot[2][b] = ez[b];
  }
}

//==========================================================
void Ionosphere::init(const Real dipoleMomentCode[3],
                      const Real bodyCenterCode[3], Real rPlanetSI, Real si2NoL,
                      Real si2NoE, Real si2NoB, Real bEqSI) {
  (void)si2NoB;
  if (params_.nThetaHem < 5)
    Abort("Ionosphere::init: nThetaHem must be >= 5.");
  if (params_.nPhi < 5)
    Abort("Ionosphere::init: nPhi must be >= 5.");
  if (rPlanetSI <= 0.0)
    Abort("Ionosphere::init: rPlanetSI must be positive.");

  rPlanetSI_ = rPlanetSI;
  bEqSI_ = bEqSI;

  for (int d = 0; d < 3; ++d)
    frame_.center[d] = bodyCenterCode[d];

  build_dipole_frame(dipoleMomentCode, frame_.rot, frame_.dipoleSign);

  frame_.rIonoSI = params_.rIonoRatio * rPlanetSI;
  frame_.rIonoCode = frame_.rIonoSI * si2NoL;
  frame_.si2NoE = si2NoE;
  frame_.no2SiL = Real(1.0) / si2NoL;

  nTheta_ = 2 * params_.nThetaHem - 1;
  nPhi_ = params_.nPhi;
  dTheta_ = dPI / static_cast<Real>(nTheta_ - 1);
  dPhi_ = Real(2.0) * dPI / static_cast<Real>(nPhi_ - 1);
  invDTheta_ = Real(1.0) / dTheta_;
  invDPhi_ = Real(1.0) / dPhi_;
  latBoundaryRad_ = params_.latBoundaryDeg * cDegToRad;

  theta_.resize(nTheta_);
  sinTheta_.resize(nTheta_);
  cosTheta_.resize(nTheta_);
  dPhiCorot_.resize(nTheta_);
  const Real corotCoef = -params_.omegaPlanetSI * bEqSI_ *
                         (rPlanetSI_ * rPlanetSI_ * rPlanetSI_) /
                         frame_.rIonoSI;
  for (int i = 0; i < nTheta_; ++i) {
    theta_[i] = i * dTheta_;
    sinTheta_[i] = std::sin(theta_[i]);
    cosTheta_[i] = std::cos(theta_[i]);
    dPhiCorot_[i] = corotCoef * sinTheta_[i] * sinTheta_[i];
  }

  phi_.resize(nPhi_);
  sinPhi_.resize(nPhi_);
  cosPhi_.resize(nPhi_);
  for (int j = 0; j < nPhi_; ++j) {
    phi_[j] = j * dPhi_;
    sinPhi_[j] = std::sin(phi_[j]);
    cosPhi_[j] = std::cos(phi_[j]);
  }

  const size_t nTot = static_cast<size_t>(nTheta_) * nPhi_;
  jr_.assign(nTot, 0.0);
  jrTmp_.assign(nTot, 0.0);
  sigmaPBase_.assign(nTot, params_.sigmaPedersen);
  sigmaHBase_.assign(nTot, params_.sigmaHall);
  sigmaP_.assign(nTot, params_.sigmaPedersen);
  sigmaH_.assign(nTot, params_.sigmaHall);
  sigmaThTh_.assign(nTot, 0.0);
  sigmaThPh_.assign(nTot, 0.0);
  sigmaPhPh_.assign(nTot, 0.0);
  dSigmaThTh_dTh_.assign(nTot, 0.0);
  dSigmaThPh_dTh_.assign(nTot, 0.0);
  dSigmaThPh_dPh_.assign(nTot, 0.0);
  dSigmaPhPh_dPh_.assign(nTot, 0.0);

  phiPot_.assign(nTot, 0.0);
  phiSolved_.assign(nTot, 0.0);
  dPhiDTheta_.assign(nTot, 0.0);
  dPhiDPhi_.assign(nTot, 0.0);
  bcastBuf_.assign(2 * nTot + 3, 0.0);
  hasPrevSolution_ = false;

  // Determine active colatitude index ranges for North and South hemispheres
  // and pre-allocate persistent operator/BiCGSTAB workspaces.
  constexpr Real epsBnd = Real(1.0e-12);
  const int iEq = params_.nThetaHem - 1;
  int iMaxN = 0;
  for (int i = 1; i < iEq; ++i) {
    if (Real(0.5) * dPI - theta_[i] > latBoundaryRad_ + epsBnd)
      iMaxN = i;
  }
  int iMinS = nTheta_ - 1;
  for (int i = nTheta_ - 2; i > iEq; --i) {
    if (theta_[i] - Real(0.5) * dPI > latBoundaryRad_ + epsBnd)
      iMinS = i;
  }

  const bool isNeumann = (params_.lowLatBC == "neumann");
  const int nPhiUsed = nPhi_ - 1;
  hemWork_[0].allocate(1, iMaxN, nPhiUsed, true, isNeumann);
  hemWork_[1].allocate(iMinS, nTheta_ - 2, nPhiUsed, false, isNeumann);

  isStaticConductance_ = (params_.conductanceModel != "auroral");
  matrixBuilt_ = false;

  compute_base_conductances();
  compute_conductances();
  compute_conductance_tensor();
  build_hemisphere_matrix(hemWork_[0]);
  build_hemisphere_matrix(hemWork_[1]);
  matrixBuilt_ = true;

  sync_device_gradients();
  ++generation_;
  initialized_ = true;
}

//==========================================================
void Ionosphere::compute_base_conductances() {
  const bool useSolar = (params_.conductanceModel == "solar_euv" ||
                         params_.conductanceModel == "auroral");

  const Real f107Term = std::pow(std::max(params_.f107, Real(1.0)), Real(0.49));
  const Real f107Hall = std::pow(std::max(params_.f107, Real(1.0)), Real(0.53));

  for (int i = 0; i < nTheta_; ++i) {
    const Real sn = sinTheta_[i];
    const Real zSmg = cosTheta_[i];
    for (int j = 0; j < nPhi_; ++j) {
      const int k = idx2d(i, j, nPhi_);
      Real sp = params_.sigmaPedersen;
      Real sh = params_.sigmaHall;

      if (useSolar) {
        const Real xSmg = sn * cosPhi_[j];
        const Real ySmg = sn * sinPhi_[j];
        const Real cosChi = frame_.rot[0][0] * xSmg + frame_.rot[1][0] * ySmg +
                            frame_.rot[2][0] * zSmg;
        if (cosChi > Real(0.0)) {
          const Real sqrtCos = std::sqrt(cosChi);
          const Real spEuv =
              f107Term * (Real(0.34) * cosChi + Real(0.93) * sqrtCos);
          const Real shEuv =
              f107Hall * (Real(0.81) * cosChi + Real(0.54) * sqrtCos);
          sp = std::sqrt(sp * sp + spEuv * spEuv);
          sh = std::sqrt(sh * sh + shEuv * shEuv);
        }
      }

      sigmaPBase_[k] = std::max(sp, Real(0.1));
      sigmaHBase_[k] = sh;
    }
  }
}

//==========================================================
void Ionosphere::compute_conductances() {
  if (isStaticConductance_) {
    sigmaP_ = sigmaPBase_;
    sigmaH_ = sigmaHBase_;
    return;
  }

  // Auroral conductance: start from precomputed solar-EUV/background base and
  // add upward-FAC enhancement only where jr > 0.
  const size_t nTot = sigmaPBase_.size();
  const Real invFacScale = (params_.auroralFacScale > Real(0.0))
                               ? (Real(1.0) / params_.auroralFacScale)
                               : Real(0.0);
  const Real spMax = params_.auroralSigmaPMax;
  const Real hallRatio = params_.auroralHallRatio;

  for (size_t k = 0; k < nTot; ++k) {
    Real sp = sigmaPBase_[k];
    Real sh = sigmaHBase_[k];
    const Real jUp = jr_[k];
    if (jUp > Real(0.0) && invFacScale > Real(0.0)) {
      const Real spAur = spMax * std::tanh(std::sqrt(jUp * invFacScale));
      const Real shAur = hallRatio * spAur;
      sp = std::sqrt(sp * sp + spAur * spAur);
      sh = std::sqrt(sh * sh + shAur * shAur);
    }
    sigmaP_[k] = std::max(sp, Real(0.1));
    sigmaH_[k] = sh;
  }
}

//==========================================================
void Ionosphere::compute_conductance_tensor() {
  const Real s0 = std::max(params_.sigmaParallel, Real(10.0));
  const Real hallSign = -frame_.dipoleSign;

  for (int i = 0; i < nTheta_; ++i) {
    const Real sn = sinTheta_[i];
    const Real cs = cosTheta_[i];
    const Real sn2 = sn * sn;
    const Real cs2 = cs * cs;
    const Real cs3 = Real(1.0) + Real(3.0) * cs2;
    const Real cs4 = std::sqrt(cs3);

    for (int j = 0; j < nPhi_; ++j) {
      const int k = idx2d(i, j, nPhi_);
      const Real sp = sigmaP_[k];
      const Real sh = sigmaH_[k];
      const Real invC = Real(1.0) / (Real(4.0) * s0 * cs2 + sp * sn2);

      sigmaThTh_[k] = s0 * sp * cs3 * invC;
      sigmaThPh_[k] = hallSign * Real(2.0) * s0 * sh * cs * cs4 * invC;
      sigmaPhPh_[k] = sp + (sh * sh * sn2) * invC;
    }
  }

  const Real inv2dTh = Real(0.5) * invDTheta_;
  const Real inv2dPh = Real(0.5) * invDPhi_;
  const int nPhiUsed = nPhi_ - 1;

  const Real *AMREX_RESTRICT sTT = sigmaThTh_.data();
  const Real *AMREX_RESTRICT sTP = sigmaThPh_.data();
  const Real *AMREX_RESTRICT sPP = sigmaPhPh_.data();
  Real *AMREX_RESTRICT dTTdTh = dSigmaThTh_dTh_.data();
  Real *AMREX_RESTRICT dTPdTh = dSigmaThPh_dTh_.data();
  Real *AMREX_RESTRICT dTPdPh = dSigmaThPh_dPh_.data();
  Real *AMREX_RESTRICT dPPdPh = dSigmaPhPh_dPh_.data();

  // Theta derivatives: peel poles (i = 0, i = nTheta_ - 1) so interior rows
  // vectorize with unit-stride access.
  std::fill(dTTdTh, dTTdTh + nPhi_, Real(0.0));
  std::fill(dTPdTh, dTPdTh + nPhi_, Real(0.0));
  if (nTheta_ > 1) {
    const int lastOff = (nTheta_ - 1) * nPhi_;
    std::fill(dTTdTh + lastOff, dTTdTh + lastOff + nPhi_, Real(0.0));
    std::fill(dTPdTh + lastOff, dTPdTh + lastOff + nPhi_, Real(0.0));
  }
  for (int i = 1; i < nTheta_ - 1; ++i) {
    const Real *AMREX_RESTRICT ttM = sTT + (i - 1) * nPhi_;
    const Real *AMREX_RESTRICT ttP = sTT + (i + 1) * nPhi_;
    const Real *AMREX_RESTRICT tpM = sTP + (i - 1) * nPhi_;
    const Real *AMREX_RESTRICT tpP = sTP + (i + 1) * nPhi_;
    Real *AMREX_RESTRICT dstTT = dTTdTh + i * nPhi_;
    Real *AMREX_RESTRICT dstTP = dTPdTh + i * nPhi_;
#pragma omp simd
    for (int j = 0; j < nPhi_; ++j) {
      dstTT[j] = (ttP[j] - ttM[j]) * inv2dTh;
      dstTP[j] = (tpP[j] - tpM[j]) * inv2dTh;
    }
  }

  // Periodic phi derivatives: peel j = 0 and j = nPhiUsed - 1 (plus wrap).
  for (int i = 0; i < nTheta_; ++i) {
    const int rowOff = i * nPhi_;
    const Real *AMREX_RESTRICT tpRow = sTP + rowOff;
    const Real *AMREX_RESTRICT ppRow = sPP + rowOff;
    Real *AMREX_RESTRICT dstTP = dTPdPh + rowOff;
    Real *AMREX_RESTRICT dstPP = dPPdPh + rowOff;

    dstTP[0] = (tpRow[1] - tpRow[nPhiUsed - 1]) * inv2dPh;
    dstPP[0] = (ppRow[1] - ppRow[nPhiUsed - 1]) * inv2dPh;

#pragma omp simd
    for (int j = 1; j < nPhiUsed - 1; ++j) {
      dstTP[j] = (tpRow[j + 1] - tpRow[j - 1]) * inv2dPh;
      dstPP[j] = (ppRow[j + 1] - ppRow[j - 1]) * inv2dPh;
    }

    if (nPhiUsed > 1) {
      const int jLast = nPhiUsed - 1;
      dstTP[jLast] = (tpRow[0] - tpRow[jLast - 1]) * inv2dPh;
      dstPP[jLast] = (ppRow[0] - ppRow[jLast - 1]) * inv2dPh;
    }
    dstTP[nPhi_ - 1] = dstTP[0];
    dstPP[nPhi_ - 1] = dstPP[0];
  }
}

//==========================================================
void Ionosphere::build_hemisphere_matrix(HemisphereWorkspace &hw) {
  if (hw.nTh <= 0 || hw.nPh <= 0)
    return;

  const Real invDTh2 = invDTheta_ * invDTheta_;
  const Real invDPh2 = invDPhi_ * invDPhi_;
  const Real inv2DTh = Real(0.5) * invDTheta_;
  const Real inv2DPh = Real(0.5) * invDPhi_;

  for (int j = 0; j < hw.nPh; ++j) {
    for (int iLoc = 0; iLoc < hw.nTh; ++iLoc) {
      const int i = hw.iMin + iLoc;
      const int kGrid = idx2d(i, j, nPhi_);
      const int k = hw.wf_index(iLoc, j);

      const Real sn = sinTheta_[i];
      const Real cs = cosTheta_[i];
      const Real sn2 = sn * sn;

      const Real termTh2 = sigmaThTh_[kGrid] * sn2 * invDTh2;
      const Real termTh1 =
          (sigmaThTh_[kGrid] * sn * cs + dSigmaThTh_dTh_[kGrid] * sn2 -
           dSigmaThPh_dPh_[kGrid] * sn) *
          inv2DTh;
      const Real termPh2 = sigmaPhPh_[kGrid] * invDPh2;
      const Real termPh1 =
          (dSigmaThPh_dTh_[kGrid] * sn + dSigmaPhPh_dPh_[kGrid]) * inv2DPh;

      hw.cA[k] = -Real(2.0) * (termTh2 + termPh2);
      hw.cB[k] = termTh2 - termTh1;
      hw.cC[k] = termTh2 + termTh1;
      hw.cD[k] = termPh2 - termPh1;
      hw.cE[k] = termPh2 + termPh1;
    }
  }

  if (hw.isNeumann) {
    const int kPin = hw.kPin;
    hw.cA[kPin] = 1.0;
    hw.cB[kPin] = hw.cC[kPin] = hw.cD[kPin] = hw.cE[kPin] = 0.0;
  }

  hw.build_ilu0();
}

//==========================================================
void Ionosphere::smooth_jr() {
  if (params_.nSmoothFac <= 0)
    return;

  const int nPhiUsed = nPhi_ - 1;
  Real *AMREX_RESTRICT jrPtr = jr_.data();
  Real *AMREX_RESTRICT tmpPtr = jrTmp_.data();

  for (int pass = 0; pass < params_.nSmoothFac; ++pass) {
    for (int i = 1; i < nTheta_ - 1; ++i) {
      const int row = i * nPhi_;
      const int rowM = (i - 1) * nPhi_;
      const int rowP = (i + 1) * nPhi_;

      // j = 0 (wraps to jPrev = nPhiUsed - 1)
      tmpPtr[row] = Real(0.5) * jrPtr[row] +
                    Real(0.125) * (jrPtr[rowM] + jrPtr[rowP] +
                                   jrPtr[row + nPhiUsed - 1] + jrPtr[row + 1]);

#pragma omp simd
      for (int j = 1; j < nPhiUsed - 1; ++j) {
        tmpPtr[row + j] =
            Real(0.5) * jrPtr[row + j] +
            Real(0.125) * (jrPtr[rowM + j] + jrPtr[rowP + j] +
                           jrPtr[row + j - 1] + jrPtr[row + j + 1]);
      }

      // j = nPhiUsed - 1 (wraps to jNext = 0)
      if (nPhiUsed > 1) {
        const int jLast = nPhiUsed - 1;
        tmpPtr[row + jLast] =
            Real(0.5) * jrPtr[row + jLast] +
            Real(0.125) * (jrPtr[rowM + jLast] + jrPtr[rowP + jLast] +
                           jrPtr[row + jLast - 1] + jrPtr[row]);
      }

      tmpPtr[row + nPhi_ - 1] = tmpPtr[row];
    }

    // Poles: average of adjacent ring.
    Real northAvg = 0.0, southAvg = 0.0;
    const int rowNorthAdj = nPhi_;
    const int rowSouthAdj = (nTheta_ - 2) * nPhi_;
    for (int j = 0; j < nPhiUsed; ++j) {
      northAvg += tmpPtr[rowNorthAdj + j];
      southAvg += tmpPtr[rowSouthAdj + j];
    }
    const Real invN = Real(1.0) / static_cast<Real>(nPhiUsed);
    northAvg *= invN;
    southAvg *= invN;
    const int rowSouthPole = (nTheta_ - 1) * nPhi_;
    for (int j = 0; j < nPhi_; ++j) {
      tmpPtr[j] = northAvg;
      tmpPtr[rowSouthPole + j] = southAvg;
    }

    jr_.swap(jrTmp_);
    jrPtr = jr_.data();
    tmpPtr = jrTmp_.data();
  }
}

//==========================================================
bool Ionosphere::solve_hemisphere(HemisphereWorkspace &hw, int &nIterOut,
                                  Real &residOut) {
  if (hw.nTh <= 0 || hw.nPh <= 0) {
    nIterOut = 0;
    residOut = 0.0;
    return true;
  }

  const int nPhiUsed = nPhi_ - 1;
  const int iEq = params_.nThetaHem - 1;
  const Real rIono2 = frame_.rIonoSI * frame_.rIonoSI;

  Real jrShift = 0.0;
  if (hw.isNeumann) {
    Real num = 0.0, den = 0.0;
    const Real *AMREX_RESTRICT jrPtr = jr_.data();
    for (int i = hw.iMin; i <= hw.iMax; ++i) {
      const Real *AMREX_RESTRICT row = jrPtr + i * nPhi_;
      Real rowSum = 0.0;
#pragma omp simd reduction(+ : rowSum)
      for (int j = 0; j < nPhiUsed; ++j) {
        rowSum += row[j];
      }
      const Real w = sinTheta_[i];
      num += rowSum * w;
      den += w * static_cast<Real>(nPhiUsed);
    }
    if (den > Real(0.0))
      jrShift = num / den;
  }

  for (int j = 0; j < hw.nPh; ++j) {
    for (int iLoc = 0; iLoc < hw.nTh; ++iLoc) {
      const int i = hw.iMin + iLoc;
      const int kGrid = idx2d(i, j, nPhi_);
      const int k = hw.wf_index(iLoc, j);
      const Real sn = sinTheta_[i];
      hw.rhs[k] = (jr_[kGrid] - jrShift) * rIono2 * sn * sn;
      hw.x[k] = hasPrevSolution_ ? phiSolved_[kGrid] : Real(0.0);
    }
  }

  if (hw.isNeumann) {
    hw.rhs[hw.kPin] = 0.0;
  }

  const bool converged = solve_bicgstab(
      hw, params_.solverTol, params_.solverMaxIter, nIterOut, residOut);

  // Scatter solution into this hemisphere's disjoint rows of phiSolved_.
  const int iStartClear = hw.isNorth ? 0 : (iEq + 1);
  const int iEndClear = hw.isNorth ? iEq : (nTheta_ - 1);
  for (int i = iStartClear; i <= iEndClear; ++i) {
    std::fill(phiSolved_.begin() + i * nPhi_,
              phiSolved_.begin() + (i + 1) * nPhi_, Real(0.0));
  }

  const int nSol = static_cast<int>(hw.x.size());
  Real phiMin = hw.x[0], phiMax = hw.x[0];
  const Real *AMREX_RESTRICT xSol = hw.x.data();
#pragma omp simd reduction(min : phiMin) reduction(max : phiMax)
  for (int k = 0; k < nSol; ++k) {
    const Real val = xSol[k];
    phiMin = std::min(phiMin, val);
    phiMax = std::max(phiMax, val);
  }

  for (int iLoc = 0; iLoc < hw.nTh; ++iLoc) {
    const int i = hw.iMin + iLoc;
    Real *AMREX_RESTRICT dstRow = phiSolved_.data() + i * nPhi_;
    for (int j = 0; j < hw.nPh; ++j) {
      dstRow[j] = hw.x[hw.wf_index(iLoc, j)];
    }
  }

  // Pole average.
  const int iPole = hw.isNorth ? 0 : (nTheta_ - 1);
  const int iAdj = hw.isNorth ? 1 : (nTheta_ - 2);
  Real poleVal = 0.0;
  for (int j = 0; j < nPhiUsed; ++j)
    poleVal += phiSolved_[idx2d(iAdj, j, nPhi_)];
  poleVal /= static_cast<Real>(nPhiUsed);
  for (int j = 0; j < nPhi_; ++j)
    phiSolved_[idx2d(iPole, j, nPhi_)] = poleVal;

  // Periodic wrap in longitude.
  for (int i = iStartClear; i <= iEndClear; ++i) {
    phiSolved_[idx2d(i, nPhi_ - 1, nPhi_)] = phiSolved_[idx2d(i, 0, nPhi_)];
  }

  if (hw.isNorth) {
    cpcpNorthKV_ = (phiMax - phiMin) * Real(1.0e-3);
  } else {
    cpcpSouthKV_ = (phiMax - phiMin) * Real(1.0e-3);
  }

  return converged;
}

//==========================================================
void Ionosphere::add_corotation_potential() {
  if (!params_.useCorotation)
    return;

  for (int i = 0; i < nTheta_; ++i) {
    const Real dPhiCorot = dPhiCorot_[i];
    Real *AMREX_RESTRICT row = phiPot_.data() + i * nPhi_;
#pragma omp simd
    for (int j = 0; j < nPhi_; ++j) {
      row[j] += dPhiCorot;
    }
  }
}

//==========================================================
void Ionosphere::compute_gradients() {
  const int nPhiUsed = nPhi_ - 1;
  const Real inv2dTh = Real(0.5) * invDTheta_;
  const Real inv2dPh = Real(0.5) * invDPhi_;
  const Real *AMREX_RESTRICT pPot = phiPot_.data();
  Real *AMREX_RESTRICT dTh = dPhiDTheta_.data();
  Real *AMREX_RESTRICT dPh = dPhiDPhi_.data();

  // Interior central differences in theta.
  for (int i = 1; i < nTheta_ - 1; ++i) {
    const Real *AMREX_RESTRICT rowM = pPot + (i - 1) * nPhi_;
    const Real *AMREX_RESTRICT rowP = pPot + (i + 1) * nPhi_;
    Real *AMREX_RESTRICT dst = dTh + i * nPhi_;
#pragma omp simd
    for (int j = 0; j < nPhi_; ++j) {
      dst[j] = (rowP[j] - rowM[j]) * inv2dTh;
    }
  }

  // Second-order one-sided differences at the poles.
  {
    const Real *AMREX_RESTRICT r0 = pPot;
    const Real *AMREX_RESTRICT r1 = pPot + nPhi_;
    const Real *AMREX_RESTRICT r2 = pPot + 2 * nPhi_;
    Real *AMREX_RESTRICT dstN = dTh;
    const int iLast = nTheta_ - 1;
    const Real *AMREX_RESTRICT rs0 = pPot + iLast * nPhi_;
    const Real *AMREX_RESTRICT rs1 = pPot + (iLast - 1) * nPhi_;
    const Real *AMREX_RESTRICT rs2 = pPot + (iLast - 2) * nPhi_;
    Real *AMREX_RESTRICT dstS = dTh + iLast * nPhi_;
#pragma omp simd
    for (int j = 0; j < nPhi_; ++j) {
      dstN[j] = (Real(4.0) * r1[j] - Real(3.0) * r0[j] - r2[j]) * inv2dTh;
      dstS[j] = (Real(3.0) * rs0[j] - Real(4.0) * rs1[j] + rs2[j]) * inv2dTh;
    }
  }

  // Periodic central differences in phi.
  for (int i = 0; i < nTheta_; ++i) {
    const Real *AMREX_RESTRICT row = pPot + i * nPhi_;
    Real *AMREX_RESTRICT dst = dPh + i * nPhi_;
    dst[0] = (row[1] - row[nPhiUsed - 1]) * inv2dPh;
#pragma omp simd
    for (int j = 1; j < nPhiUsed; ++j) {
      dst[j] = (row[j + 1] - row[j - 1]) * inv2dPh;
    }
    if (nPhiUsed > 1) {
      dst[nPhiUsed - 1] = (row[0] - row[nPhiUsed - 2]) * inv2dPh;
    }
    dst[nPhi_ - 1] = dst[0];
  }
}

//==========================================================
void Ionosphere::sync_device_gradients() {
  const size_t nTot = dPhiDTheta_.size();
  dPhiDThetaDev_.resize(nTot);
  dPhiDPhiDev_.resize(nTot);
  if (nTot > 0) {
    Gpu::copyAsync(Gpu::hostToDevice, dPhiDTheta_.begin(), dPhiDTheta_.end(),
                   dPhiDThetaDev_.begin());
    Gpu::copyAsync(Gpu::hostToDevice, dPhiDPhi_.begin(), dPhiDPhi_.end(),
                   dPhiDPhiDev_.begin());
    Gpu::streamSynchronize();
  }
}

//==========================================================
bool Ionosphere::solve(const std::vector<Real> &jrOutwardSI, bool doBroadcast) {
  if (!initialized_)
    Abort("Ionosphere::solve called before init().");

  const size_t nTot = static_cast<size_t>(nTheta_) * nPhi_;
  if (jrOutwardSI.size() != nTot)
    Abort("Ionosphere::solve: jrOutwardSI size mismatch.");

  int convFlag = 1;

  if (!doBroadcast || ParallelDescriptor::IOProcessor()) {
    jr_ = jrOutwardSI;
    smooth_jr();

    if (!isStaticConductance_ || !matrixBuilt_) {
      compute_conductances();
      compute_conductance_tensor();
      build_hemisphere_matrix(hemWork_[0]);
      build_hemisphere_matrix(hemWork_[1]);
      matrixBuilt_ = true;
    }

    bool convN = true;
    bool convS = true;
    parallel_sections(
        [&]() {
          convN =
              solve_hemisphere(hemWork_[0], lastIterNorth_, lastResidNorth_);
        },
        [&]() {
          convS =
              solve_hemisphere(hemWork_[1], lastIterSouth_, lastResidSouth_);
        });
    convFlag = (convN && convS) ? 1 : 0;

    // Temporal blending with previous solution.
    const Real alpha = std::clamp(params_.relaxAlpha, Real(0.0), Real(1.0));
    if (!hasPrevSolution_ || alpha >= Real(1.0)) {
      phiPot_ = phiSolved_;
    } else {
      const Real oneMinusAlpha = Real(1.0) - alpha;
      Real *AMREX_RESTRICT pPot = phiPot_.data();
      const Real *AMREX_RESTRICT pSol = phiSolved_.data();
#pragma omp simd
      for (size_t k = 0; k < nTot; ++k) {
        pPot[k] = oneMinusAlpha * pPot[k] + alpha * pSol[k];
      }
    }
    hasPrevSolution_ = true;

    add_corotation_potential();
    compute_gradients();
  }

  if (doBroadcast && ParallelDescriptor::NProcs() > 1) {
    const int root = ParallelDescriptor::IOProcessorNumber();
    if (ParallelDescriptor::IOProcessor()) {
      std::copy(dPhiDTheta_.begin(), dPhiDTheta_.end(), bcastBuf_.begin());
      std::copy(dPhiDPhi_.begin(), dPhiDPhi_.end(), bcastBuf_.begin() + nTot);
      bcastBuf_[2 * nTot + 0] = cpcpNorthKV_;
      bcastBuf_[2 * nTot + 1] = cpcpSouthKV_;
      bcastBuf_[2 * nTot + 2] = static_cast<Real>(convFlag);
    }
    ParallelDescriptor::Bcast(bcastBuf_.data(), bcastBuf_.size(), root);
    if (!ParallelDescriptor::IOProcessor()) {
      std::copy(bcastBuf_.begin(), bcastBuf_.begin() + nTot,
                dPhiDTheta_.begin());
      std::copy(bcastBuf_.begin() + nTot, bcastBuf_.begin() + 2 * nTot,
                dPhiDPhi_.begin());
      cpcpNorthKV_ = bcastBuf_[2 * nTot + 0];
      cpcpSouthKV_ = bcastBuf_[2 * nTot + 1];
      convFlag = (bcastBuf_[2 * nTot + 2] > Real(0.5)) ? 1 : 0;
    }
    hasPrevSolution_ = true;
  }

  sync_device_gradients();
  ++generation_;
  return convFlag != 0;
}

//==========================================================
void Ionosphere::set_potential(const std::vector<Real> &phiIn) {
  const size_t nTot = static_cast<size_t>(nTheta_) * nPhi_;
  if (phiIn.size() != nTot)
    Abort("Ionosphere::set_potential: size mismatch.");
  phiPot_ = phiIn;
  phiSolved_ = phiIn;
  hasPrevSolution_ = true;
  compute_gradients();
  sync_device_gradients();
  ++generation_;
}

//==========================================================
IonosphereDeviceView Ionosphere::device_view() const {
  IonosphereDeviceView v;
  v.frame = frame_;
  v.nTheta = nTheta_;
  v.nPhi = nPhi_;
  v.invDTheta = invDTheta_;
  v.invDPhi = invDPhi_;
  v.eScale = -frame_.si2NoE / frame_.no2SiL;
  v.latBoundaryRad = latBoundaryRad_;
  v.dPhiDTheta = dPhiDThetaDev_.data();
  v.dPhiDPhi = dPhiDPhiDev_.data();
  return v;
}

//==========================================================
void Ionosphere::write_diagnostics(const std::string &fileName, Real timeSI,
                                   int cycle) const {
  if (!ParallelDescriptor::IOProcessor())
    return;

  std::ofstream out(fileName);
  if (!out.is_open()) {
    Print() << "WARNING: Ionosphere::write_diagnostics could not open "
            << fileName << "\n";
    return;
  }

  out << "# FLEKS Ionosphere Diagnostics\n";
  out << "# cycle = " << cycle << "  time_SI = " << std::scientific
      << std::setprecision(6) << timeSI << "  cpcpN_kV = " << cpcpNorthKV_
      << "  cpcpS_kV = " << cpcpSouthKV_ << "\n";
  out << "# nTheta = " << nTheta_ << "  nPhi = " << nPhi_
      << "  rIono_m = " << frame_.rIonoSI << "\n";
  out << "# theta_deg phi_deg Jr_uAm2 SigmaP_S SigmaH_S Phi_kV Eth_mVm "
         "Eph_mVm\n";

  const Real invRI = Real(1.0) / frame_.rIonoSI;
  for (int i = 0; i < nTheta_; ++i) {
    const Real thDeg = theta_[i] / cDegToRad;
    const Real sn = sinTheta_[i];
    for (int j = 0; j < nPhi_; ++j) {
      const int k = idx2d(i, j, nPhi_);
      const Real phDeg = phi_[j] / cDegToRad;
      const Real eth = -dPhiDTheta_[k] * invRI * Real(1.0e3); // mV/m
      const Real eph = (sn > Real(1.0e-8))
                           ? (-dPhiDPhi_[k] * invRI / sn * Real(1.0e3))
                           : Real(0.0);
      out << std::fixed << std::setprecision(4) << thDeg << " " << phDeg << " "
          << std::scientific << std::setprecision(6) << jr_[k] * Real(1.0e6)
          << " " << sigmaP_[k] << " " << sigmaH_[k] << " "
          << phiPot_[k] * Real(1.0e-3) << " " << eth << " " << eph << "\n";
    }
  }
}

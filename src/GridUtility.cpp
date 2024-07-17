#include <iostream>

#include "Constants.h"
#include "GridUtility.h"
#include "Timer.h"

using namespace amrex;

void lap_node_to_node(const MultiFab& srcMF, MultiFab& dstMF,
                      const DistributionMapping dm, const Geometry& gm) {
  const Real* invDx = gm.InvCellSize();

  BoxArray centerBA =
      convert(srcMF.boxArray(), IntVect{ AMREX_D_DECL(0, 0, 0) });

  // Need and just need 1 ghost cell layer.
  MultiFab centerMF(centerBA, dm, 3, 1);
  centerMF.setVal(0.0);

  for (int i = 0; i < srcMF.nComp(); ++i) {
    MultiFab srcAliasMF(srcMF, make_alias, i, 1);
    grad_node_to_center(srcAliasMF, centerMF, invDx);

    // centerMF.FillBoundary(gm.periodicity());
    MultiFab dstAliasMF(dstMF, make_alias, i, 1);
    div_center_to_node(centerMF, dstAliasMF, invDx);
  }
}

void grad_node_to_center(const MultiFab& nodeMF, MultiFab& centerMF,
                         const Real* invDx) {
  timing_func("grad_node_to_center");

  for (MFIter mfi(centerMF, doTiling); mfi.isValid(); ++mfi) {
    Box box = mfi.validbox();
    box.grow(1);

    const Array4<Real>& center = centerMF[mfi].array();
    const Array4<Real const>& node = nodeMF[mfi].array();

    ParallelFor(box, [&](int i, int j, int k) {
      int kp1 = nDim > 2 ? k + 1 : k;
      center(i, j, k, ix_) =
          0.25 * invDx[ix_] *
          (node(i + 1, j, k) - node(i, j, k) + node(i + 1, j, kp1) -
           node(i, j, kp1) + node(i + 1, j + 1, k) - node(i, j + 1, k) +
           node(i + 1, j + 1, kp1) - node(i, j + 1, kp1));
      center(i, j, k, iy_) =
          0.25 * invDx[iy_] *
          (node(i, j + 1, k) - node(i, j, k) + node(i, j + 1, kp1) -
           node(i, j, kp1) + node(i + 1, j + 1, k) - node(i + 1, j, k) +
           node(i + 1, j + 1, kp1) - node(i + 1, j, kp1));
      center(i, j, k, iz_) =
          0.25 * invDx[iz_] *
          (node(i, j, kp1) - node(i, j, k) + node(i + 1, j, kp1) -
           node(i + 1, j, k) + node(i, j + 1, kp1) - node(i, j + 1, k) +
           node(i + 1, j + 1, kp1) - node(i + 1, j + 1, k));
    });
  }
}

void grad_center_to_node(const MultiFab& centerMF, MultiFab& nodeMF,
                         const Real* invDx) {

  for (MFIter mfi(nodeMF, doTiling); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();

    const Array4<Real>& node = nodeMF[mfi].array();
    const Array4<Real const>& center = centerMF[mfi].array();

    ParallelFor(box, [&](int i, int j, int k) {
      int km1 = nDim > 2 ? k - 1 : k;

      node(i, j, k, ix_) = 0.25 * invDx[ix_] *
                           (center(i, j, k) - center(i - 1, j, k) +
                            center(i, j, km1) - center(i - 1, j, km1) +
                            center(i, j - 1, k) - center(i - 1, j - 1, k) +
                            center(i, j - 1, km1) - center(i - 1, j - 1, km1));
      node(i, j, k, iy_) = 0.25 * invDx[iy_] *
                           (center(i, j, k) - center(i, j - 1, k) +
                            center(i, j, km1) - center(i, j - 1, km1) +
                            center(i - 1, j, k) - center(i - 1, j - 1, k) +
                            center(i - 1, j, km1) - center(i - 1, j - 1, km1));
      node(i, j, k, iz_) =
          0.25 * invDx[iz_] *
          (center(i, j, k) - center(i, j, km1) + center(i - 1, j, k) -
           center(i - 1, j, km1) + center(i, j - 1, k) - center(i, j - 1, km1) +
           center(i - 1, j - 1, k) - center(i - 1, j - 1, km1));
    });
  }
}

void div_center_to_node(const MultiFab& centerMF, MultiFab& nodeMF,
                        const Real* invDx) {

  for (MFIter mfi(nodeMF, doTiling); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();

    const Array4<Real const>& center = centerMF[mfi].array();
    const Array4<Real>& node = nodeMF[mfi].array();

    ParallelFor(box, [&](int i, int j, int k) {
      int km1 = nDim > 2 ? k - 1 : k;
      const Real compX =
          0.25 * invDx[ix_] *
          (center(i, j, k, ix_) - center(i - 1, j, k, ix_) +
           center(i, j, km1, ix_) - center(i - 1, j, km1, ix_) +
           center(i, j - 1, k, ix_) - center(i - 1, j - 1, k, ix_) +
           center(i, j - 1, km1, ix_) - center(i - 1, j - 1, km1, ix_));

      const Real compY =
          0.25 * invDx[iy_] *
          (center(i, j, k, iy_) - center(i, j - 1, k, iy_) +
           center(i, j, km1, iy_) - center(i, j - 1, km1, iy_) +
           center(i - 1, j, k, iy_) - center(i - 1, j - 1, k, iy_) +
           center(i - 1, j, km1, iy_) - center(i - 1, j - 1, km1, iy_));

      const Real compZ =
          0.25 * invDx[iz_] *
          (center(i, j, k, iz_) - center(i, j, km1, iz_) +
           center(i - 1, j, k, iz_) - center(i - 1, j, km1, iz_) +
           center(i, j - 1, k, iz_) - center(i, j - 1, km1, iz_) +
           center(i - 1, j - 1, k, iz_) - center(i - 1, j - 1, km1, iz_));
      node(i, j, k) = compX + compY + compZ;
    });
  }
}

void div_node_to_center(const MultiFab& nodeMF, MultiFab& centerMF,
                        const Real* invDx) {

  for (MFIter mfi(centerMF, doTiling); mfi.isValid(); ++mfi) {
    const Box& box = mfi.fabbox();

    const Array4<Real const>& node = nodeMF[mfi].array();
    const Array4<Real>& center = centerMF[mfi].array();

    ParallelFor(box, [&](int i, int j, int k) {
      int kp1 = nDim > 2 ? k + 1 : k;

      const Real compX =
          0.25 * invDx[ix_] *
          (node(i + 1, j, k, ix_) - node(i, j, k, ix_) +
           node(i + 1, j, kp1, ix_) - node(i, j, kp1, ix_) +
           node(i + 1, j + 1, k, ix_) - node(i, j + 1, k, ix_) +
           node(i + 1, j + 1, kp1, ix_) - node(i, j + 1, kp1, ix_));

      const Real compY =
          0.25 * invDx[iy_] *
          (node(i, j + 1, k, iy_) - node(i, j, k, iy_) +
           node(i, j + 1, kp1, iy_) - node(i, j, kp1, iy_) +
           node(i + 1, j + 1, k, iy_) - node(i + 1, j, k, iy_) +
           node(i + 1, j + 1, kp1, iy_) - node(i + 1, j, kp1, iy_));

      const Real compZ =
          0.25 * invDx[iz_] *
          (node(i, j, kp1, iz_) - node(i, j, k, iz_) +
           node(i + 1, j, kp1, iz_) - node(i + 1, j, k, iz_) +
           node(i, j + 1, kp1, iz_) - node(i, j + 1, k, iz_) +
           node(i + 1, j + 1, kp1, iz_) - node(i + 1, j + 1, k, iz_));
      center(i, j, k) = compX + compY + compZ;
    });
  }
}

void div_center_to_center(const MultiFab& srcMF, MultiFab& dstMF,
                          const Real* invDx) {

  for (MFIter mfi(dstMF, doTiling); mfi.isValid(); ++mfi) {
    Box box = mfi.validbox();
    box.grow(1);

    const Array4<Real const>& srcArr = srcMF[mfi].array();
    const Array4<Real>& dstArr = dstMF[mfi].array();

    ParallelFor(box, [&](int i, int j, int k) {
      int km1 = nDim > 2 ? k - 1 : k;
      int kp1 = nDim > 2 ? k + 1 : k;
      Real compX = 0;
      for (int jj = -1; jj < 2; ++jj)
        for (int kk = -1; kk < 2; ++kk) {
          int k0 = nDim > 2 ? k + kk : 0;
          compX +=
              srcArr(i + 1, j + jj, k0, ix_) - srcArr(i - 1, j + jj, k0, ix_);
        }
      compX *= 0.5 * invDx[ix_];

      Real compY = 0;
      for (int ii = -1; ii < 2; ++ii)
        for (int kk = -1; kk < 2; ++kk) {
          int k0 = nDim > 2 ? k + kk : 0;
          compY +=
              srcArr(i + ii, j + 1, k0, iy_) - srcArr(i + ii, j - 1, k0, iy_);
        }
      compY *= 0.5 * invDx[iy_];

      Real compZ = 0;
      for (int ii = -1; ii < 2; ++ii)
        for (int jj = -1; jj < 2; ++jj) {
          compZ += srcArr(i + ii, j + jj, kp1, iz_) -
                   srcArr(i + ii, j + jj, km1, iz_);
        }
      compZ *= 0.5 * invDx[iz_];

      dstArr(i, j, k) = (compX + compY + compZ) / 9;
    });
  }
}

void print_MultiFab(const MultiFab& data, std::string tag, Geometry& gm,
                    int nshift) {
  AllPrint() << "-----" << tag << " begin-----" << std::endl;
  Real sum = 0;
  Real sum2 = 0;

  bool isCenter = data.ixType().cellCentered();
  bool isNode = !isCenter;

  const Box& gbx = convert(gm.Domain(), data.boxArray().ixType());

  for (MFIter mfi(data); mfi.isValid(); ++mfi) {
    const FArrayBox& fab = data[mfi];
    const Box& box = mfi.validbox();
    const Array4<const Real>& data = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    // Avoid double counting the share edges.
    int iMax = hi.x, jMax = hi.y, kMax = hi.z;
    int iMin = lo.x, jMin = lo.y, kMin = lo.z;

    if (isNode) {
      // Avoid double counting the shared edges.
      iMax--;
      jMax--;
      kMax--;

      if ((!gm.isPeriodic(ix_)) && gbx.bigEnd(ix_) == hi.x)
        iMax++;
      if ((!gm.isPeriodic(iy_)) && gbx.bigEnd(iy_) == hi.y)
        jMax++;
      if ((!gm.isPeriodic(iz_)) && gbx.bigEnd(iz_) == hi.z)
        kMax++;
    }

    if (!gm.isPeriodic(ix_) && gbx.bigEnd(ix_) == hi.x) {
      iMax += nshift;
    }

    if (!gm.isPeriodic(iy_) && gbx.bigEnd(iy_) == hi.y) {
      jMax += nshift;
    }

    if (!gm.isPeriodic(iz_) && gbx.bigEnd(iz_) == hi.z) {
      kMax += nshift;
    }

    if (!gm.isPeriodic(ix_) && gbx.smallEnd(ix_) == lo.x) {
      iMin -= nshift;
    }

    if (!gm.isPeriodic(iy_) && gbx.smallEnd(iy_) == lo.y) {
      jMin -= nshift;
    }

    if (!gm.isPeriodic(iz_) && gbx.smallEnd(iz_) == lo.z) {
      kMin -= nshift;
    }

    for (int i = iMin; i <= iMax; ++i)
      for (int j = jMin; j <= jMax; ++j)
        for (int k = kMin; k <= kMax; ++k)
          for (int iVar = 0; iVar < data.nComp(); iVar++) {
            AllPrint() << " i = " << i << " j = " << j << " k = " << k
                       << " iVar = " << iVar
                       << " data = " << data(i, j, k, iVar) << std::endl;
            sum += data(i, j, k, iVar);
            sum2 += pow(data(i, j, k, iVar), 2);
          }
  }
  AllPrint() << "sum = " << sum << " sum2 = " << sqrt(sum2)
             << " on proc = " << ParallelDescriptor::MyProc() << std::endl;
  AllPrint() << "-----" << tag << " end-----" << std::endl;
}

void print_MultiFab(const MultiFab& data, std::string tag, const int iVarStart,
                    const int iVarEnd, int nshift) {
  AllPrint() << "-----" << tag << " begin-----" << std::endl;
  Real sum = 0;
  Real sum2 = 0;

  for (MFIter mfi(data); mfi.isValid(); ++mfi) {
    const FArrayBox& fab = data[mfi];
    const Box& box = mfi.validbox();
    const Array4<const Real>& data = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    AllPrint() << "--------------------------------" << std::endl;
    for (int i = lo.x - nshift; i <= hi.x + nshift; ++i)
      for (int j = lo.y - nshift; j <= hi.y + nshift; ++j)
        for (int k = lo.z; k <= hi.z; ++k)
          for (int iVar = iVarStart; iVar < iVarEnd; iVar++) {
            AllPrint() << " i = " << i << " j = " << j << " k = " << k
                       << " iVar = " << iVar
                       << " data = " << data(i, j, k, iVar) << std::endl;
            sum += data(i, j, k, iVar);
            sum2 += pow(data(i, j, k, iVar), 2);
          }
  }
  AllPrint() << "sum = " << sum << " sum2 = " << sqrt(sum2)
             << " on proc = " << ParallelDescriptor::MyProc() << std::endl;
  AllPrint() << "-----" << tag << " end-----" << std::endl;
}

void print_MultiFab(const iMultiFab& data, std::string tag, int nshift) {
  AllPrint() << "-----" << tag << " begin-----" << std::endl;
  Real sum = 0;
  Real sum2 = 0;

  for (MFIter mfi(data); mfi.isValid(); ++mfi) {
    const IArrayBox& fab = data[mfi];
    const Box& box = mfi.validbox();
    const Array4<const int>& dataArr = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    AllPrint() << "box = " << box << std::endl;

    for (int i = lo.x - nshift; i <= hi.x + nshift; ++i)
      for (int j = lo.y - nshift; j <= hi.y + nshift; ++j)
        for (int k = lo.z - nshift; k <= hi.z + nshift; ++k)
          for (int iVar = 0; iVar < data.nComp(); iVar++) {
            AllPrint() << " i = " << i << " j = " << j << " k = " << k
                       << " iVar = " << iVar
                       << " data = " << dataArr(i, j, k, iVar) << std::endl;
            sum += dataArr(i, j, k, iVar);
            sum2 += pow(dataArr(i, j, k, iVar), 2);
          }
  }
  AllPrint() << "sum = " << sum << " sum2 = " << sqrt(sum2)
             << " on proc = " << ParallelDescriptor::MyProc() << std::endl;
  AllPrint() << "-----" << tag << " end-----" << std::endl;
}

void curl_center_to_node(const MultiFab& centerMF, MultiFab& nodeMF,
                         const Real* invDx) {
  Real cZDY, cYDZ;

  Real cXDZ, cZDX;
  Real cYDX, cXDY;

  for (MFIter mfi(nodeMF, doTiling); mfi.isValid(); ++mfi) {
    Box box = mfi.fabbox();
    box.grow(-1);

    const Array4<Real>& nodeArr = nodeMF[mfi].array();
    const Array4<Real const>& centerArr = centerMF[mfi].array();

    ParallelFor(box, [&](int i, int j, int k) {
      int km1 = nDim > 2 ? k - 1 : k;
      cZDY =
          0.25 * invDx[iy_] *
          (centerArr(i, j, k, iz_) - centerArr(i, j - 1, k, iz_) +
           centerArr(i, j, km1, iz_) - centerArr(i, j - 1, km1, iz_) +
           centerArr(i - 1, j, k, iz_) - centerArr(i - 1, j - 1, k, iz_) +
           centerArr(i - 1, j, km1, iz_) - centerArr(i - 1, j - 1, km1, iz_));
      cYDZ =
          0.25 * invDx[iz_] *
          (centerArr(i, j, k, iy_) - centerArr(i, j, km1, iy_) +
           centerArr(i - 1, j, k, iy_) - centerArr(i - 1, j, km1, iy_) +
           centerArr(i, j - 1, k, iy_) - centerArr(i, j - 1, km1, iy_) +
           centerArr(i - 1, j - 1, k, iy_) - centerArr(i - 1, j - 1, km1, iy_));
      // curl - Y
      cXDZ =
          0.25 * invDx[iz_] *
          (centerArr(i, j, k, ix_) - centerArr(i, j, km1, ix_) +
           centerArr(i - 1, j, k, ix_) - centerArr(i - 1, j, km1, ix_) +
           centerArr(i, j - 1, k, ix_) - centerArr(i, j - 1, km1, ix_) +
           centerArr(i - 1, j - 1, k, ix_) - centerArr(i - 1, j - 1, km1, ix_));

      cZDX =
          0.25 * invDx[ix_] *
          (centerArr(i, j, k, iz_) - centerArr(i - 1, j, k, iz_) +
           centerArr(i, j, km1, iz_) - centerArr(i - 1, j, km1, iz_) +
           centerArr(i, j - 1, k, iz_) - centerArr(i - 1, j - 1, k, iz_) +
           centerArr(i, j - 1, km1, iz_) - centerArr(i - 1, j - 1, km1, iz_));

      // curl - Z
      cYDX =
          0.25 * invDx[ix_] *
          (centerArr(i, j, k, iy_) - centerArr(i - 1, j, k, iy_) +
           centerArr(i, j, km1, iy_) - centerArr(i - 1, j, km1, iy_) +
           centerArr(i, j - 1, k, iy_) - centerArr(i - 1, j - 1, k, iy_) +
           centerArr(i, j - 1, km1, iy_) - centerArr(i - 1, j - 1, km1, iy_));

      cXDY =
          0.25 * invDx[iy_] *
          (centerArr(i, j, k, ix_) - centerArr(i, j - 1, k, ix_) +
           centerArr(i, j, km1, ix_) - centerArr(i, j - 1, km1, ix_) +
           centerArr(i - 1, j, k, ix_) - centerArr(i - 1, j - 1, k, ix_) +
           centerArr(i - 1, j, km1, ix_) - centerArr(i - 1, j - 1, km1, ix_));

      nodeArr(i, j, k, ix_) = cZDY - cYDZ;
      nodeArr(i, j, k, iy_) = cXDZ - cZDX;
      nodeArr(i, j, k, iz_) = cYDX - cXDY;
    });
  }
}

void curl_node_to_center(const MultiFab& nodeMF, MultiFab& centerMF,
                         const Real* invDx) {
  Real cZDY, cYDZ;

  Real cXDZ, cZDX;
  Real cYDX, cXDY;

  for (MFIter mfi(centerMF, doTiling); mfi.isValid(); ++mfi) {
    const Box& box = mfi.fabbox();
    const Array4<Real>& centerArr = centerMF[mfi].array();
    const Array4<Real const>& nodeArr = nodeMF[mfi].array();

    ParallelFor(box, [&](int i, int j, int k) {
      int kp1 = nDim > 2 ? k + 1 : k;
      cZDY = 0.25 * invDx[iy_] *
             (nodeArr(i, j + 1, k, iz_) - nodeArr(i, j, k, iz_) +
              nodeArr(i, j + 1, kp1, iz_) - nodeArr(i, j, kp1, iz_) +
              nodeArr(i + 1, j + 1, k, iz_) - nodeArr(i + 1, j, k, iz_) +
              nodeArr(i + 1, j + 1, kp1, iz_) - nodeArr(i + 1, j, kp1, iz_));
      cYDZ = 0.25 * invDx[iz_] *
             (nodeArr(i, j, kp1, iy_) - nodeArr(i, j, k, iy_) +
              nodeArr(i + 1, j, kp1, iy_) - nodeArr(i + 1, j, k, iy_) +
              nodeArr(i, j + 1, kp1, iy_) - nodeArr(i, j + 1, k, iy_) +
              nodeArr(i + 1, j + 1, kp1, iy_) - nodeArr(i + 1, j + 1, k, iy_));
      // curl - Y
      cXDZ = 0.25 * invDx[iz_] *
             (nodeArr(i, j, kp1, ix_) - nodeArr(i, j, k, ix_) +
              nodeArr(i + 1, j, kp1, ix_) - nodeArr(i + 1, j, k, ix_) +
              nodeArr(i, j + 1, kp1, ix_) - nodeArr(i, j + 1, k, ix_) +
              nodeArr(i + 1, j + 1, kp1, ix_) - nodeArr(i + 1, j + 1, k, ix_));

      cZDX = 0.25 * invDx[ix_] *
             (nodeArr(i + 1, j, k, iz_) - nodeArr(i, j, k, iz_) +
              nodeArr(i + 1, j, kp1, iz_) - nodeArr(i, j, kp1, iz_) +
              nodeArr(i + 1, j + 1, k, iz_) - nodeArr(i, j + 1, k, iz_) +
              nodeArr(i + 1, j + 1, kp1, iz_) - nodeArr(i, j + 1, kp1, iz_));

      // curl - Z
      cYDX = 0.25 * invDx[ix_] *
             (nodeArr(i + 1, j, k, iy_) - nodeArr(i, j, k, iy_) +
              nodeArr(i + 1, j, kp1, iy_) - nodeArr(i, j, kp1, iy_) +
              nodeArr(i + 1, j + 1, k, iy_) - nodeArr(i, j + 1, k, iy_) +
              nodeArr(i + 1, j + 1, kp1, iy_) - nodeArr(i, j + 1, kp1, iy_));

      cXDY = 0.25 * invDx[iy_] *
             (nodeArr(i, j + 1, k, ix_) - nodeArr(i, j, k, ix_) +
              nodeArr(i, j + 1, kp1, ix_) - nodeArr(i, j, kp1, ix_) +
              nodeArr(i + 1, j + 1, k, ix_) - nodeArr(i + 1, j, k, ix_) +
              nodeArr(i + 1, j + 1, kp1, ix_) - nodeArr(i + 1, j, kp1, ix_));

      centerArr(i, j, k, ix_) = cZDY - cYDZ;
      centerArr(i, j, k, iy_) = cXDZ - cZDX;
      centerArr(i, j, k, iz_) = cYDX - cXDY;
    });
  }
}

void average_center_to_node(const MultiFab& centerMF, MultiFab& nodeMF) {
  for (MFIter mfi(nodeMF, doTiling); mfi.isValid(); ++mfi) {
    Box box = mfi.fabbox();
    box.grow(-1);

    const Array4<Real>& nodeArr = nodeMF[mfi].array();
    const Array4<Real const>& centerArr = centerMF[mfi].array();

    ParallelFor(box, centerMF.nComp(), [&](int i, int j, int k, int iVar) {
      int km1 = nDim > 2 ? k - 1 : k;
      nodeArr(i, j, k, iVar) =
          0.125 *
          (centerArr(i - 1, j - 1, km1, iVar) +
           centerArr(i - 1, j - 1, k, iVar) + centerArr(i - 1, j, km1, iVar) +
           centerArr(i - 1, j, k, iVar) + centerArr(i, j - 1, km1, iVar) +
           centerArr(i, j - 1, k, iVar) + centerArr(i, j, km1, iVar) +
           centerArr(i, j, k, iVar));
    });
  }
}
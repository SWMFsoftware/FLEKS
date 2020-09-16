#include <iostream>

#include "Utility.h"
#include "Constants.h"
#include "Timer.h"

using namespace std;
using namespace amrex;

void lap_node_to_node(const amrex::MultiFab& srcMF, amrex::MultiFab& dstMF,
                      const amrex::DistributionMapping dm,
                      const amrex::Geometry& geom,
                      const amrex::iMultiFab& status) {
  const amrex::Real* invDx = geom.InvCellSize();

  BoxArray centerBA =
      convert(srcMF.boxArray(), IntVect{ AMREX_D_DECL(0, 0, 0) });

  // Need and just need 1 ghost cell layer.
  MultiFab centerMF(centerBA, dm, 3, 1);
  centerMF.setVal(0.0);

  for (int i = 0; i < srcMF.nComp(); i++) {
    MultiFab srcAliasMF(srcMF, amrex::make_alias, i, 1);
    grad_node_to_center(srcAliasMF, centerMF, invDx, status);

    //centerMF.FillBoundary(geom.periodicity());
    MultiFab dstAliasMF(dstMF, amrex::make_alias, i, 1);
    div_center_to_node(centerMF, dstAliasMF, invDx);
  }
}

void grad_node_to_center(const amrex::MultiFab& nodeMF,
                         amrex::MultiFab& centerMF, const amrex::Real* invDx,
                         const amrex::iMultiFab& status) {
  timing_func("grad_node_to_center");

  for (amrex::MFIter mfi(centerMF, doTiling); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int imin = lo.x - 1, jmin = lo.y - 1, kmin = lo.z - 1;
    int imax = hi.x + 1, jmax = hi.y + 1, kmax = hi.z + 1;
    
    const amrex::Array4<amrex::Real>& center = centerMF[mfi].array();
    const amrex::Array4<amrex::Real const>& node = nodeMF[mfi].array();

    for (int k = kmin; k <= kmax; ++k)
      for (int j = jmin; j <= jmax; ++j)
        for (int i = imin; i <= imax; ++i) {
          center(i, j, k, ix_) =
              0.25 * invDx[ix_] *
              (node(i + 1, j, k) - node(i, j, k) + node(i + 1, j, k + 1) -
               node(i, j, k + 1) + node(i + 1, j + 1, k) - node(i, j + 1, k) +
               node(i + 1, j + 1, k + 1) - node(i, j + 1, k + 1));
          center(i, j, k, iy_) =
              0.25 * invDx[iy_] *
              (node(i, j + 1, k) - node(i, j, k) + node(i, j + 1, k + 1) -
               node(i, j, k + 1) + node(i + 1, j + 1, k) - node(i + 1, j, k) +
               node(i + 1, j + 1, k + 1) - node(i + 1, j, k + 1));
          center(i, j, k, iz_) =
              0.25 * invDx[iz_] *
              (node(i, j, k + 1) - node(i, j, k) + node(i + 1, j, k + 1) -
               node(i + 1, j, k) + node(i, j + 1, k + 1) - node(i, j + 1, k) +
               node(i + 1, j + 1, k + 1) - node(i + 1, j + 1, k));
        }
  }
}

void grad_center_to_node(const amrex::MultiFab& centerMF,
                         amrex::MultiFab& nodeMF, const amrex::Real* invDx) {

  for (amrex::MFIter mfi(nodeMF, doTiling); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const amrex::Array4<amrex::Real>& node = nodeMF[mfi].array();
    const amrex::Array4<amrex::Real const>& center = centerMF[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {

          node(i, j, k, ix_) =
              0.25 * invDx[ix_] *
              (center(i, j, k) - center(i - 1, j, k) + center(i, j, k - 1) -
               center(i - 1, j, k - 1) + center(i, j - 1, k) -
               center(i - 1, j - 1, k) + center(i, j - 1, k - 1) -
               center(i - 1, j - 1, k - 1));
          node(i, j, k, iy_) =
              0.25 * invDx[iy_] *
              (center(i, j, k) - center(i, j - 1, k) + center(i, j, k - 1) -
               center(i, j - 1, k - 1) + center(i - 1, j, k) -
               center(i - 1, j - 1, k) + center(i - 1, j, k - 1) -
               center(i - 1, j - 1, k - 1));
          node(i, j, k, iz_) =
              0.25 * invDx[iz_] *
              (center(i, j, k) - center(i, j, k - 1) + center(i - 1, j, k) -
               center(i - 1, j, k - 1) + center(i, j - 1, k) -
               center(i, j - 1, k - 1) + center(i - 1, j - 1, k) -
               center(i - 1, j - 1, k - 1));
        }
  }
}

void div_center_to_node(const amrex::MultiFab& centerMF,
                        amrex::MultiFab& nodeMF, const amrex::Real* invDx) {

  for (amrex::MFIter mfi(nodeMF, doTiling); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const amrex::Array4<amrex::Real const>& center = centerMF[mfi].array();
    const amrex::Array4<amrex::Real>& node = nodeMF[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {

          const Real compX =
              0.25 * invDx[ix_] *
              (center(i, j, k, ix_) - center(i - 1, j, k, ix_) +
               center(i, j, k - 1, ix_) - center(i - 1, j, k - 1, ix_) +
               center(i, j - 1, k, ix_) - center(i - 1, j - 1, k, ix_) +
               center(i, j - 1, k - 1, ix_) - center(i - 1, j - 1, k - 1, ix_));

          const Real compY =
              0.25 * invDx[iy_] *
              (center(i, j, k, iy_) - center(i, j - 1, k, iy_) +
               center(i, j, k - 1, iy_) - center(i, j - 1, k - 1, iy_) +
               center(i - 1, j, k, iy_) - center(i - 1, j - 1, k, iy_) +
               center(i - 1, j, k - 1, iy_) - center(i - 1, j - 1, k - 1, iy_));

          const Real compZ =
              0.25 * invDx[iz_] *
              (center(i, j, k, iz_) - center(i, j, k - 1, iz_) +
               center(i - 1, j, k, iz_) - center(i - 1, j, k - 1, iz_) +
               center(i, j - 1, k, iz_) - center(i, j - 1, k - 1, iz_) +
               center(i - 1, j - 1, k, iz_) - center(i - 1, j - 1, k - 1, iz_));
          node(i, j, k) = compX + compY + compZ;
        }
  }
}

void div_node_to_center(const amrex::MultiFab& nodeMF,
                        amrex::MultiFab& centerMF, const amrex::Real* invDx) {

  for (amrex::MFIter mfi(centerMF, doTiling); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.fabbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const amrex::Array4<amrex::Real const>& node = nodeMF[mfi].array();
    const amrex::Array4<amrex::Real>& center = centerMF[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {

          const Real compX =
              0.25 * invDx[ix_] *
              (node(i + 1, j, k, ix_) - node(i, j, k, ix_) +
               node(i + 1, j, k + 1, ix_) - node(i, j, k + 1, ix_) +
               node(i + 1, j + 1, k, ix_) - node(i, j + 1, k, ix_) +
               node(i + 1, j + 1, k + 1, ix_) - node(i, j + 1, k + 1, ix_));

          const Real compY =
              0.25 * invDx[iy_] *
              (node(i, j + 1, k, iy_) - node(i, j, k, iy_) +
               node(i, j + 1, k + 1, iy_) - node(i, j, k + 1, iy_) +
               node(i + 1, j + 1, k, iy_) - node(i + 1, j, k, iy_) +
               node(i + 1, j + 1, k + 1, iy_) - node(i + 1, j, k + 1, iy_));

          const Real compZ =
              0.25 * invDx[iz_] *
              (node(i, j, k + 1, iz_) - node(i, j, k, iz_) +
               node(i + 1, j, k + 1, iz_) - node(i + 1, j, k, iz_) +
               node(i, j + 1, k + 1, iz_) - node(i, j + 1, k, iz_) +
               node(i + 1, j + 1, k + 1, iz_) - node(i + 1, j + 1, k, iz_));
          center(i, j, k) = compX + compY + compZ;
        }
  }
}

void div_center_to_center(const amrex::MultiFab& srcMF, amrex::MultiFab& dstMF,
                          const amrex::Real* invDx) {

  for (amrex::MFIter mfi(dstMF, doTiling); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const amrex::Array4<amrex::Real const>& srcArr = srcMF[mfi].array();
    const amrex::Array4<amrex::Real>& dstArr = dstMF[mfi].array();

    for (int k = lo.z - 1; k <= hi.z + 1; ++k)
      for (int j = lo.y - 1; j <= hi.y + 1; ++j)
        for (int i = lo.x - 1; i <= hi.x + 1; ++i) {
          Real compX = 0;
          for (int jj = -1; jj < 2; jj++)
            for (int kk = -1; kk < 2; kk++) {
              compX += srcArr(i + 1, j + jj, k + kk, ix_) -
                       srcArr(i - 1, j + jj, k + kk, ix_);
            }
          compX *= 0.5 * invDx[ix_];

          Real compY = 0;
          for (int ii = -1; ii < 2; ii++)
            for (int kk = -1; kk < 2; kk++) {
              compY += srcArr(i + ii, j + 1, k + kk, iy_) -
                       srcArr(i + ii, j - 1, k + kk, iy_);
            }
          compY *= 0.5 * invDx[iy_];

          Real compZ = 0;
          for (int ii = -1; ii < 2; ii++)
            for (int jj = -1; jj < 2; jj++) {
              compZ += srcArr(i + ii, j + jj, k + 1, iz_) -
                       srcArr(i + ii, j + jj, k - 1, iz_);
            }
          compZ *= 0.5 * invDx[iz_];

          dstArr(i, j, k) = (compX + compY + compZ) / 9;
        }
  }
}

void print_MultiFab(const amrex::MultiFab& data, std::string tag,
                    Geometry& geom, int nshift) {
  AllPrint() << "-----" << tag << " begin-----" << std::endl;
  Real sum = 0;
  Real sum2 = 0;

  bool isCenter = data.ixType().cellCentered();
  bool isNode = !isCenter;

  const Box& gbx = convert(geom.Domain(), data.boxArray().ixType());

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

      if ((!geom.isPeriodic(ix_)) && gbx.bigEnd(ix_) == hi.x)
        iMax++;
      if ((!geom.isPeriodic(iy_)) && gbx.bigEnd(iy_) == hi.y)
        jMax++;
      if ((!geom.isPeriodic(iz_)) && gbx.bigEnd(iz_) == hi.z)
        kMax++;
    }

    if (!geom.isPeriodic(ix_) && gbx.bigEnd(ix_) == hi.x) {
      iMax += nshift;
    }

    if (!geom.isPeriodic(iy_) && gbx.bigEnd(iy_) == hi.y) {
      jMax += nshift;
    }

    if (!geom.isPeriodic(iz_) && gbx.bigEnd(iz_) == hi.z) {
      kMax += nshift;
    }

    if (!geom.isPeriodic(ix_) && gbx.smallEnd(ix_) == lo.x) {
      iMin -= nshift;
    }

    if (!geom.isPeriodic(iy_) && gbx.smallEnd(iy_) == lo.y) {
      jMin -= nshift;
    }

    if (!geom.isPeriodic(iz_) && gbx.smallEnd(iz_) == lo.z) {
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

void print_MultiFab(const amrex::MultiFab& data, std::string tag, int nshift) {
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

void print_MultiFab(const amrex::iMultiFab& data, std::string tag, int nshift) {
  AllPrint() << "-----" << tag << " begin-----" << std::endl;
  Real sum = 0;
  Real sum2 = 0;

  for (MFIter mfi(data); mfi.isValid(); ++mfi) {
    const IArrayBox& fab = data[mfi];
    const Box& box = mfi.validbox();
    const Array4<const int>& dataArr = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    AllPrint()<<"box = "<<box<<std::endl;

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
    const Box& box = mfi.fabbox();
    const Array4<Real>& nodeArr = nodeMF[mfi].array();
    const Array4<Real const>& centerArr = centerMF[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z + 1; k <= hi.z - 1; ++k)
      for (int j = lo.y + 1; j <= hi.y - 1; ++j)
        for (int i = lo.x + 1; i <= hi.x - 1; ++i) {
          cZDY =
              0.25 * invDx[iy_] *
              (centerArr(i, j, k, iz_) - centerArr(i, j - 1, k, iz_) +
               centerArr(i, j, k - 1, iz_) - centerArr(i, j - 1, k - 1, iz_) +
               centerArr(i - 1, j, k, iz_) - centerArr(i - 1, j - 1, k, iz_) +
               centerArr(i - 1, j, k - 1, iz_) -
               centerArr(i - 1, j - 1, k - 1, iz_));
          cYDZ =
              0.25 * invDx[iz_] *
              (centerArr(i, j, k, iy_) - centerArr(i, j, k - 1, iy_) +
               centerArr(i - 1, j, k, iy_) - centerArr(i - 1, j, k - 1, iy_) +
               centerArr(i, j - 1, k, iy_) - centerArr(i, j - 1, k - 1, iy_) +
               centerArr(i - 1, j - 1, k, iy_) -
               centerArr(i - 1, j - 1, k - 1, iy_));
          // curl - Y
          cXDZ =
              0.25 * invDx[iz_] *
              (centerArr(i, j, k, ix_) - centerArr(i, j, k - 1, ix_) +
               centerArr(i - 1, j, k, ix_) - centerArr(i - 1, j, k - 1, ix_) +
               centerArr(i, j - 1, k, ix_) - centerArr(i, j - 1, k - 1, ix_) +
               centerArr(i - 1, j - 1, k, ix_) -
               centerArr(i - 1, j - 1, k - 1, ix_));

          cZDX =
              0.25 * invDx[ix_] *
              (centerArr(i, j, k, iz_) - centerArr(i - 1, j, k, iz_) +
               centerArr(i, j, k - 1, iz_) - centerArr(i - 1, j, k - 1, iz_) +
               centerArr(i, j - 1, k, iz_) - centerArr(i - 1, j - 1, k, iz_) +
               centerArr(i, j - 1, k - 1, iz_) -
               centerArr(i - 1, j - 1, k - 1, iz_));

          // curl - Z
          cYDX =
              0.25 * invDx[ix_] *
              (centerArr(i, j, k, iy_) - centerArr(i - 1, j, k, iy_) +
               centerArr(i, j, k - 1, iy_) - centerArr(i - 1, j, k - 1, iy_) +
               centerArr(i, j - 1, k, iy_) - centerArr(i - 1, j - 1, k, iy_) +
               centerArr(i, j - 1, k - 1, iy_) -
               centerArr(i - 1, j - 1, k - 1, iy_));

          cXDY =
              0.25 * invDx[iy_] *
              (centerArr(i, j, k, ix_) - centerArr(i, j - 1, k, ix_) +
               centerArr(i, j, k - 1, ix_) - centerArr(i, j - 1, k - 1, ix_) +
               centerArr(i - 1, j, k, ix_) - centerArr(i - 1, j - 1, k, ix_) +
               centerArr(i - 1, j, k - 1, ix_) -
               centerArr(i - 1, j - 1, k - 1, ix_));

          nodeArr(i, j, k, ix_) = cZDY - cYDZ;
          nodeArr(i, j, k, iy_) = cXDZ - cZDX;
          nodeArr(i, j, k, iz_) = cYDX - cXDY;
        }
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

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          // Needs to be improved. --Yuxi
          cZDY = 0.25 * invDx[iy_] *
                 (nodeArr(i, j + 1, k, iz_) - nodeArr(i, j, k, iz_) +
                  nodeArr(i, j + 1, k + 1, iz_) - nodeArr(i, j, k + 1, iz_) +
                  nodeArr(i + 1, j + 1, k, iz_) - nodeArr(i + 1, j, k, iz_) +
                  nodeArr(i + 1, j + 1, k + 1, iz_) -
                  nodeArr(i + 1, j, k + 1, iz_));
          cYDZ = 0.25 * invDx[iz_] *
                 (nodeArr(i, j, k + 1, iy_) - nodeArr(i, j, k, iy_) +
                  nodeArr(i + 1, j, k + 1, iy_) - nodeArr(i + 1, j, k, iy_) +
                  nodeArr(i, j + 1, k + 1, iy_) - nodeArr(i, j + 1, k, iy_) +
                  nodeArr(i + 1, j + 1, k + 1, iy_) -
                  nodeArr(i + 1, j + 1, k, iy_));
          // curl - Y
          cXDZ = 0.25 * invDx[iz_] *
                 (nodeArr(i, j, k + 1, ix_) - nodeArr(i, j, k, ix_) +
                  nodeArr(i + 1, j, k + 1, ix_) - nodeArr(i + 1, j, k, ix_) +
                  nodeArr(i, j + 1, k + 1, ix_) - nodeArr(i, j + 1, k, ix_) +
                  nodeArr(i + 1, j + 1, k + 1, ix_) -
                  nodeArr(i + 1, j + 1, k, ix_));

          cZDX = 0.25 * invDx[ix_] *
                 (nodeArr(i + 1, j, k, iz_) - nodeArr(i, j, k, iz_) +
                  nodeArr(i + 1, j, k + 1, iz_) - nodeArr(i, j, k + 1, iz_) +
                  nodeArr(i + 1, j + 1, k, iz_) - nodeArr(i, j + 1, k, iz_) +
                  nodeArr(i + 1, j + 1, k + 1, iz_) -
                  nodeArr(i, j + 1, k + 1, iz_));

          // curl - Z
          cYDX = 0.25 * invDx[ix_] *
                 (nodeArr(i + 1, j, k, iy_) - nodeArr(i, j, k, iy_) +
                  nodeArr(i + 1, j, k + 1, iy_) - nodeArr(i, j, k + 1, iy_) +
                  nodeArr(i + 1, j + 1, k, iy_) - nodeArr(i, j + 1, k, iy_) +
                  nodeArr(i + 1, j + 1, k + 1, iy_) -
                  nodeArr(i, j + 1, k + 1, iy_));

          cXDY = 0.25 * invDx[iy_] *
                 (nodeArr(i, j + 1, k, ix_) - nodeArr(i, j, k, ix_) +
                  nodeArr(i, j + 1, k + 1, ix_) - nodeArr(i, j, k + 1, ix_) +
                  nodeArr(i + 1, j + 1, k, ix_) - nodeArr(i + 1, j, k, ix_) +
                  nodeArr(i + 1, j + 1, k + 1, ix_) -
                  nodeArr(i + 1, j, k + 1, ix_));

          centerArr(i, j, k, ix_) = cZDY - cYDZ;
          centerArr(i, j, k, iy_) = cXDZ - cZDX;
          centerArr(i, j, k, iz_) = cYDX - cXDY;
        }
  }
}

void average_center_to_node(const amrex::MultiFab& centerMF,
                            amrex::MultiFab& nodeMF) {
  for (MFIter mfi(nodeMF, doTiling); mfi.isValid(); ++mfi) {
    const Box& box = mfi.fabbox();
    const Array4<Real>& nodeArr = nodeMF[mfi].array();
    const Array4<Real const>& centerArr = centerMF[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int iVar = 0; iVar < centerMF.nComp(); iVar++)
      for (int i = lo.x + 1; i <= hi.x - 1; ++i)
        for (int j = lo.y + 1; j <= hi.y - 1; ++j)
          for (int k = lo.z + 1; k <= hi.z - 1; ++k) {
            nodeArr(i, j, k, iVar) =
                0.125 *
                (centerArr(i - 1, j - 1, k - 1, iVar) +
                 centerArr(i - 1, j - 1, k, iVar) +
                 centerArr(i - 1, j, k - 1, iVar) +
                 centerArr(i - 1, j, k, iVar) +
                 centerArr(i, j - 1, k - 1, iVar) +
                 centerArr(i, j - 1, k, iVar) + centerArr(i, j, k - 1, iVar) +
                 centerArr(i, j, k, iVar));
          }
  }
}


double read_mem_usage(){                                                                                                    
  // This function returns the resident set size (RSS) of                                                                   
  // this processor in unit MB.                                                                                             
                                                                                                                            
  // From wiki:                                                                                                             
  // RSS is the portion of memory occupied by a process that is                                                             
  // held in main memory (RAM).                                                                                             
                                                                                                                            
  double rssMB = 0.0;                                                                                                       
                                                                                                                            
  ifstream stat_stream("/proc/self/stat", ios_base::in);                                                                    
                                                                                                                            
  if (!stat_stream.fail()) {                                                                                                
    // Dummy vars for leading entries in stat that we don't care about                                                      
    string pid, comm, state, ppid, pgrp, session, tty_nr;                                                                   
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;                                                                  
    string utime, stime, cutime, cstime, priority, nice;                                                                    
    string O, itrealvalue, starttime;                                                                                       
                                                                                                                            
    // Two values we want                                                                                                   
    unsigned long vsize;                                                                                                    
    unsigned long rss;                                                                                                      
                                                                                                                            
    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >>                                             
      tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime >>                                                  
      stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >>                                                  
      starttime >> vsize >> rss; // Ignore the rest                                                                         
    stat_stream.close();                                                                                                    
                                                                                                                            
    rssMB = rss*sysconf(_SC_PAGE_SIZE)/1024.0/1024.0;                                                                       
  }                                                                                                                         
                                                                                                                            
  return rssMB;                                                                                                             
}  
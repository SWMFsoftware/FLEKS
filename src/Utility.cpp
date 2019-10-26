#include "Utility.h"
#include "Constants.h"

using namespace amrex;

void lap_node_to_node(const amrex::MultiFab& srcMF, amrex::MultiFab& dstMF,
                      const amrex::DistributionMapping dm,
                      const amrex::Geometry& geom) {
  const amrex::Real* invDx = geom.InvCellSize();

  BoxArray centerBA =
      convert(srcMF.boxArray(), IntVect{ AMREX_D_DECL(0, 0, 0) });

  // Need and just need 1 ghost cell layer.
  MultiFab centerMF(centerBA, dm, 3, 1);

  for (int i = 0; i < srcMF.nComp(); i++) {
    MultiFab srcAliasMF(srcMF, amrex::make_alias, i, 1);
    grad_node_to_center(srcAliasMF, centerMF, invDx);

    centerMF.FillBoundary(geom.periodicity());

    MultiFab dstAliasMF(dstMF, amrex::make_alias, i, 1);
    div_center_to_node(centerMF, dstAliasMF, invDx);
  }
}

void grad_node_to_center(const amrex::MultiFab& nodeMF,
                         amrex::MultiFab& centerMF, const amrex::Real* invDx) {

  Real invdx = invDx[ix_];
  Real invdy = invDx[iy_];
  Real invdz = invDx[iz_];

  for (amrex::MFIter mfi(centerMF); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const amrex::Array4<amrex::Real>& center = centerMF[mfi].array();
    const amrex::Array4<amrex::Real const>& node = nodeMF[mfi].array();

    for (int i = lo.x; i <= hi.x; ++i)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int k = lo.z; k <= hi.z; ++k) {
          center(i, j, k, ix_) =
              .25 * (node(i + 1, j, k) - node(i, j, k)) * invdx +
              .25 * (node(i + 1, j, k + 1) - node(i, j, k + 1)) * invdx +
              .25 * (node(i + 1, j + 1, k) - node(i, j + 1, k)) * invdx +
              .25 * (node(i + 1, j + 1, k + 1) - node(i, j + 1, k + 1)) * invdx;
          center(i, j, k, iy_) =
              .25 * (node(i, j + 1, k) - node(i, j, k)) * invdy +
              .25 * (node(i, j + 1, k + 1) - node(i, j, k + 1)) * invdy +
              .25 * (node(i + 1, j + 1, k) - node(i + 1, j, k)) * invdy +
              .25 * (node(i + 1, j + 1, k + 1) - node(i + 1, j, k + 1)) * invdy;
          center(i, j, k, iz_) =
              .25 * (node(i, j, k + 1) - node(i, j, k)) * invdz +
              .25 * (node(i + 1, j, k + 1) - node(i + 1, j, k)) * invdz +
              .25 * (node(i, j + 1, k + 1) - node(i, j + 1, k)) * invdz +
              .25 * (node(i + 1, j + 1, k + 1) - node(i + 1, j + 1, k)) * invdz;
        }
  }
}

void div_center_to_node(const amrex::MultiFab& centerMF,
                        amrex::MultiFab& nodeMF, const amrex::Real* invDx) {

  Real invdx = invDx[ix_];
  Real invdy = invDx[iy_];
  Real invdz = invDx[iz_];

  for (amrex::MFIter mfi(nodeMF); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const amrex::Array4<amrex::Real const>& center = centerMF[mfi].array();
    const amrex::Array4<amrex::Real>& node = nodeMF[mfi].array();

    for (int i = lo.x; i <= hi.x; ++i)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int k = lo.z; k <= hi.z; ++k) {

          const Real compX =
              .25 * (center(i, j, k, ix_) - center(i - 1, j, k, ix_)) * invdx +
              .25 * (center(i, j, k - 1, ix_) - center(i - 1, j, k - 1, ix_)) *
                  invdx +
              .25 * (center(i, j - 1, k, ix_) - center(i - 1, j - 1, k, ix_)) *
                  invdx +
              .25 *
                  (center(i, j - 1, k - 1, ix_) -
                   center(i - 1, j - 1, k - 1, ix_)) *
                  invdx;
          const Real compY =
              .25 * (center(i, j, k, iy_) - center(i, j - 1, k, iy_)) * invdy +
              .25 * (center(i, j, k - 1, iy_) - center(i, j - 1, k - 1, iy_)) *
                  invdy +
              .25 * (center(i - 1, j, k, iy_) - center(i - 1, j - 1, k, iy_)) *
                  invdy +
              .25 *
                  (center(i - 1, j, k - 1, iy_) -
                   center(i - 1, j - 1, k - 1, iy_)) *
                  invdy;
          const Real compZ =
              .25 * (center(i, j, k, iz_) - center(i, j, k - 1, iz_)) * invdz +
              .25 * (center(i - 1, j, k, iz_) - center(i - 1, j, k - 1, iz_)) *
                  invdz +
              .25 * (center(i, j - 1, k, iz_) - center(i, j - 1, k - 1, iz_)) *
                  invdz +
              .25 *
                  (center(i - 1, j - 1, k, iz_) -
                   center(i - 1, j - 1, k - 1, iz_)) *
                  invdz;
          node(i, j, k) = compX + compY + compZ;
        }
  }
}

void convert_1d_to_3d(const double* const p, amrex::MultiFab& MF,
                      amrex::Geometry& geom) {
  int iCount = 0;
  for (amrex::MFIter mfi(MF); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const amrex::Array4<amrex::Real>& arr = MF[mfi].array();

    // Avoid double counting the share edges.
    int iMax = hi.x - 1, jMax = hi.y - 1, kMax = hi.z - 1;

    if ((!geom.isPeriodic(ix_)) && box.bigEnd(ix_) == hi.x)
      iMax++;
    if ((!geom.isPeriodic(iy_)) && box.bigEnd(iy_) == hi.y)
      jMax++;
    if ((!geom.isPeriodic(iz_)) && box.bigEnd(iz_) == hi.z)
      kMax++;

    for (int i = lo.x; i <= iMax; ++i)
      for (int j = lo.y; j <= jMax; ++j)
        for (int k = lo.z; k <= kMax; ++k)
          for (int iVar = 0; iVar < MF.nComp(); iVar++) {
            arr(i, j, k, iVar) = p[iCount++];
          }
  }
}

void convert_3d_to_1d(const amrex::MultiFab& MF, double* const p,
                      amrex::Geometry& geom) {
  int iCount = 0;
  for (amrex::MFIter mfi(MF); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const amrex::Array4<amrex::Real const>& arr = MF[mfi].array();

    // Avoid double counting the share edges.
    int iMax = hi.x - 1, jMax = hi.y - 1, kMax = hi.z - 1;

    if ((!geom.isPeriodic(ix_)) && box.bigEnd(ix_) == hi.x)
      iMax++;
    if ((!geom.isPeriodic(iy_)) && box.bigEnd(iy_) == hi.y)
      jMax++;
    if ((!geom.isPeriodic(iz_)) && box.bigEnd(iz_) == hi.z)
      kMax++;

    for (int i = lo.x; i <= iMax; ++i)
      for (int j = lo.y; j <= jMax; ++j)
        for (int k = lo.z; k <= kMax; ++k)
          for (int iVar = 0; iVar < MF.nComp(); iVar++) {
            p[iCount++] = arr(i, j, k, iVar);
          }
  }
}

void print_MultiFab(amrex::MultiFab& data, std::string tag) {
  AllPrint() << "-----" << tag << " begin-----" << std::endl;
  Real sum = 0;
  for (MFIter mfi(data); mfi.isValid(); ++mfi) {
    FArrayBox& fab = data[mfi];
    const Box& box = mfi.validbox();
    Array4<Real> const& data = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int i = lo.x; i <= hi.x; ++i)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int k = lo.z; k <= hi.z; ++k)
          for (int iVar = 0; iVar < data.nComp(); iVar++) {
            AllPrint() << " i = " << i << " j = " << j << " k = " << k
                       << " iVar = " << iVar
                       << " data = " << data(i, j, k, iVar) << std::endl;
            sum += data(i, j, k, iVar);
          }
  }
  AllPrint()<<"sum = "<<sum<<" on proc = "<<ParallelDescriptor::MyProc()<<std::endl;
  AllPrint() << "-----" << tag << " end-----" << std::endl;
}

void curl_center_to_node(const MultiFab& centerMF, MultiFab& nodeMF,
                         const Real* invDx) {
  Real compZDY, compYDZ;

  Real compXDZ, compZDX;
  Real compYDX, compXDY;

  for (MFIter mfi(nodeMF); mfi.isValid(); ++mfi) // Loop over grids
  {
    const Box& box = mfi.validbox();
    const Array4<Real>& nodeArr = nodeMF[mfi].array();
    const Array4<Real const>& centerArr = centerMF[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          // Needs to be improved. --Yuxi
          compZDY =
              .25 * (centerArr(i, j, k, iz_) - centerArr(i, j - 1, k, iz_)) *
                  invDx[iy_] +
              .25 *
                  (centerArr(i, j, k - 1, iz_) -
                   centerArr(i, j - 1, k - 1, iz_)) *
                  invDx[iy_] +
              .25 *
                  (centerArr(i - 1, j, k, iz_) -
                   centerArr(i - 1, j - 1, k, iz_)) *
                  invDx[iy_] +
              .25 *
                  (centerArr(i - 1, j, k - 1, iz_) -
                   centerArr(i - 1, j - 1, k - 1, iz_)) *
                  invDx[iy_];
          compYDZ =
              .25 * (centerArr(i, j, k, iy_) - centerArr(i, j, k - 1, iy_)) *
                  invDx[iz_] +
              .25 *
                  (centerArr(i - 1, j, k, iy_) -
                   centerArr(i - 1, j, k - 1, iy_)) *
                  invDx[iz_] +
              .25 *
                  (centerArr(i, j - 1, k, iy_) -
                   centerArr(i, j - 1, k - 1, iy_)) *
                  invDx[iz_] +
              .25 *
                  (centerArr(i - 1, j - 1, k, iy_) -
                   centerArr(i - 1, j - 1, k - 1, iy_)) *
                  invDx[iz_];
          // curl - Y
          compXDZ =
              .25 * (centerArr(i, j, k, ix_) - centerArr(i, j, k - 1, ix_)) *
                  invDx[iz_] +
              .25 *
                  (centerArr(i - 1, j, k, ix_) -
                   centerArr(i - 1, j, k - 1, ix_)) *
                  invDx[iz_] +
              .25 *
                  (centerArr(i, j - 1, k, ix_) -
                   centerArr(i, j - 1, k - 1, ix_)) *
                  invDx[iz_] +
              .25 *
                  (centerArr(i - 1, j - 1, k, ix_) -
                   centerArr(i - 1, j - 1, k - 1, ix_)) *
                  invDx[iz_];
          compZDX =
              .25 * (centerArr(i, j, k, iz_) - centerArr(i - 1, j, k, iz_)) *
                  invDx[ix_] +
              .25 *
                  (centerArr(i, j, k - 1, iz_) -
                   centerArr(i - 1, j, k - 1, iz_)) *
                  invDx[ix_] +
              .25 *
                  (centerArr(i, j - 1, k, iz_) -
                   centerArr(i - 1, j - 1, k, iz_)) *
                  invDx[ix_] +
              .25 *
                  (centerArr(i, j - 1, k - 1, iz_) -
                   centerArr(i - 1, j - 1, k - 1, iz_)) *
                  invDx[ix_];
          // curl - Z
          compYDX =
              .25 * (centerArr(i, j, k, iy_) - centerArr(i - 1, j, k, iy_)) *
                  invDx[ix_] +
              .25 *
                  (centerArr(i, j, k - 1, iy_) -
                   centerArr(i - 1, j, k - 1, iy_)) *
                  invDx[ix_] +
              .25 *
                  (centerArr(i, j - 1, k, iy_) -
                   centerArr(i - 1, j - 1, k, iy_)) *
                  invDx[ix_] +
              .25 *
                  (centerArr(i, j - 1, k - 1, iy_) -
                   centerArr(i - 1, j - 1, k - 1, iy_)) *
                  invDx[ix_];
          compXDY =
              .25 * (centerArr(i, j, k, ix_) - centerArr(i, j - 1, k, ix_)) *
                  invDx[iy_] +
              .25 *
                  (centerArr(i, j, k - 1, ix_) -
                   centerArr(i, j - 1, k - 1, ix_)) *
                  invDx[iy_] +
              .25 *
                  (centerArr(i - 1, j, k, ix_) -
                   centerArr(i - 1, j - 1, k, ix_)) *
                  invDx[iy_] +
              .25 *
                  (centerArr(i - 1, j, k - 1, ix_) -
                   centerArr(i - 1, j - 1, k - 1, ix_)) *
                  invDx[iy_];

          nodeArr(i, j, k, ix_) = compZDY - compYDZ;
          nodeArr(i, j, k, iy_) = compXDZ - compZDX;
          nodeArr(i, j, k, iz_) = compYDX - compXDY;

          Print() << " i = " << i << " j = " << j << " k = " << k
                  << " nodeX = " << nodeArr(i, j, k, ix_)
                  << " nodeY = " << nodeArr(i, j, k, iy_)
                  << " nodeZ = " << nodeArr(i, j, k, iz_) << std::endl;
        }
  }
}

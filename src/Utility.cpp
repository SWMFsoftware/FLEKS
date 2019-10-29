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

  for (amrex::MFIter mfi(centerMF, doTiling); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.tilebox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const amrex::Array4<amrex::Real>& center = centerMF[mfi].array();
    const amrex::Array4<amrex::Real const>& node = nodeMF[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
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
    const amrex::Box& box = mfi.tilebox();
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
    const amrex::Box& box = mfi.tilebox();
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
    const amrex::Box& box = mfi.tilebox();
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

void convert_1d_to_3d(const double* const p, amrex::MultiFab& MF,
                      amrex::Geometry& geom) {
  int iCount = 0;
  for (amrex::MFIter mfi(MF, doTiling); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.tilebox();
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

    for (int iVar = 0; iVar < MF.nComp(); iVar++)
      for (int k = lo.z; k <= kMax; ++k)
        for (int j = lo.y; j <= jMax; ++j)
          for (int i = lo.x; i <= iMax; ++i) {
            arr(i, j, k, iVar) = p[iCount++];
          }
  }
}

void convert_3d_to_1d(const amrex::MultiFab& MF, double* const p,
                      amrex::Geometry& geom) {
  int iCount = 0;
  for (amrex::MFIter mfi(MF, doTiling); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.tilebox();
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

    for (int iVar = 0; iVar < MF.nComp(); iVar++)
      for (int k = lo.z; k <= kMax; ++k)
        for (int j = lo.y; j <= jMax; ++j)
          for (int i = lo.x; i <= iMax; ++i) {
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
  AllPrint() << "sum = " << sum << " on proc = " << ParallelDescriptor::MyProc()
             << std::endl;
  AllPrint() << "-----" << tag << " end-----" << std::endl;
}

void curl_center_to_node(const MultiFab& centerMF, MultiFab& nodeMF,
                         const Real* invDx) {
  Real cZDY, cYDZ;

  Real cXDZ, cZDX;
  Real cYDX, cXDY;

  for (MFIter mfi(nodeMF, doTiling); mfi.isValid(); ++mfi) {
    const Box& box = mfi.tilebox();
    const Array4<Real>& nodeArr = nodeMF[mfi].array();
    const Array4<Real const>& centerArr = centerMF[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
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
    const Box& box = mfi.tilebox();
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
    const Box& box = mfi.tilebox();
    const Array4<Real>& nodeArr = nodeMF[mfi].array();
    const Array4<Real const>& centerArr = centerMF[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int iVar = 0; iVar < centerMF.nComp(); iVar++)
      for (int i = lo.x; i <= hi.x; ++i)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int k = lo.z; k <= hi.z; ++k) {
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

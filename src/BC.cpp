
#include <AMReX_BC_TYPES.H>
#include <AMReX_PhysBCFunct.H>
#include <Constants.h>

#include "BC.h"

using namespace amrex;

// // Use FabArray.setDomainBoundary()!!!!!!!!!!!!!!!!!!!
// void zero_boundary_cpu(amrex::Box const& bx,
//                        amrex::Array4<amrex::Real> const& arr, const int
//                        iStart, const int nComp, amrex::GeometryData const&
//                        geom, const amrex::Real time, const amrex::BCRec* bcr,
//                        const int bcomp, const int orig_comp) {
//   //"g" means "global".
//   const Box& gbx = geom.Domain();
//   const auto glo = gbx.loVect(); // Do not include ghost cells.
//   const auto ghi = gbx.hiVect();

//   int igMin = glo[ix_], igMax = ghi[ix_];
//   int jgMin = glo[iy_], jgMax = ghi[iy_];
//   int kgMin = glo[iz_], kgMax = ghi[iz_];

//   // Include ghost cells.
//   const auto lo = bx.loVect();
//   const auto hi = bx.hiVect();

//   int iMin = lo[ix_], iMax = hi[ix_];
//   int jMin = lo[iy_], jMax = hi[iy_];
//   int kMin = lo[iz_], kMax = hi[iz_];

//   const Real val = 0;

//   // x left
//   if (bcr->lo(ix_) == BCType::ext_dir) {
//     for (int iVar = iStart; iVar < nComp; iVar++)
//       for (int k = kMin; k <= kMax; k++)
//         for (int j = jMin; j <= jMax; j++)
//           for (int i = iMin; i <= igMin - 1 + nVirGst; i++) {
//             arr(i, j, k, iVar) = val;
//           }
//   }

//   // x right
//   if (bcr->hi(ix_) == BCType::ext_dir) {
//     for (int iVar = iStart; iVar < nComp; iVar++)
//       for (int k = kMin; k <= kMax; k++)
//         for (int j = jMin; j <= jMax; j++)
//           for (int i = igMax + 1 - nVirGst; i <= iMax; i++) {
//             arr(i, j, k, iVar) = val;
//           }
//   }

//   // y left
//   if (bcr->lo(iy_) == BCType::ext_dir) {
//     for (int iVar = iStart; iVar < nComp; iVar++)
//       for (int k = kMin; k <= kMax; k++)
//         for (int j = jMin; j <= jgMin - 1 + nVirGst; j++)
//           for (int i = iMin; i <= iMax; i++) {
//             arr(i, j, k, iVar) = val;
//           }
//   }

//   // y right
//   if (bcr->lo(iy_) == BCType::ext_dir) {
//     for (int iVar = iStart; iVar < nComp; iVar++)
//       for (int k = kMin; k <= kMax; k++)
//         for (int j = jgMax + 1 - nVirGst; j <= jMax; j++)
//           for (int i = iMin; i <= iMax; i++) {
//             arr(i, j, k, iVar) = val;
//           }
//   }

//   // z left
//   if (bcr->lo(iz_) == BCType::ext_dir) {
//     for (int iVar = iStart; iVar < nComp; iVar++)
//       for (int k = kMin; k <= kgMin - 1 + nVirGst; k++)
//         for (int j = jMin; j <= jMax; j++)
//           for (int i = iMin; i <= iMax; i++) {
//             arr(i, j, k, iVar) = val;
//           }
//   }

//   // z right
//   if (bcr->hi(iz_) == BCType::ext_dir) {
//     for (int iVar = iStart; iVar < nComp; iVar++)
//       for (int k = kgMax + 1 - nVirGst; k <= kMax; k++)
//         for (int j = jMin; j <= jMax; j++)
//           for (int i = iMin; i <= iMax; i++) {
//             arr(i, j, k, iVar) = val;
//           }
//   }
// }

void apply_float_boundary(MultiFab& mf, const Geometry& geom, const int iStart,
                          const int nComp, const int nshift) {

  if (geom.isAllPeriodic())
    return;
  if (mf.nGrow() == 0)
    return;

  //! create a grown domain box containing valid + periodic cells
  const Box& domain = geom.Domain();
  Box gdomain = amrex::convert(domain, mf.boxArray().ixType());
  const IntVect& ngrow = mf.nGrowVect();
  for (int i = 0; i < nDimMax; ++i) {
    if (geom.isPeriodic(i)) {
      gdomain.grow(i, ngrow[i]);
    }
  }

  {

    Vector<BCRec> bcDomain(1);
    for (int idim = 0; idim < nDimMax; ++idim) {
      bcDomain[0].setLo(idim, BCType::foextrap);
      bcDomain[0].setHi(idim, BCType::foextrap);
    }

    //"g" means "global".
    const Box& gbx = convert(geom.Domain(), mf.boxArray().ixType());
    const auto glo = gbx.loVect(); // Do not include ghost cells.
    const auto ghi = gbx.hiVect();

    int igMin = glo[ix_], igMax = ghi[ix_];
    int jgMin = glo[iy_], jgMax = ghi[iy_];
    int kgMin = glo[iz_], kgMax = ghi[iz_];

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.fabbox();

      //! if there are cells not in the valid + periodic grown box
      //! we need to fill them here
      //!
      if (!gdomain.contains(bx)) {
        //! Based on bcDomain for the domain, we need to make bcr for this Box
        Vector<BCRec> bcr(1);
        amrex::setBC(bx, domain, 0, 0, 1, bcDomain, bcr);

        amrex::Array4<amrex::Real> const& arr = mf[mfi].array();

        // Include ghost cells.
        const auto lo = bx.loVect();
        const auto hi = bx.hiVect();

        int iMin = lo[ix_], iMax = hi[ix_];
        int jMin = lo[iy_], jMax = hi[iy_];
        int kMin = lo[iz_], kMax = hi[iz_];

        // x left
        if (bcr[0].lo(ix_) == BCType::foextrap) {
          for (int i = iMin; i <= igMin - 1 + nshift; i++)
            for (int j = jMin; j <= jMax; j++)
              for (int k = kMin; k <= kMax; k++)
                for (int iVar = iStart; iVar < nComp; iVar++)

                {
                  // Print() << "1 i = " << i << " j = " << j << " k = " << k
                  //         << " ivar = " << iVar
                  //         << " val = " << arr(i, j, k, iVar) << std::endl;

                  arr(i, j, k, iVar) = arr(igMin + nshift, j, k, iVar);

                  // Print() << "x-left i = " << i << " j = " << j << " k = " <<
                  // k
                  //         << " ivar = " << iVar
                  //         << " val = " << arr(i, j, k, iVar) << std::endl;
                }
        }

        // x right
        if (bcr[0].hi(ix_) == BCType::foextrap) {
          for (int iVar = iStart; iVar < nComp; iVar++)
            for (int k = kMin; k <= kMax; k++)
              for (int j = jMin; j <= jMax; j++)
                for (int i = igMax + 1 - nshift; i <= iMax; i++) {
                  arr(i, j, k, iVar) = arr(igMax - nshift, j, k, iVar);

                  // Print() << "x-right i = " << i << " j = " << j << " k = "
                  // << k
                  //         << " ivar = " << iVar
                  //         << " val = " << arr(i, j, k, iVar) << std::endl;
                }
        }

        // y left
        if (bcr[0].lo(iy_) == BCType::foextrap) {
          for (int iVar = iStart; iVar < nComp; iVar++)
            for (int k = kMin; k <= kMax; k++)
              for (int j = jMin; j <= jgMin - 1 + nshift; j++)
                for (int i = iMin; i <= iMax; i++) {
                  arr(i, j, k, iVar) = arr(i, jgMin + nshift, k, iVar);
                  // Print() << "y-left i = " << i << " j = " << j << " k = " <<
                  // k
                  //         << " ivar = " << iVar
                  //         << " val = " << arr(i, j, k, iVar) << std::endl;
                }
        }

        // y right
        if (bcr[0].hi(iy_) == BCType::foextrap) {
          for (int iVar = iStart; iVar < nComp; iVar++)
            for (int k = kMin; k <= kMax; k++)
              for (int j = jgMax + 1 - nshift; j <= jMax; j++)
                for (int i = iMin; i <= iMax; i++) {
                  arr(i, j, k, iVar) = arr(i, jgMax - nshift, k, iVar);

                  // Print() << "y-right i = " << i << " j = " << j << " k = "
                  // << k
                  //         << " ivar = " << iVar
                  //         << " val = " << arr(i, j, k, iVar) << std::endl;
                }
        }

        // z left
        if (bcr[0].lo(iz_) == BCType::foextrap) {
          for (int iVar = iStart; iVar < nComp; iVar++)
            for (int k = kMin; k <= kgMin - 1 + nshift; k++)
              for (int j = jMin; j <= jMax; j++)
                for (int i = iMin; i <= iMax; i++) {
                  arr(i, j, k, iVar) = arr(i, j, kgMin + nshift, iVar);
                  // Print() << "z-left i = " << i << " j = " << j << " k = " <<
                  // k
                  //         << " ivar = " << iVar
                  //         << " val = " << arr(i, j, k, iVar) << std::endl;
                }
        }

        // z right
        if (bcr[0].hi(iz_) == BCType::foextrap) {
          for (int iVar = iStart; iVar < nComp; iVar++)
            for (int k = kgMax + 1 - nshift; k <= kMax; k++)
              for (int j = jMin; j <= jMax; j++)
                for (int i = iMin; i <= iMax; i++) {
                  arr(i, j, k, iVar) = arr(i, j, kgMax - nshift, iVar);
                  // Print() << "z-right i = " << i << " j = " << j << " k = "
                  // << k
                  //         << " ivar = " << iVar
                  //         << " val = " << arr(i, j, k, iVar) << std::endl;
                }
        }
      }
    }
  }
}

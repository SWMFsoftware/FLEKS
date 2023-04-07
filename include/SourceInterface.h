#ifndef _SOURCEINTERFACE_H_
#define _SOURCEINTERFACE_H_

#include "FluidInterface.h"

// #define _MERCURY_EXOSPHERE_

#ifdef _MERCURY_EXOSPHERE_
extern "C" {
void get_source_wrapper(double xyzSI[3], double sourceSI[6]);
}
#else
inline void get_source_wrapper(double xyzSI[3], double sourceSI[6]) {}
#endif

class SourceInterface : public FluidInterface {
public:
  SourceInterface(const FluidInterface& other, int id, std::string tag,
                  FluidType typeIn = SourceFluid)
      : FluidInterface(other, id, tag, typeIn) {
    initFromSWMF = false;
  }

  // Set nodeFluid[0] from get_source_wrapper.
  void get_source_from_fluid(const FluidInterface& other) {
    std::string nameFunc = "FS:get_source_from_fluid";
    amrex::Print() << nameFunc << " is called.";

    FluidInterface::set_node_fluid(other);

    const amrex::Real* dx = Geom(0).CellSize();
    const auto plo = Geom(0).ProbLo();

    // Global NODE box.
    const amrex::Box gbx = convert(Geom(0).Domain(), { 1, 1, 1 });

    const double no2siL = get_No2SiL();

    if (!nodeFluid[0].empty())
      for (amrex::MFIter mfi(nodeFluid[0]); mfi.isValid(); ++mfi) {
        // For each block, looping through all nodes, including ghost nodes.
        const amrex::Box& box = mfi.fabbox();
        const auto lo = lbound(box);
        const auto hi = ubound(box);

        const amrex::Array4<amrex::Real>& arr = nodeFluid[0][mfi].array();

        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              int idx[nDimMax] = { i, j, k };
              for (int iDim = 0; iDim < nDimMax; iDim++) {
                if (Geom(0).isPeriodic(iDim)) {
                  idx[iDim] = shift_periodic_index(
                      idx[iDim], gbx.smallEnd(iDim), gbx.bigEnd(iDim));
                }
              }
              double xyz[3];
              for (int iDim = 0; iDim < get_fluid_dimension(); iDim++) {
                xyz[iDim] = (idx[iDim] * dx[iDim] + plo[iDim]) * no2siL;
              }

              double source[6];
              get_source_wrapper(xyz, source);

              for (int iFluid = 0; iFluid < nFluid; iFluid++) {
                arr(i, j, k, iRho_I[iFluid]) = 0;
                arr(i, j, k, iUx_I[iFluid]) = 0;
                arr(i, j, k, iUy_I[iFluid]) = 0;
                arr(i, j, k, iUz_I[iFluid]) = 0;
                arr(i, j, k, iP_I[iFluid]) = 0;
              }

              if (source[0] > 0) {
                const int iNa = 1;
                arr(i, j, k, iRho_I[iNa]) = source[0] * Si2NoRho / get_Si2NoT();
                arr(i, j, k, iUx_I[iNa]) = source[1] / source[0] * Si2NoV;
                arr(i, j, k, iUy_I[iNa]) = source[2] / source[0] * Si2NoV;
                arr(i, j, k, iUz_I[iNa]) = source[3] / source[0] * Si2NoV;
                arr(i, j, k, iP_I[iNa]) = source[4] * Si2NoP / get_Si2NoT();
                arr(i, j, k, iPe) = source[5] * Si2NoP / get_Si2NoT();
              }

              amrex::Real r = 0;
              for (int i = 0; i < 3; i++) {
                xyz[i] *= no2siL / rPlanetSi;
                r += xyz[i] * xyz[i];
              }
              r = sqrt(r);
              if (r <= 1) {
                printf("Warning: x=%e, y=%e, z=%e, r=%e < 1.0 !\n", xyz[0],
                       xyz[1], xyz[2], r);
              }
            } // for k
      }

    if (!isGridEmpty && useCurrent) {
      // The current stored in nodeFluid[0] are used to initialize particle
      // velocities. Since the source particles only contribute little to the
      // total current, set current to zero.
      amrex::MultiFab currentMF(nodeFluid[0], amrex::make_alias, iJx, nDimMax);
      currentMF.setVal(0, currentMF.nGrow());
    }
  }
};

#endif
#ifndef _FLUIDSOURCE_H_
#define _FLUIDSOURCE_H_

#include "FluidInterface.h"

//#define _MERCURY_EXOSPHERE_

#ifdef _MERCURY_EXOSPHERE_
extern "C" {
void get_source_wrapper(double xyzSI[3], double sourceSI[6]);
}
#else
inline void get_source_wrapper(double xyzSI[3], double sourceSI[6]) {}
#endif

class FluidSource : public FluidInterface {
public:
  FluidSource(const FluidInterface& other, int id, std::string tag,
              FluidType typeIn = SourceFluid)
      : FluidInterface(other, id, tag, typeIn) {
    initFromSWMF = false;
  }

  void set_node_fluid(const FluidInterface& other) override {
    std::string nameFunc = "FS:set_node_fluid";
    amrex::Print() << nameFunc << " is called.";

    FluidInterface::set_node_fluid(other);

    const amrex::Real* dx = Geom(0).CellSize();
    const auto plo = Geom(0).ProbLo();

    // Global NODE box.
    const amrex::Box gbx = convert(Geom(0).Domain(), { 1, 1, 1 });

    const double no2siL = get_No2SiL();

    if (!nodeFluid.empty())
      for (amrex::MFIter mfi(nodeFluid); mfi.isValid(); ++mfi) {
        // For each block, looping through all nodes, including ghost nodes.
        const amrex::Box& box = mfi.fabbox();
        const auto lo = lbound(box);
        const auto hi = ubound(box);

        const amrex::Array4<amrex::Real>& arr = nodeFluid[mfi].array();

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

            } // for k
      }

    {
      // The current stored in nodeFluid are used to initialize particle
      // velocities. Since the source particles only contribute little the total
      // current, set current to zero.
      if (isGridEmpty)
        return;
      if (nVarCoupling == nVarFluid)
        return;

      amrex::MultiFab currentMF(nodeFluid, amrex::make_alias, iJx, nDimMax);
      currentMF.setVal(0, currentMF.nGrow());
    }

    save_amrex_file();
  }
};

#endif
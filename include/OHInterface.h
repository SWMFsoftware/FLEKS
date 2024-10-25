#ifndef _OHINTERFACE_H_
#define _OHINTERFACE_H_

#include "FluidInterface.h"
#include "SWMFInterface.h"

class OHInterface : public FluidInterface {
public:
  OHInterface(const FluidInterface &other, int id, std::string tag,
              FluidType typeIn = SourceFluid)
      : FluidInterface(other, id, tag, typeIn) {
    initFromSWMF = false;

    if (myType != PICFluid) {
      nS = 0;
      MoMi_S.clear();
      QoQi_S.clear();
      useCurrent = false;
    } else {
      useCurrent = true;
    }

    varNames.clear();
    calc_normalization_units();
  };

  virtual int get_neu_source_region(const amrex::MFIter &mfi,
                                    const amrex::IntVect ijk,
                                    const int iLev) const override {
    if (!isnodeFluidReady)
      return -1;

    const amrex::Real gamma = 5. / 3;
    auto get_fluid_moments =
        [this, &gamma](const amrex::MFIter &mfi, const amrex::IntVect ijk,
                       const int iFluid, const int iLev, amrex::Real &n,
                       amrex::Real &u2, amrex::Real &mach2, amrex::Real &T,
                       amrex::Real &p) {
          // amu/m^3
          n = get_fluid_mass_density(mfi, ijk, iFluid, iLev) * get_No2SiRho() /
              cProtonMassSI;

          // km/s
          const amrex::Real ux =
              get_fluid_ux(mfi, ijk, iFluid, iLev) * get_No2SiV() * 1e-3;
          const amrex::Real uy =
              get_fluid_uy(mfi, ijk, iFluid, iLev) * get_No2SiV() * 1e-3;
          const amrex::Real uz =
              get_fluid_uz(mfi, ijk, iFluid, iLev) * get_No2SiV() * 1e-3;
          u2 = ux * ux + uy * uy + uz * uz;

          // Pa
          p = get_fluid_p(mfi, ijk, iFluid, iLev) * get_No2SiP();
          T = p / n / cBoltzmannSI;

          // (km/s)^2
          const amrex::Real cs2 = gamma * p / (n * cProtonMassSI) * 1e-6;

          mach2 = u2 / cs2;
        };

    const int iSW = 0, iPUI3 = 1, iPUI2 = 2, iTotal = 3, nSize = 4;

    amrex::Real n[nSize] = { 0 }, u2[nSize] = { 0 }, m2[nSize] = { 0 },
                mach2[nSize] = { 0 }, T[nSize] = { 0 }, p[nSize] = { 0 };

    for (int iFluid = 0; iFluid < nFluid; ++iFluid) {
      get_fluid_moments(mfi, ijk, iFluid, iLev, n[iFluid], u2[iFluid],
                        mach2[iFluid], T[iFluid], p[iFluid]);
      n[iTotal] += n[iFluid];
      m2[iTotal] += n[iFluid] * u2[iFluid];
      p[iTotal] += p[iFluid];
    }
    // km/s
    u2[iTotal] = m2[iTotal] / n[iTotal];
    T[iTotal] = p[iTotal] / n[iTotal] / cBoltzmannSI;

    // (km/s)^2
    const amrex::Real cs2 =
        gamma * p[iTotal] / n[iTotal] / cProtonMassSI * 1e-6;
    mach2[iTotal] = u2[iTotal] / cs2;

    const amrex::Real x = Geom(iLev).CellCenter(ijk[ix_], ix_);
    const amrex::Real y = Geom(iLev).CellCenter(ijk[iy_], iy_);
    const amrex::Real z = Geom(iLev).CellCenter(ijk[iz_], iz_);

    // cAU
    amrex::Real r = sqrt(x * x + y * y + z * z) * get_No2SiL() / cAUSI;

    int iRegion = -1;

    amrex::Real levHP = get_value(mfi, ijk, iLevSet, iLev);

    OH_get_charge_exchange_region(&iRegion, &r, &n[iTotal], &u2[iTotal],
                                  &T[iTotal], &mach2[iTotal], &levHP);

    return iRegion;
  }
};

#endif
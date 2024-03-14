#ifndef _OHINTERFACE_H_
#define _OHINTERFACE_H_

#include "FluidInterface.h"
#include "SWMFInterface.h"

class OHInterface : public FluidInterface {
public:
  OHInterface(const FluidInterface& other, int id, std::string tag,
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

    useMultiFluid = false;
    useMultiSpecies = false;
    useElectronFluid = false;
    useMhdPe = false;
    useAnisoP = false;

    // For BATSRUS MHD equations (+ LevelHP_)
    nVarFluid = 5 + 3 + 1;

    nFluid = 1;
    nSpeciesFluid = 1;
    nIon = 1;

    iRho_I.clear();
    iRhoUx_I.clear();
    iRhoUy_I.clear();
    iRhoUz_I.clear();
    iPpar_I.clear();
    iP_I.clear();
    iUx_I.clear();
    iUy_I.clear();
    iUz_I.clear();

    int idx = 0;
    iRho_I.push_back(idx++);
    varNames.push_back("rho");
    iRhoUx_I.push_back(idx);
    iUx_I.push_back(idx++);
    varNames.push_back("ux");
    iRhoUy_I.push_back(idx);
    iUy_I.push_back(idx++);
    varNames.push_back("uy");
    iRhoUz_I.push_back(idx);
    iUz_I.push_back(idx++);
    varNames.push_back("uz");
    iBx = idx++;
    varNames.push_back("bx");
    iBy = idx++;
    varNames.push_back("by");
    iBz = idx++;
    varNames.push_back("bz");
    
    iLevSet = idx++;
    varNames.push_back("levHP");

    iP_I.push_back(idx);
    iPpar_I.push_back(idx++);
    varNames.push_back("p");

    if (useCurrent) {
      iJx = idx++;
      varNames.push_back("jx");
      iJy = idx++;
      varNames.push_back("jy");
      iJz = idx++;
      varNames.push_back("jz");
    }

    calc_normalization_units();

    calc_conversion_units();
  };

  virtual int get_neu_source_region(const amrex::MFIter& mfi,
                                    const amrex::IntVect ijk, const int iFluid,
                                    const int iLev) const override {
    if (!isnodeFluidReady)
      return -1;

    // amu/m^3
    amrex::Real n = get_fluid_mass_density(mfi, ijk, iFluid, iLev) *
                    get_No2SiRho() / cProtonMassSI;

    // km/s
    const amrex::Real ux =
        get_fluid_ux(mfi, ijk, iFluid, iLev) * get_No2SiV() * 1e-3;
    const amrex::Real uy =
        get_fluid_uy(mfi, ijk, iFluid, iLev) * get_No2SiV() * 1e-3;
    const amrex::Real uz =
        get_fluid_uz(mfi, ijk, iFluid, iLev) * get_No2SiV() * 1e-3;
    amrex::Real u2 = ux * ux + uy * uy + uz * uz;

    // Pa
    const amrex::Real p = get_fluid_p(mfi, ijk, iFluid, iLev) * get_No2SiP();
    amrex::Real T = p / n / cBoltzmannSI;

    const amrex::Real gamma = 5. / 3;

    // (km/s)^2
    const amrex::Real cs2 = gamma * p / (n * cProtonMassSI) * 1e-6;

    amrex::Real mach2 = u2 / cs2;

    const amrex::Real x = Geom(iLev).CellCenter(ijk[ix_], ix_);
    const amrex::Real y = Geom(iLev).CellCenter(ijk[iy_], iy_);
    const amrex::Real z = Geom(iLev).CellCenter(ijk[iz_], iz_);

    // cAU
    amrex::Real r = sqrt(x * x + y * y + z * z) * get_No2SiL() / cAUSI;

    int iRegion = -1;

    amrex::Real levHP = get_value(mfi, ijk, iLevSet, iLev);

    OH_get_charge_exchange_region(&iRegion, &r, &n, &u2, &T, &mach2, &levHP);

    return iRegion;
  }
};

#endif
#ifndef _OHINTERFACE_H_
#define _OHINTERFACE_H_

#include "FluidInterface.h"

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
    }

    useMultiFluid = false;
    useMultiSpecies = false;
    useElectronFluid = false;
    useMhdPe = false;
    useAnisoP = false;

    // For BATSRUS MHD equations
    nVarFluid = 5 + 3;

    nVarCoupling = nVarFluid;
    if (myType == PICFluid)
      nVarCoupling += 3;

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

    iP_I.push_back(idx);
    iPpar_I.push_back(idx++);
    varNames.push_back("p");

    if (nVarCoupling > nVarFluid) {
      iJx = idx++;
      varNames.push_back("jx");
      iJy = idx++;
      varNames.push_back("jy");
      iJz = idx++;
      varNames.push_back("jz");
    }

    calc_normalized_units();
  };
};

#endif
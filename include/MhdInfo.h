#ifndef _MHDINFO_H_
#define _MHDINFO_H_

#include <AMReX_Vector.H>
#include <string>

class MhdInfo {
public:
  int nDimFluid;

  // Number of variables passing between MHD and PIC.
  int nVarFluid;

  // Number of fluid at the MHD side. One 'fluid' has its own density,
  // velocity and pressure. Electron can be one fluid.
  int nFluid;

  // Number of species at the MHD side. One 'species' only has its own density.
  int nSpeciesFluid = 0;

  // Total number of ion/electron species exist in the fluid code.
  int nIon = -1;

  // These default flags are set for stand-alone PIC initialization
  bool useMultiSpecies = false;
  bool useMultiFluid = false;
  bool useElectronFluid = true;
  bool useAnisoP = false;
  bool useMhdPe = false;

  amrex::Vector<int> iRho_I, iRhoUx_I, iRhoUy_I, iRhoUz_I, iPpar_I, iP_I, iUx_I,
      iUy_I, iUz_I;

  int iBx, iBy, iBz, iEx, iEy, iEz, iPe, iRhoTotal, iLevSet;

  // Variable names of nodeFluid.
  amrex::Vector<std::string> varNames;

  void update_oh_info() {
    useMultiFluid = nFluid > 1;
    useMultiSpecies = false;
    useElectronFluid = false;
    useMhdPe = false;
    useAnisoP = false;

    // For BATSRUS MHD equations (+ LevelHP_)
    nVarFluid = 5 * nFluid + 3 + 1;

    // TODO: should this be 0? --YC
    nSpeciesFluid = 0;
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

    if (nFluid > 1) {
      iRho_I.push_back(idx++);
      varNames.push_back("pu3rho");
      iRhoUx_I.push_back(idx);
      iUx_I.push_back(idx++);
      varNames.push_back("pu3ux");
      iRhoUy_I.push_back(idx);
      iUy_I.push_back(idx++);
      varNames.push_back("pu3uy");
      iRhoUz_I.push_back(idx);
      iUz_I.push_back(idx++);
      varNames.push_back("pu3uz");

      iP_I.push_back(idx);
      iPpar_I.push_back(idx++);
      varNames.push_back("pu3p");
    }
  }
};

#endif // _MHDINFO_H_
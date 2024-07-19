#ifndef _MHDINFO_H_
#define _MHDINFO_H_

#include <AMReX_Print.H>
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

  int iBx = -1, iBy = -1, iBz = -1;
  int iEx = -1, iEy = -1, iEz = -1;
  int iPe = -1, iRhoTotal = -1, iLevSet = -1;

  // Variable names of nodeFluid.
  amrex::Vector<std::string> varNames;

  int get_nFluid() const { return nFluid; }

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

    varNames.clear();

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
  void print_mhd_info(std::string tag = "") {
    amrex::Print() << "===============MHD Info: " << tag
                   << " ===============" << std::endl;
    amrex::Print() << "nDimFluid: " << nDimFluid << std::endl;
    amrex::Print() << "nVarFluid: " << nVarFluid << std::endl;
    amrex::Print() << "nFluid: " << nFluid << std::endl;
    amrex::Print() << "nSpeciesFluid: " << nSpeciesFluid << std::endl;
    amrex::Print() << "nIon: " << nIon << std::endl;

    amrex::Print() << "useMultiSpecies: " << (useMultiSpecies ? "T" : "F")
                   << std::endl;
    amrex::Print() << "useMultiFluid   : " << (useMultiFluid ? "T" : "F")
                   << std::endl;
    amrex::Print() << "useElectronFluid: " << (useElectronFluid ? "T" : "F")
                   << std::endl;
    amrex::Print() << "useAnisoP       : " << (useAnisoP ? "T" : "F")
                   << std::endl;
    amrex::Print() << "useMhdPe        : " << (useMhdPe ? "T" : "F")
                   << std::endl;
    amrex::Print() << "iRho_I: ";
    for (const auto& i : iRho_I) {
      amrex::Print() << i << " ";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "iRhoUx_I: ";
    for (const auto& i : iRhoUx_I) {
      amrex::Print() << i << " ";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "iRhoUy_I: ";
    for (const auto& i : iRhoUy_I) {
      amrex::Print() << i << " ";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "iRhoUz_I: ";
    for (const auto& i : iRhoUz_I) {
      amrex::Print() << i << " ";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "iPpar_I: ";
    for (const auto& i : iPpar_I) {
      amrex::Print() << i << " ";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "iP_I: ";
    for (const auto& i : iP_I) {
      amrex::Print() << i << " ";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "iUx_I: ";
    for (const auto& i : iUx_I) {
      amrex::Print() << i << " ";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "iUy_I: ";
    for (const auto& i : iUy_I) {
      amrex::Print() << i << " ";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "iUz_I: ";
    for (const auto& i : iUz_I) {
      amrex::Print() << i << " ";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "iBx: " << iBx << std::endl;
    amrex::Print() << "iBy: " << iBy << std::endl;
    amrex::Print() << "iBz: " << iBz << std::endl;
    amrex::Print() << "iEx: " << iEx << std::endl;
    amrex::Print() << "iEy: " << iEy << std::endl;
    amrex::Print() << "iEz: " << iEz << std::endl;
    amrex::Print() << "iPe: " << iPe << std::endl;
    amrex::Print() << "iRhoTotal: " << iRhoTotal << std::endl;
    amrex::Print() << "iLevSet: " << iLevSet << std::endl;
    amrex::Print() << "varNames: ";
    for (const auto& var : varNames) {
      amrex::Print() << var << " ";
    }
    amrex::Print() << std::endl;
  }
};

#endif // _MHDINFO_H_

#include <fstream>
#include <iostream>
#include <vector>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "DataContainer.h"
#include "GridUtility.h"

using namespace amrex;

void AMReXDataContainer::read_header(std::string& headerName, int& nVar,
                                     int& nDim, Real& time, int& finest_level,
                                     RealBox& domain, Box& cellBox,
                                     Vector<std::string>& varNames) {
  std::string funcName = "AMReXDataContainer::read_header()";
  BL_PROFILE(funcName);

  std::ifstream HeaderFile;
  HeaderFile.open(headerName, std::ifstream::in);

  HeaderFile.precision(17);

  std::string versionName;

  HeaderFile >> versionName;
  HeaderFile >> nVar;

  varNames.clear();
  for (int ivar = 0; ivar < nVar; ++ivar) {
    std::string var;
    HeaderFile >> var;
    varNames.push_back(var);
  }

  HeaderFile >> nDim;
  HeaderFile >> time;
  HeaderFile >> finest_level;

  for (int i = 0; i < nDim; ++i) {
    Real lo;
    HeaderFile >> lo;
    domain.setLo(i, lo);
  }

  for (int i = 0; i < nDim; ++i) {
    Real hi;
    HeaderFile >> hi;
    domain.setHi(i, hi);
  }

  for (int i = 0; i < finest_level; ++i) {
    int refRatio;
    HeaderFile >> refRatio;
  }

  if (finest_level == 0) {
    HeaderFile.ignore(100000, '\n');
  }

  // Read the base grid and ignore refined leve grids.
  HeaderFile >> cellBox;
  HeaderFile.ignore(100000, '\n');

  HeaderFile.close();
}

void AMReXDataContainer::read_header() {

  std::string headerName = dirIn + "/Header";

  int finestLev;
  AMReXDataContainer::read_header(headerName, nVar, nDim, time, finestLev,
                                  domain, cellBox, varNames);

  SetFinestLevel(finestLev);
}

void AMReXDataContainer::read() {
  std::string funcName = "AMReXDataContainer::read()";
  BL_PROFILE(funcName);

  Print() << "Reading in " << dirIn << std::endl;

  read_header();

  Grid grid(Geom(0), get_amr_info(), nGst);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    VisMF::Read(mf[iLev], dirIn + "/Level_" + std::to_string(iLev) + "/Cell");

    grid.SetBoxArray(iLev, mf[iLev].boxArray());
    grid.SetDistributionMap(iLev, mf[iLev].DistributionMap());
  }
  grid.SetFinestLevel(n_lev() - 1);

  // Print() << "grids.ba0 = " << grid.boxArray(0) << std::endl;

  regrid(grid.boxArray(0), &grid);
}

size_t AMReXDataContainer::loop_cell(bool doStore, Vector<float>& vars,
                                     bool doStoreLoc) {
  std::string funcName = "AMReXDataContainer::loop_cell()";
  BL_PROFILE(funcName);

  if (doStore)
    vars.clear();

  size_t iCount = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    const int ncomp = mf[iLev].nComp();
    for (MFIter mfi(iCell[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto& status = cell_status(iLev)[mfi].array();
      const Array4<Real>& cell = iCell[iLev][mfi].array();
      const Array4<Real>& data = mf[iLev][mfi].array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i)
            if (!bit::is_refined(status(i, j, k))) {
              iCount++;
              cell(i, j, k) = iCount;
              if (doStore) {
                if (doStoreLoc) {
                  vars.push_back(Geom(iLev).CellCenter(i, ix_));
                  vars.push_back(Geom(iLev).CellCenter(j, iy_));
                  vars.push_back(Geom(iLev).CellCenter(k, iz_));
                } else {
                  for (int iVar = 0; iVar < ncomp; iVar++)
                    vars.push_back(data(i, j, k, iVar));
                }
              }
            }
    }

    iCell[iLev].FillBoundary();
  }

  for (int iLev = n_lev() - 2; iLev >= 0; iLev--) {
    fill_fine_lev_bny_cell_from_coarse(
        iCell[iLev], iCell[iLev + 1], 0, iCell[iLev].nComp(), ref_ratio[iLev],
        Geom(iLev), Geom(iLev + 1), cell_status(iLev + 1));
  }
  return iCount;
}

size_t AMReXDataContainer::loop_zone(bool doStore, Vector<size_t>& zones) {
  std::string funcName = "AMReXDataContainer::loop_zone()";
  BL_PROFILE(funcName);

  size_t iBrick = 0;

  if (doStore)
    zones.clear();

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi(iCell[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& cell = iCell[iLev][mfi].array();
      const auto& status = cell_status(iLev)[mfi].array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      // Loop over valid cells + low end ghost cells
      for (int k = lo.z - 1; k <= hi.z; ++k)
        for (int j = lo.y - 1; j <= hi.y; ++j)
          for (int i = lo.x - 1; i <= hi.x; ++i)
            if ((iLev > 0 && bit::is_lev_boundary(status(i, j, k))) ||
                (box.contains(i, j, k) && !bit::is_refined(status(i, j, k)))) {
              bool isBrick = true;

              for (int kk = k; kk <= k + 1; kk++)
                for (int jj = j; jj <= j + 1; jj++)
                  for (int ii = i; ii <= i + 1; ii++) {
                    if (cell(ii, jj, kk) == 0)
                      isBrick = false;
                  }

              if (isBrick) {
                iBrick++;
                if (doStore) {
                  zones.push_back((int)cell(i, j, k));
                  zones.push_back((int)cell(i + 1, j, k));
                  zones.push_back((int)cell(i + 1, j + 1, k));
                  zones.push_back((int)cell(i, j + 1, k));
                  zones.push_back((int)cell(i, j, k + 1));
                  zones.push_back((int)cell(i + 1, j, k + 1));
                  zones.push_back((int)cell(i + 1, j + 1, k + 1));
                  zones.push_back((int)cell(i, j + 1, k + 1));
                }
              }
            }
    }
  }
  return iBrick;
}

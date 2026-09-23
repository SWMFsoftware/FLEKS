
#include <fstream>
#include <limits>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "DataContainer.h"
#include "GridUtility.h"

using namespace amrex;

void AMReXDataContainer::read_header(const std::string& headerName, int& nVar,
                                     int& nDim, Real& time, int& finest_level,
                                     RealBox& domain, Box& cellBox,
                                     Vector<std::string>& varNames) {
  BL_PROFILE("AMReXDataContainer::read_header");

  std::ifstream HeaderFile(headerName, std::ifstream::in);
  if (!HeaderFile.is_open()) {
    Abort("Error: cannot open header file: " + headerName);
  }

  HeaderFile.precision(17);

  std::string versionName;

  HeaderFile >> versionName;
  HeaderFile >> nVar;

  varNames.clear();
  varNames.reserve(nVar);
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
    HeaderFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  // Read the base grid and ignore refined level grids.
  HeaderFile >> cellBox;
  HeaderFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

void AMReXDataContainer::read_header() {
  const std::string headerName = filename + "/Header";

  int finestLev;
  AMReXDataContainer::read_header(headerName, nVar, nDim, time, finestLev,
                                  domain, cellBox, varNames);

  SetFinestLevel(finestLev);
}

int AMReXDataContainer::read() {
  BL_PROFILE("AMReXDataContainer::read");

  Print() << "Reading in " << filename << std::endl;

  read_header();

  nCell = 0;
  nBrick = 0;
  isCellNumbered = false;

  Grid grid(Geom(0), get_amr_info(), nGst);

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    VisMF::Read(mf[iLev],
                filename + "/Level_" + std::to_string(iLev) + "/Cell");

    grid.SetBoxArray(iLev, mf[iLev].boxArray());
    grid.SetDistributionMap(iLev, mf[iLev].DistributionMap());
  }
  grid.SetFinestLevel(n_lev() - 1);

  regrid(grid.boxArray(0), &grid);

  return iSuccess;
}

size_t AMReXDataContainer::loop_cell(bool doStore, Vector<float>& vars,
                                     bool doStoreLoc) {
  BL_PROFILE("AMReXDataContainer::loop_cell");

  if (doStore) {
    vars.clear();
    if (nCell > 0) {
      vars.reserve(nCell *
                   (doStoreLoc ? 3 : (n_lev() > 0 ? mf[0].nComp() : 1)));
    }
  }

  const bool needNumbering = !isCellNumbered;
  size_t iCount = 0;

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    const int ncomp = mf[iLev].nComp();
    const auto& geom = Geom(iLev);
    const auto dx = geom.CellSize();
    const auto probLo = geom.ProbLo();

    for (MFIter mfi(iCell[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto& status = cell_status(iLev)[mfi].array();
      const Array4<Real>& cell = iCell[iLev][mfi].array();
      const Array4<Real>& data = mf[iLev][mfi].array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k) {
#if AMREX_SPACEDIM > 2
        const float z = static_cast<float>(probLo[iz_] + (k + 0.5) * dx[iz_]);
#else
        constexpr float z = 0.0f;
#endif
        for (int j = lo.y; j <= hi.y; ++j) {
          const float y = static_cast<float>(probLo[iy_] + (j + 0.5) * dx[iy_]);
          for (int i = lo.x; i <= hi.x; ++i) {
            if (!bit::is_refined(status(i, j, k))) {
              iCount++;
              if (needNumbering) {
                cell(i, j, k) = static_cast<Real>(iCount);
              }
              if (doStore) {
                if (doStoreLoc) {
                  const float x =
                      static_cast<float>(probLo[ix_] + (i + 0.5) * dx[ix_]);
                  vars.push_back(x);
                  vars.push_back(y);
                  vars.push_back(z);
                } else {
                  for (int iVar = 0; iVar < ncomp; iVar++)
                    vars.push_back(static_cast<float>(data(i, j, k, iVar)));
                }
              }
            }
          }
        }
      }
    }

    if (needNumbering) {
      iCell[iLev].FillBoundary();
    }
  }

  if (needNumbering) {
    for (int iLev = n_lev() - 2; iLev >= 0; --iLev) {
      fill_fine_lev_bny_from_coarse(
          iCell[iLev], iCell[iLev + 1], 0, iCell[iLev].nComp(), ref_ratio[iLev],
          Geom(iLev), Geom(iLev + 1), cell_status(iLev + 1), pc_interp);
    }
    isCellNumbered = true;
  }

  nCell = iCount;
  return iCount;
}

size_t AMReXDataContainer::loop_zone(bool doStore, Vector<size_t>& zones) {
  BL_PROFILE("AMReXDataContainer::loop_zone");

  size_t iBrick = 0;

  if (doStore) {
    zones.clear();
    if (nBrick > 0) {
      zones.reserve(nBrick * 8);
    }
  }

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    for (MFIter mfi(iCell[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& cell = iCell[iLev][mfi].array();
      const auto& status = cell_status(iLev)[mfi].array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      // Loop over valid cells + lower end ghost cells
      for (int k = lo.z - 1; k <= hi.z; ++k) {
        for (int j = lo.y - 1; j <= hi.y; ++j) {
          for (int i = lo.x - 1; i <= hi.x; ++i) {
            const bool inValidBox = (i >= lo.x && j >= lo.y && k >= lo.z);
            if ((iLev > 0 && bit::is_lev_boundary(status(i, j, k))) ||
                (inValidBox && !bit::is_refined(status(i, j, k)))) {
              bool isBrick = true;

              for (int kk = k; kk <= k + 1 && isBrick; ++kk) {
                for (int jj = j; jj <= j + 1 && isBrick; ++jj) {
                  for (int ii = i; ii <= i + 1; ++ii) {
                    if (cell(ii, jj, kk) == 0) {
                      isBrick = false;
                      break;
                    }
                  }
                }
              }

              if (isBrick) {
                iBrick++;
                if (doStore) {
                  zones.push_back(static_cast<size_t>(cell(i, j, k)));
                  zones.push_back(static_cast<size_t>(cell(i + 1, j, k)));
                  zones.push_back(static_cast<size_t>(cell(i + 1, j + 1, k)));
                  zones.push_back(static_cast<size_t>(cell(i, j + 1, k)));
                  zones.push_back(static_cast<size_t>(cell(i, j, k + 1)));
                  zones.push_back(static_cast<size_t>(cell(i + 1, j, k + 1)));
                  zones.push_back(
                      static_cast<size_t>(cell(i + 1, j + 1, k + 1)));
                  zones.push_back(static_cast<size_t>(cell(i, j + 1, k + 1)));
                }
              }
            }
          }
        }
      }
    }
  }

  nBrick = iBrick;
  return iBrick;
}

void AMReXDataContainer::smooth(int nSmooth) {
  BL_PROFILE("AMReXDataContainer::smooth");

  constexpr Real coef = 0.5;
  constexpr Real weightSelf = 1.0 - coef;
  constexpr Real weightNei = coef / 2.0;

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    const int ng = 1;
    MultiFab mfOld(mf[iLev].boxArray(), mf[iLev].DistributionMap(),
                   mf[iLev].nComp(), ng);

    for (int iSmooth = 0; iSmooth < nSmooth; iSmooth++) {
      auto smooth_dir = [&](int iDir) {
        MultiFab::Copy(mfOld, mf[iLev], 0, 0, mf[iLev].nComp(), 0);
        mfOld.FillBoundary(Geom(iLev).periodicity());

        int dIdx[3] = { 0, 0, 0 };
        dIdx[iDir] = 1;
        const int di = dIdx[ix_];
        const int dj = dIdx[iy_];
        const int dk = dIdx[iz_];

        const int ncomp = mf[iLev].nComp();

        for (MFIter mfi(mf[iLev]); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();

          const auto& status = cell_status(iLev)[mfi].array();
          Array4<Real> const& arr = mf[iLev][mfi].array();
          Array4<Real> const& tmp = mfOld[mfi].array();

          const auto lo = lbound(bx);
          const auto hi = ubound(bx);

          for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
              for (int i = lo.x; i <= hi.x; ++i) {
                if (bit::is_lev_edge(status(i, j, k))) {
                  continue;
                }

                auto is_bnd_or_refined = [&](int ci, int cj, int ck) {
                  const auto st = status(ci, cj, ck);
                  return bit::is_lev_boundary(st) || bit::is_refined(st);
                };

                if (is_bnd_or_refined(i - di, j - dj, k - dk) ||
                    is_bnd_or_refined(i + di, j + dj, k + dk)) {
                  continue;
                }

                for (int iVar = 0; iVar < ncomp; ++iVar) {
                  const Real neiSum = tmp(i - di, j - dj, k - dk, iVar) +
                                      tmp(i + di, j + dj, k + dk, iVar);
                  arr(i, j, k, iVar) =
                      weightSelf * arr(i, j, k, iVar) + weightNei * neiSum;
                }
              }
            }
          }
        }

        mf[iLev].FillBoundary(Geom(iLev).periodicity());
      };

      smooth_dir(ix_);
      smooth_dir(iy_);
      smooth_dir(iz_);
    }
  }
}

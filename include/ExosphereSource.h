#ifndef _EXOSPHERESOURCE_H_
#define _EXOSPHERESOURCE_H_

#include "SourceInterface.h"
#include "Particles.h"

class ExosphereSource : public SourceInterface {
private:
  std::vector<ExosphereInfo> exoParams;
  TestCase testCase = RegularSimulation;
  double pickup_xMin = -1.0;
  double pickup_xMax = 1.0;

public:
  ExosphereSource(const FluidInterface& other, int id, std::string tag,
                  FluidType typeIn = SourceFluid)
      : SourceInterface(other, id, tag, typeIn) {
    info = "Exosphere Photoionization Source Class";
  }

  void add_exosphere_params(const std::vector<ExosphereInfo>& infos) {
    exoParams = infos;
  }

  void set_pickup_params(TestCase tCase, double xMin, double xMax) {
    testCase = tCase;
    pickup_xMin = xMin;
    pickup_xMax = xMax;
  }

  void post_regrid() override {
    distribute_arrays();
    set_source_standalone();
  }

  void set_source(const FluidInterface& other) override {
    set_node_fluid(other);
    set_source_standalone();
  }

  void set_source_standalone() {
    amrex::Print() << "ExosphereSource::set_source_standalone() is called!" << std::endl;

    // 1. Find electron species if any
    int iElec = -1;
    for (int i = 0; i < nFluid; ++i) {
      if (QoQi_S[i] < 0.0) {
        iElec = i;
        break;
      }
    }

    // 2. Initialize all nodeFluid rates to 0
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      if (nodeFluid[iLev].empty()) continue;
      nodeFluid[iLev].setVal(0.0);
    }

    // 3. Write un-normalized analytical density and pressure shapes to nodes
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      if (nodeFluid[iLev].empty()) continue;

      const amrex::Real* dx = Geom(iLev).CellSize();
      const auto plo = Geom(iLev).ProbLo();
      const amrex::Box gbx = convert(Geom(0).Domain(), { AMREX_D_DECL(1, 1, 1) });

      for (amrex::MFIter mfi(nodeFluid[iLev]); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.fabbox();
        const auto lo = lbound(box);
        const auto hi = ubound(box);
        const amrex::Array4<amrex::Real>& arr = nodeFluid[iLev][mfi].array();

        for (int k = lo.z; k <= hi.z; ++k) {
          for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
              
              amrex::IntVect idx = { AMREX_D_DECL(i, j, k) };
              for (int iDim = 0; iDim < nDim; iDim++) {
                if (Geom(iLev).isPeriodic(iDim)) {
                  idx[iDim] = shift_periodic_index(
                      idx[iDim], gbx.smallEnd(iDim), gbx.bigEnd(iDim));
                }
              }

              // Calculate node position in PIC units
              double xyz[3] = { 0, 0, 0 };
              for (int iDim = 0; iDim < get_fluid_dimension(); iDim++) {
                xyz[iDim] = (idx[iDim] * dx[iDim] + plo[iDim]);
              }

              // Compute contribution from each species profile
              for (size_t iInfo = 0; iInfo < exoParams.size(); ++iInfo) {
                const auto& param = exoParams[iInfo];
                int iSp = param.iSpecies - 1;
                if (iSp < 0 || iSp >= nFluid) continue;

                double r = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
                if (r < param.exobaseRadius) continue;

                if (testCase == Pickup && (xyz[0] < pickup_xMin || xyz[0] > pickup_xMax)) continue;

                if (param.shadowRadius > 0 && xyz[0] < 0) {
                  double perp2 = xyz[1]*xyz[1] + (get_fluid_dimension() > 2 ? xyz[2]*xyz[2] : 0.0);
                  if (perp2 < param.shadowRadius * param.shadowRadius) continue;
                }

                double dens = 0.0;
                int n0_size = param.n0.size();
                if (param.neutralProfile == "exponential") {
                  for (int idx_p = 0; idx_p < n0_size; idx_p++) {
                    dens += param.n0[idx_p] * exp(-(r - param.r0) / param.H0[idx_p]);
                  }
                } else if (param.neutralProfile == "power-law" || param.neutralProfile == "PowerLaw") {
                  for (int idx_p = 0; idx_p < n0_size; idx_p++) {
                    dens += param.n0[idx_p] * pow(param.r0 / r, param.k0[idx_p]);
                  }
                } else if (param.neutralProfile == "ChamberlainH") {
                  for (int idx_p = 0; idx_p < n0_size; idx_p++) {
                    dens += param.n0[idx_p] * exp(-param.H0[idx_p] * (1.0 / param.r0 - 1.0 / r));
                  }
                }

                arr(i, j, k, iRho_I[iSp]) += dens;
                
                double T0_K = param.T0.empty() ? 0.0 : param.T0[0];
                double mass_kg = MoMi_S[iSp] * get_No2SiM();
                double uth_SI = (T0_K > 0 && mass_kg > 0) ? sqrt(cBoltzmannSI * T0_K / mass_kg) : 0.0;
                double uth = uth_SI * get_Si2NoV();
                
                arr(i, j, k, iP_I[iSp]) += dens * uth * uth;
                if (useMhdPe) {
                  arr(i, j, k, iPe) += dens * uth * uth * 1e-2;
                }
              }
            }
          }
        }
      }
    }

    // 4. Compute exact cell-center interpolated global volume integrals
    std::vector<double> sumGlobal(exoParams.size(), 0.0);
    for (size_t iInfo = 0; iInfo < exoParams.size(); ++iInfo) {
      int iSp = exoParams[iInfo].iSpecies - 1;
      if (iSp < 0 || iSp >= nFluid) continue;

      double sumLocal = 0.0;
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        if (nodeFluid[iLev].empty()) continue;

        const amrex::Real* dx = Geom(iLev).CellSize();
        const double vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

        for (amrex::MFIter mfi(nodeFluid[iLev]); mfi.isValid(); ++mfi) {
          amrex::Box cell_box = mfi.validbox();
          cell_box.convert(amrex::IndexType::TheCellType());

          const auto lo = lbound(cell_box);
          const auto hi = ubound(cell_box);
          const amrex::Array4<amrex::Real>& arr = nodeFluid[iLev][mfi].array();

          for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
              for (int i = lo.x; i <= hi.x; ++i) {
                // Direct trilinear/bilinear node average to cell center to match evaluation exactly
                double sum_node = 0.0;
                int i_max = i + 1;
                int j_max = (get_fluid_dimension() > 1) ? j + 1 : j;
                int k_max = (get_fluid_dimension() > 2) ? k + 1 : k;

                for (int nk = k; nk <= k_max; nk++) {
                  for (int nj = j; nj <= j_max; nj++) {
                    for (int ni = i; ni <= i_max; ni++) {
                      sum_node += arr(ni, nj, nk, iRho_I[iSp]);
                    }
                  }
                }
                
                int n_corners = (i_max - i + 1) * (j_max - j + 1) * (k_max - k + 1);
                double dens = sum_node / n_corners;
                sumLocal += dens * vol;
              }
            }
          }
        }
      }

      sumGlobal[iInfo] = sumLocal;
      amrex::ParallelDescriptor::ReduceRealSum(sumGlobal[iInfo]);
    }

    // 5. Scale components by norm_factor
    for (size_t iInfo = 0; iInfo < exoParams.size(); ++iInfo) {
      int iSp = exoParams[iInfo].iSpecies - 1;
      if (iSp < 0 || iSp >= nFluid) continue;

      if (sumGlobal[iInfo] > 0.0) {
        double norm_factor = (exoParams[iInfo].totalProductionRate / sumGlobal[iInfo])
                             * (MoMi_S[iSp] / get_Si2NoT());

        amrex::Print() << "  Species Info " << iInfo 
                       << ": sumGlobal = " << sumGlobal[iInfo] 
                       << ", norm_factor = " << norm_factor << std::endl;

        for (int iLev = 0; iLev < n_lev(); iLev++) {
          if (nodeFluid[iLev].empty()) continue;
          nodeFluid[iLev].mult(norm_factor, iRho_I[iSp], 1, 0);
          nodeFluid[iLev].mult(norm_factor, iP_I[iSp], 1, 0);
          if (useMhdPe) {
            nodeFluid[iLev].mult(norm_factor, iPe, 1, 0);
          }
        }
      }
    }

    // 6. Compute neutralizing electron rates for charge neutrality
    if (iElec >= 0) {
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        if (nodeFluid[iLev].empty()) continue;

        for (amrex::MFIter mfi(nodeFluid[iLev]); mfi.isValid(); ++mfi) {
          const amrex::Box& box = mfi.fabbox();
          const auto lo = amrex::lbound(box);
          const auto hi = amrex::ubound(box);
          const amrex::Array4<amrex::Real>& arr = nodeFluid[iLev][mfi].array();

          for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
              for (int i = lo.x; i <= hi.x; ++i) {
                double sum_num_dens = 0.0;
                for (int iSp = 0; iSp < nFluid; iSp++) {
                  if (iSp != iElec) {
                    sum_num_dens += arr(i, j, k, iRho_I[iSp]) / MoMi_S[iSp];
                  }
                }
                arr(i, j, k, iRho_I[iElec]) = sum_num_dens * MoMi_S[iElec];

                double sum_p = 0.0;
                for (int iSp = 0; iSp < nFluid; iSp++) {
                  if (iSp != iElec) {
                    double T0_K = 0.0;
                    for (size_t iInfo = 0; iInfo < exoParams.size(); ++iInfo) {
                      if (exoParams[iInfo].iSpecies - 1 == iSp) {
                        T0_K = exoParams[iInfo].T0.empty() ? 0.0 : exoParams[iInfo].T0[0];
                        break;
                      }
                    }
                    double mass_kg = MoMi_S[iElec] * get_No2SiM();
                    double uth_SI = (T0_K > 0 && mass_kg > 0) ? sqrt(cBoltzmannSI * T0_K / mass_kg) : 0.0;
                    double uth = uth_SI * get_Si2NoV();
                    sum_p += (arr(i, j, k, iRho_I[iSp]) / MoMi_S[iSp]) * MoMi_S[iElec] * uth * uth;
                  }
                }
                arr(i, j, k, iP_I[iElec]) = sum_p;
              }
            }
          }
        }
      }
    }

    if (!isGridEmpty && useCurrent) {
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        amrex::MultiFab currentMF(nodeFluid[iLev], amrex::make_alias, iJx, 3);
        currentMF.setVal(0, currentMF.nGrow());
      }
    }

    for (int iSp = 0; iSp < nFluid; iSp++) {
      double sumVal = 0.0;
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        if (!nodeFluid[iLev].empty()) {
          sumVal += nodeFluid[iLev].sum(iRho_I[iSp]);
        }
      }
      amrex::Print() << "  Species " << iSp << " nodeFluid mass rate sum = " << sumVal << std::endl;
    }
  }
};

#endif

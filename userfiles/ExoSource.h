#ifndef _EXOSOURCE_H_
#define _EXOSOURCE_H_

#include "SourceInterface.h"

class UserSource : public SourceInterface {
public:
  UserSource(Grid& gridIn, const FluidInterface& other, int id, std::string tag,
             FluidType typeIn)
      : SourceInterface(gridIn, other, id, tag, typeIn) {
    info = "Exosphere Source";
    useFluidSource = true;
  }

  void register_parameter_commands(
      std::vector<std::string>& commands) const override {
    commands.push_back("#PHOTOIONIZATION");
    commands.push_back("#ELECTRONIMPACT");
    commands.push_back("#CHARGEEXCHANGE");
    commands.push_back("#SHADOWCYLINDER");
    commands.push_back("#OPTICALDEPTH");
    commands.push_back("#RECOMBINATION");
    commands.push_back("#CHEMISTRY");
    commands.push_back("#PHOTOCHEMICALBOUNDARY");
  }

  // ---- Species <-> nodeFluid slot mapping ----
  // The mapping depends on how FLEKS was initialized:
  //  * useElectronFluid == true (standalone runs, or an MHD equation with an
  //    electron state variable): species 0 is the electron and owns
  //    iRho_I[0], so species iSp maps to iRho_I[iSp].
  //  * useElectronFluid == false (coupled to a multi-species or multi-fluid
  //    MHD equation without an electron state variable, e.g. MhdMars):
  //    species 0 is the quasi-neutral electron and has NO slot of its own,
  //    iRho_I[0..nIon-1] hold ion species 1..nIon, and all ion species share
  //    one momentum slot (iUx_I[0]) and one pressure slot (iP_I[0]).
  // get_number_density(), get_p() and convert_moment_to_velocity() already
  // follow this convention, so every direct nodeFluid access in this file has
  // to go through the helpers below.
  bool has_electron_species() const { return (nS > 0 && QoQi_S[0] < 0.0); }

  // True when the electron has its own density/momentum/pressure slots.
  bool electron_has_state_slot() const { return useElectronFluid; }

  // True when all ion species share one momentum and one pressure slot.
  bool shared_ion_state() const { return !electron_has_state_slot(); }

  // nodeFluid slot holding the mass density of species iSp, or -1 when the
  // species has no state slot of its own (quasi-neutral electron).
  int node_rho_slot(const int iSp) const {
    if (iSp < 0 || iSp >= nS)
      return -1;
    if (iSp == 0)
      return (electron_has_state_slot() ? iRho_I[0] : -1);
    return (electron_has_state_slot() ? iRho_I[iSp] : iRho_I[iSp - 1]);
  }

  // get_p() returns (1 - PeRatio)*iP_I[0] for the ions when the electron
  // pressure is not a separate MHD state variable.  The ion pressure written
  // into the shared slot therefore has to be pre-divided by (1 - PeRatio).
  amrex::Real ion_pressure_factor() const {
    if (useMhdPe)
      return 1.0;
    return (PeRatio < 1.0 ? 1.0 / (1.0 - PeRatio) : 1.0);
  }

  // ---- EUV attenuation (#OPTICALDEPTH) ----
  // The default attenuation is the #SHADOWCYLINDER mask: 0 inside the
  // cylinder and 1 outside.  With #OPTICALDEPTH the attenuation is exp(-tau/mu)
  // with
  //   tau = sum_c n_c(r) * sigma_c * H_c(r)     (vertical column)
  //   mu  = max(proj/r, cosSzaFloor)            (cos of the solar zenith angle)
  // i.e. the model of BATSRUS ModUserMars::set_neutral_density (UseChapman = F
  // branch).  nu0 is the unattenuated rate, following the BATSRUS Rate_I
  // convention (= RateDim_I of the selected solar condition, which matches the
  // nuPhoto0 values of the Mars decks).
  bool useOpticalDepth = false;
  bool useChapmanFunction = false;
  // Exosphere component used for the Chapman slant column; < 0 selects the
  // component with the largest vertical optical depth (CO2 for Mars).
  int chapmanComponent = -1;
  amrex::Vector<amrex::Real> optCrossSection; // [m^2], per exosphere component
  amrex::Vector<amrex::Real> optScaleHeight;  // [m]; < 0 -> from the profile
  amrex::Real optSolarDir[3] = { 1.0, 0.0, 0.0 }; // unit vector toward the Sun
  amrex::Real optMinProduction = 1.0e-5;          // attenuation floor
  bool optMinProductionSet = false;    // user overrode optMinProduction
  amrex::Real optTauFloor = 6.0e-3;    // floor on the vertical optical depth
  amrex::Real optTauCutoff = 11.5;     // tau above which the floor is used
  amrex::Real optCosSzaFloor = 5.0e-4; // floor on mu
  amrex::Real optChapmanTauMax = 13.8; // Chapman branch: opaque threshold

  // Local scale height [m] of exosphere component iC at radius r, i.e. the
  // exact H = -1/(d ln n / dr) of the profile selected by #EXOSPHERE, unless
  // #OPTICALDEPTH supplied an explicit value.
  amrex::Real exosphere_scale_height(amrex::Real r, int iC) const {
    if (iC < 0 || iC >= nExoComponent)
      return 0.0;
    if (iC < static_cast<int>(optScaleHeight.size()) &&
        optScaleHeight[iC] > 0.0)
      return optScaleHeight[iC];
    if (r <= 0.0)
      return 0.0;
    const std::string& type = exoTypeC[iC];
    if (type == "Exponential")
      return (exoH0[iC] > 0.0 ? exoH0[iC] : 0.0);
    if (type == "Chamberlain")
      return (exoH0[iC] > 0.0 ? r * r / exoH0[iC] : 0.0);
    if (type == "Power-Law")
      return (exoK0[iC] > 0.0 ? r / exoK0[iC] : 0.0);
    return 0.0;
  }

  // Vertical optical depth sum_c n_c sigma_c H_c.  With onlyComp >= 0 only
  // that component contributes (the Chapman branch of BATSRUS uses CO2 only).
  amrex::Real vertical_optical_depth(amrex::Real r, int onlyComp) const {
    amrex::Real tau = 0.0;
    for (int iC = 0; iC < nExoComponent; ++iC) {
      if (onlyComp >= 0 && iC != onlyComp)
        continue;
      const amrex::Real sigma = (iC < static_cast<int>(optCrossSection.size()))
                                    ? optCrossSection[iC]
                                    : 0.0;
      if (sigma <= 0.0)
        continue;
      const amrex::Real h = exosphere_scale_height(r, iC);
      if (h <= 0.0)
        continue;
      tau += get_exosphere_component_density(r, iC) * sigma * h;
    }
    return tau;
  }

  // Optical depth along a ray that grazes the planet at the tangent radius p
  // [m]: the Chapman grazing-incidence column sum_c sigma_c n_c(p) H_c
  // sqrt(2 pi p / H_c), all of which lies between the Sun and a point behind
  // the terminator plane.
  amrex::Real grazing_optical_depth(amrex::Real p) const {
    amrex::Real tau = 0.0;
    for (int iC = 0; iC < nExoComponent; ++iC) {
      const amrex::Real sigma = (iC < static_cast<int>(optCrossSection.size()))
                                    ? optCrossSection[iC]
                                    : 0.0;
      if (sigma <= 0.0)
        continue;
      const amrex::Real h = exosphere_scale_height(p, iC);
      if (h <= 0.0)
        continue;
      tau += get_exosphere_component_density(p, iC) * sigma * h *
             std::sqrt(2.0 * M_PI * p / h);
    }
    return tau;
  }

  // Chapman function for a curved atmosphere (Smith & Smith 1972 fits, the
  // same expressions as BATSRUS ModUserMars).  Xp = r/H.
  static amrex::Real chapman_function(amrex::Real Xp, amrex::Real cosSZA) {
    const amrex::Real chap_y = sqrt(0.5 * Xp) * std::abs(cosSZA);
    if (cosSZA > 0.0) {
      // SZA < 90 deg (equation 13)
      if (chap_y < 8.0)
        return sqrt(0.5 * M_PI * Xp) * (1.0606963 + 0.5564383 * chap_y) /
               (1.0619896 + 1.7245609 * chap_y + chap_y * chap_y);
      if (chap_y < 100.0)
        return sqrt(0.5 * M_PI * Xp) * 0.56498823 / (0.6651874 + chap_y);
      return 0.0;
    }
    // 180 > SZA > 90 deg (equation 15)
    const amrex::Real sinSZA = sqrt(std::max(1.0 - cosSZA * cosSZA, 0.0));
    if (chap_y < 8.0)
      return sqrt(2.0 * M_PI * Xp) *
             (sqrt(sinSZA) * exp(Xp * (1.0 - sinSZA)) -
              0.5 * (1.0606963 + 0.5564383 * chap_y) /
                  (1.0619896 + 1.7245609 * chap_y + chap_y * chap_y));
    if (chap_y < 100.0)
      return sqrt(2.0 * M_PI * Xp) *
             (sqrt(sinSZA) * exp(std::min(100.0, Xp * (1.0 - sinSZA))) -
              0.5 * 0.56498823 / (0.6651874 + chap_y));
    return 0.0;
  }

  int chapman_component() const {
    return (chapmanComponent >= 0 && chapmanComponent < nExoComponent)
               ? chapmanComponent
               : 0;
  }

  // Attenuation factor (0..1) applied to the unattenuated photoionization
  // rate of a neutral species at position xyz (r = |xyz|).
  amrex::Real photo_attenuation(const amrex::Real xyz[3], amrex::Real r) const {
    if (!useOpticalDepth)
      return is_in_shadow(xyz[0], xyz[1], xyz[2]) ? 0.0 : 1.0;
    if (is_in_shadow(xyz[0], xyz[1], xyz[2]))
      return optMinProduction;

    const amrex::Real r_safe = std::max(r, 1.0e-3);
    const amrex::Real proj = xyz[0] * optSolarDir[0] + xyz[1] * optSolarDir[1] +
                             xyz[2] * optSolarDir[2];
    const amrex::Real cosSZA = proj / r_safe;
    const amrex::Real mu = std::max(cosSZA, optCosSzaFloor);

    if (useChapmanFunction) {
      const int iC = chapman_component();
      const amrex::Real tauV = vertical_optical_depth(r, iC);
      if (tauV > optChapmanTauMax)
        return optMinProduction * std::max(mu, optMinProduction);
      const amrex::Real h = exosphere_scale_height(r, iC);
      if (h <= 0.0)
        return optMinProduction;
      const amrex::Real chap = chapman_function(r / h, cosSZA);
      // The Smith & Smith fit turns negative deep on the nightside, where
      // sin(SZA) -> 0 and the slant column is in fact optically thick, so
      // exp(-tau*chap) would exceed unity (an unphysical flux *gain*).
      // BATSRUS never evaluates it because its own tau exceeds 13.8 there;
      // for a thin column that guard does not trigger, so treat the point as
      // opaque explicitly.
      if (chap <= 0.0)
        return optMinProduction * std::max(mu, optMinProduction);
      return std::max(exp(-tauV * chap), optMinProduction);
    }

    amrex::Real tau;
    if (proj <= 0.0) {
      // Behind the terminator plane (BATSRUS: dark). That holds for the
      // thermosphere, but a corona extends to several planetary radii: only
      // the geometric shadow of the planet (radius from #SHADOWCYLINDER,
      // otherwise #PLANETRADIUS) is dark, and elsewhere the solar ray grazes
      // the planet at the tangent radius p.
      const amrex::Real p = std::sqrt(std::max(r * r - proj * proj, 0.0));
      const amrex::Real rShadow =
          useShadowCylinder ? shadowCylinderRadius : get_rPlanet_SI();
      if (p <= rShadow)
        return optMinProduction;
      tau = std::max(grazing_optical_depth(p), optTauFloor);
    } else {
      tau = std::max(vertical_optical_depth(r, -1), optTauFloor) / mu;
    }
    if (tau >= optTauCutoff)
      return optMinProduction;
    return std::max(exp(-tau), optMinProduction);
  }

  // ---- Dual-cell averages of the neutral exosphere ----
  // The neutral profiles vary on the scale height (6.5 km for CO2 at Mars),
  // far below the cell size, so a value sampled at a node aliases with the
  // grid: the total production would depend on where the nodes happen to
  // fall. Every node instead carries the average over its dual cell (the
  // cell of size dx centered on the node) of n_c and of n_c * attenuation,
  // plus the fraction of the dual cell outside the inner body. The profiles
  // are static, so the averages are computed once per grid.
  amrex::Real innerRadiusSI = 0.0; // radius of the inner body [m]
  // #INNERBOUNDARY photochemical: below this radius [m] the ions are produced and
  // lost in local photochemical equilibrium, so the ionization there does not
  // feed the plasma (0: no such boundary).
  amrex::Real photochemRadiusSI = 0.0;
  // Per node: [ext_c, extProdAtt_c, extProd_c, frac], see
  // build_neutral_average().
  amrex::Vector<amrex::MultiFab> neuAvg;
  amrex::Vector<amrex::BoxArray> neuAvgGrids;
  static constexpr int maxCubeDepth = 8;

  void set_inner_radius_si(amrex::Real r) override {
    innerRadiusSI = r;
    neuAvg.clear();
    neuAvgGrids.clear();
  }

  void set_photochem_radius_si(amrex::Real r) override {
    photochemRadiusSI = r;
    neuAvg.clear();
    neuAvgGrids.clear();
  }

  // Scale height [m] of the steepest component that is not negligible
  // anywhere in the radial range [rMin, rMax]; a large value if none is.
  amrex::Real steepest_scale_height(amrex::Real rMin) const {
    amrex::Real hMin = 1e300;
    const amrex::Real rRef = std::max(rMin, get_rPlanet_SI());
    for (int iC = 0; iC < nExoComponent; ++iC) {
      const amrex::Real n = get_exosphere_component_density(rRef, iC);
      if (n <= 1e-20 * std::max(exoN0[iC], 1e-300))
        continue;
      const amrex::Real h = exosphere_scale_height(rRef, iC);
      if (h > 0.0)
        hMin = std::min(hMin, h);
    }
    return hMin;
  }

  // Integrate n_c and n_c * attenuation over the box centered at c with
  // sizes h (SI, h[d] = 0 for an inactive dimension), and the volume of the
  // box outside rCut. Boxes that straddle rCut or are coarse compared to the
  // local scale height are split; a leaf uses the exact box average of a
  // planar exponential, prod_d sinh(a_d)/a_d with a_d = h_d n_d / (2 H).
  void integrate_cube(const amrex::Real c[3], const amrex::Real h[3], int depth,
                      amrex::Real rCut, amrex::Real* sumN, amrex::Real* sumNA,
                      amrex::Real& vOut) const {
    const amrex::Real r = std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    const amrex::Real halfDiag =
        0.5 * std::sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
    amrex::Real vol = 1.0;
    for (int d = 0; d < 3; ++d)
      if (h[d] > 0)
        vol *= h[d];
    if (r + halfDiag <= rCut)
      return;

    const bool straddle = r - halfDiag < rCut;
    const amrex::Real hMax = std::max({ h[0], h[1], h[2] });
    const amrex::Real hScale = steepest_scale_height(r - halfDiag);
    // A straddling box is split until it is fine compared to the scale
    // height (its density is evaluated just above rCut), any other box until
    // the planar-exponential average applies.
    const bool split =
        straddle ? hMax > 0.5 * hScale && hScale < 1e299 : hMax > 2.0 * hScale;
    if (split && depth < maxCubeDepth) {
      for (int ic = 0; ic < 8; ++ic) {
        amrex::Real cc[3], hh[3];
        bool skip = false;
        for (int d = 0; d < 3; ++d) {
          const int sgn = ((ic >> d) & 1) ? 1 : -1;
          if (h[d] > 0) {
            hh[d] = 0.5 * h[d];
            cc[d] = c[d] + sgn * 0.25 * h[d];
          } else {
            if (sgn > 0)
              skip = true;
            hh[d] = 0.0;
            cc[d] = c[d];
          }
        }
        if (!skip)
          integrate_cube(cc, hh, depth + 1, rCut, sumN, sumNA, vOut);
      }
      return;
    }

    // Leaf: the part outside rCut, estimated on a 4^nDim sub-grid when the
    // box straddles it.
    amrex::Real fOut = 1.0;
    if (straddle) {
      int nIn = 0, nTot = 0;
      for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
          for (int e = 0; e < 4; ++e) {
            const amrex::Real off[3] = { (a + 0.5) / 4 - 0.5,
                                         (b + 0.5) / 4 - 0.5,
                                         (e + 0.5) / 4 - 0.5 };
            amrex::Real q2 = 0;
            for (int d = 0; d < 3; ++d) {
              const amrex::Real x = c[d] + off[d] * h[d];
              q2 += x * x;
            }
            nTot++;
            if (q2 > rCut * rCut)
              nIn++;
          }
      fOut = amrex::Real(nIn) / nTot;
    }
    if (fOut <= 0)
      return;
    vOut += fOut * vol;

    // Evaluate at the leaf center, or for a straddling leaf just above rCut.
    amrex::Real ce[3] = { c[0], c[1], c[2] };
    amrex::Real re = r;
    if (straddle && r > 0) {
      const amrex::Real rE = std::max(r, rCut + 0.25 * halfDiag);
      for (int d = 0; d < 3; ++d)
        ce[d] = c[d] * rE / r;
      re = rE;
    }
    if (re <= 0)
      return;
    const amrex::Real att = photo_attenuation(ce, re);
    for (int iC = 0; iC < nExoComponent; ++iC) {
      amrex::Real n = get_exosphere_component_density(re, iC);
      if (n <= 0)
        continue;
      const amrex::Real hC = exosphere_scale_height(re, iC);
      if (hC > 0 && !straddle) {
        for (int d = 0; d < 3; ++d) {
          const amrex::Real a = 0.5 * h[d] * std::abs(ce[d] / re) / hC;
          if (a > 1e-6)
            n *= std::sinh(a) / a;
        }
      }
      sumN[iC] += n * fOut * vol;
      sumNA[iC] += n * att * fOut * vol;
    }
  }

  // Fill neuAvg[iLev] for the nodes of nodeFluid[iLev].
  void build_neutral_average(int iLev) {
    timing_func("UserSource::build_neutral_average");
    const int nC = nExoComponent;
    if (static_cast<int>(neuAvg.size()) < n_lev_max()) {
      neuAvg.resize(n_lev_max());
      neuAvgGrids.resize(n_lev_max());
    }
    neuAvg[iLev].define(nodeFluid[iLev].boxArray(),
                        nodeFluid[iLev].DistributionMap(), 3 * nC + 1, 0);
    neuAvgGrids[iLev] = nodeFluid[iLev].boxArray();

    const amrex::Box gbx = convert(Geom(0).Domain(), { AMREX_D_DECL(1, 1, 1) });
    const amrex::Real no2SiL = get_No2SiL();
    const amrex::Real* dx = Geom(iLev).CellSize();
    const auto plo = Geom(iLev).ProbLo();
    const int nDimF = get_fluid_dimension();
    const amrex::Real rCut = std::max(innerRadiusSI, get_rPlanet_SI());
    const amrex::Real rProd = std::max(rCut, photochemRadiusSI);
    amrex::Real hBox[3] = { 0, 0, 0 };
    for (int d = 0; d < nDimF; ++d)
      hBox[d] = dx[d] * no2SiL;
    amrex::Real vBox = 1.0;
    for (int d = 0; d < nDimF; ++d)
      vBox *= hBox[d];
    std::vector<amrex::Real> sumN(nC), sumNA(nC), prodN(nC), prodNA(nC);

    for (amrex::MFIter mfi(neuAvg[iLev]); mfi.isValid(); ++mfi) {
      auto arr = neuAvg[iLev][mfi].array();
      const auto lo = lbound(mfi.validbox());
      const auto hi = ubound(mfi.validbox());
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            amrex::IntVect idx = { AMREX_D_DECL(i, j, k) };
            for (int iDim = 0; iDim < nDim; iDim++)
              if (Geom(iLev).isPeriodic(iDim))
                idx[iDim] = shift_periodic_index(idx[iDim], gbx.smallEnd(iDim),
                                                 gbx.bigEnd(iDim));
            amrex::Real c[3] = { 0, 0, 0 };
            for (int d = 0; d < nDimF; ++d)
              c[d] = (idx[d] * dx[d] + plo[d]) * no2SiL;
            std::fill(sumN.begin(), sumN.end(), 0.0);
            std::fill(sumNA.begin(), sumNA.end(), 0.0);
            amrex::Real vOut = 0;
            // The plasma volume is cut at the body only: without a body the
            // neutral-free region below #PLANETRADIUS still holds plasma.
            integrate_cube(c, hBox, 0, rCut, sumN.data(), sumNA.data(), vOut);
            amrex::Real vPlasma = vOut;
            if (innerRadiusSI <= 0) {
              vPlasma = vBox;
            } else if (rCut > innerRadiusSI) {
              // Body below #PLANETRADIUS: the plasma volume is outside the
              // body, not outside the neutral atmosphere.
              std::vector<amrex::Real> dN(nC), dNA(nC);
              vPlasma = 0;
              integrate_cube(c, hBox, 0, innerRadiusSI, dN.data(), dNA.data(),
                             vPlasma);
            }
            // The ionization that feeds the plasma: above the
            // photochemical boundary only.
            if (rProd > rCut) {
              std::fill(prodN.begin(), prodN.end(), 0.0);
              std::fill(prodNA.begin(), prodNA.end(), 0.0);
              amrex::Real vDummy = 0;
              integrate_cube(c, hBox, 0, rProd, prodN.data(), prodNA.data(),
                             vDummy);
            } else {
              prodN = sumN;
              prodNA = sumNA;
            }
            for (int iC = 0; iC < nC; ++iC) {
              arr(i, j, k, iC) = sumN[iC] / vBox;
              arr(i, j, k, nC + iC) = prodNA[iC] / vBox;
              arr(i, j, k, 2 * nC + iC) = prodN[iC] / vBox;
            }
            arr(i, j, k, 3 * nC) = vPlasma / vBox;
          }
    }
  }

  // Make sure the neutral averages match the current grid of every level.
  // The averages are expensive (minutes on a refined grid), so a load
  // balance, which only moves the boxes between processors, copies them.
  void update_neutral_average() {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      if (nodeFluid[iLev].empty())
        continue;
      const bool cached =
          iLev < static_cast<int>(neuAvg.size()) && !neuAvg[iLev].empty();
      if (cached && neuAvgGrids[iLev] == nodeFluid[iLev].boxArray()) {
        if (neuAvg[iLev].DistributionMap() !=
            nodeFluid[iLev].DistributionMap()) {
          amrex::MultiFab moved(nodeFluid[iLev].boxArray(),
                                nodeFluid[iLev].DistributionMap(),
                                neuAvg[iLev].nComp(), 0);
          moved.ParallelCopy(neuAvg[iLev]);
          neuAvg[iLev] = std::move(moved);
        }
        continue;
      }
      build_neutral_average(iLev);
    }
  }

  // ---- Exosphere density profiles ----

  amrex::Real get_exosphere_density(amrex::Real r) const override {
    amrex::Real sum = 0.0;
    for (int i = 0; i < nExoComponent; ++i)
      sum += get_exosphere_component_density(r, i);
    return sum;
  }

  amrex::Real get_exosphere_component_density(amrex::Real r,
                                              int iC) const override {
    if (exosphereType == "None")
      return 0.0;
    if (iC < 0 || iC >= nExoComponent)
      return 0.0;
    if (r < get_rPlanet_SI())
      return 0.0;

    const std::string& type = exoTypeC[iC];
    if (type == "Exponential") {
      if (exoH0[iC] > 0.0) {
        return exoN0[iC] * exp(-(r - get_rPlanet_SI()) / exoH0[iC]);
      }
    } else if (type == "Power-Law") {
      if (r > 0.0) {
        return exoN0[iC] * pow(get_rPlanet_SI() / r, exoK0[iC]);
      }
    } else if (type == "Chamberlain") {
      if (get_rPlanet_SI() > 0.0 && r > 0.0) {
        return exoN0[iC] * exp(-exoH0[iC] * (1.0 / get_rPlanet_SI() - 1.0 / r));
      }
    }
    return 0.0;
  }

  // Check whether (x, y, z) [m] relative to planet center is in shadow.
  bool is_in_shadow(amrex::Real x, amrex::Real y, amrex::Real z) const {
    if (!useShadowCylinder)
      return false;
    amrex::Real proj = x * solarDir[0] + y * solarDir[1] + z * solarDir[2];
    if (proj >= 0.0)
      return false;
    if (-proj > shadowCylinderHalfHeight)
      return false;
    amrex::Real r2 = x * x + y * y + z * z;
    amrex::Real perp2 = r2 - proj * proj;
    return perp2 <= shadowCylinderRadius * shadowCylinderRadius;
  }

  // Electron temperature in eV from PIC-normalized pressure and density.
  amrex::Real electron_temperature(amrex::Real pe, amrex::Real ne) const {
    if (ne <= 0.0)
      return 0.0;
    const amrex::Real protonMassPerCharge =
        cProtonMassSI / cUnitChargeSI;                       // [kg/C]
    const amrex::Real ur2 = get_unorm_si() * get_unorm_si(); // [m^2/s^2]
    return protonMassPerCharge * ur2 * pe / ne;
  }

  // Photoionization frequency [s^-1] for neutral component iC.
  // 'attenuation' is the factor from photo_attenuation(); pass a negative
  // value to have it computed from the position.
  amrex::Real photoionization_rate(const amrex::Real xyz[3], int iC,
                                   amrex::Real attenuation = -1.0) const {
    if (attenuation >= 0.0) {
      return photoNu0[iC] * attenuation;
    }
    amrex::Real r2 = xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2];
    amrex::Real r_m = sqrt(r2);
    if (r_m <= 0.0)
      return 0.0;
    return photoNu0[iC] * photo_attenuation(xyz, r_m);
  }

  // Electron-impact ionization frequency [s^-1] for neutral component iC.
  amrex::Real impact_ionization_rate(amrex::Real ne, amrex::Real Te_eV,
                                     int iC) const {
    if (ne <= 0.0 || Te_eV <= 0.0)
      return 0.0;
    amrex::Real ne_SI = ne / (get_Si2NoRho() * cProtonMassSI);
    amrex::Real Ei = impactEIon[iC];
    amrex::Real A = impactA[iC];
    amrex::Real K = impactK[iC];
    amrex::Real X = impactX[iC];
    amrex::Real t = Te_eV / Ei;
    amrex::Real sigma_v = A * pow(t, K) / (X + t) * exp(-Ei / Te_eV) * 1e-6;
    return ne_SI * sigma_v;
  }

  // Number density [m^-3] of species iSp in the plasma state `other`, read
  // from its own density slot (0 if it has none).
  amrex::Real species_number_density_si(const FluidInterface& other,
                                        const amrex::MFIter& mfi,
                                        const amrex::IntVect& idx, int iLev,
                                        int iSp) const {
    const int iRho = node_rho_slot(iSp);
    if (iRho < 0)
      return 0.0;
    const amrex::Real rho = other.get_value(mfi, idx, iRho, iLev);
    if (rho <= 0.0)
      return 0.0;
    return rho / get_species_mass(iSp) / (get_Si2NoRho() * cProtonMassSI);
  }

  // Charge-exchange rate coefficient sigma*|u_i| [m^3/s] of ion iSp with
  // neutral component iC (#CHARGEEXCHANGE). The ion loses its charge with the
  // frequency n_n*sigma*|u_i| and a new ion of component iC appears at the
  // neutral velocity and temperature; no electron is created or lost.
  amrex::Real cx_sigma_v(const FluidInterface& other, const amrex::MFIter& mfi,
                         const amrex::IntVect& idx, int iLev, int iSp,
                         int iC) const {
    const int iIon = iSp - 1;
    const amrex::Real ux_i = other.get_ux(mfi, idx, iSp, iLev);
    const amrex::Real uy_i = other.get_uy(mfi, idx, iSp, iLev);
    const amrex::Real uz_i = other.get_uz(mfi, idx, iSp, iLev);
    const amrex::Real u_mag_SI =
        sqrt(ux_i * ux_i + uy_i * uy_i + uz_i * uz_i) * get_unorm_si();
    // [cm^2] -> [m^2]
    return cxSigma[iC * nCXIonSpecies + iIon] * 1e-4 * u_mag_SI;
  }

  // Fraction of the linear-rate loss nu*dt that an exponential decay over dt
  // really removes, (1 - exp(-nu dt))/(nu dt). The products of a reaction are
  // scaled with it so that they carry exactly the mass the reactant particles
  // lose in apply_loss(), also for stiff rates (nu*dt >> 1).
  static amrex::Real stiff_factor(amrex::Real nuDt) {
    if (nuDt < 1e-8)
      return 1.0;
    return -std::expm1(-nuDt) / nuDt;
  }

  //-------------------------------------------------------------------
  // Require #EXOSPHERE to have been read (with nComponent > 0) before an
  // ionization command is parsed.  Domain syncs nExoComponent from fi into
  // this object just before dispatching ionization commands, so a missing or
  // out-of-order #EXOSPHERE leaves nExoComponent == 0 and is caught here.
  void require_exosphere(const std::string& command) const {
    if (nExoComponent > 0)
      return;
    amrex::Abort(printPrefix + "Error: " + command + " requires " +
                 "#EXOSPHERE with nComponent > 0 to be specified before " +
                 command + ".");
  }

  //-------------------------------------------------------------------
  // Read ionization-related parameter commands from PARAM.in.
  void read_param(const std::string& command, ReadParam& param) override {
    if (command == "#PHOTOIONIZATION") {
      usePhotoIonization = true;
      // The neutral component count is inherited from #EXOSPHERE.
      require_exosphere(command);
      photoNu0.resize(nExoComponent);
      for (int i = 0; i < nExoComponent; ++i) {
        param.read_var("nuPhoto0", photoNu0[i]);
      }
    } else if (command == "#ELECTRONIMPACT") {
      useElectronImpact = true;
      // The neutral component count is inherited from #EXOSPHERE.
      require_exosphere(command);
      impactEIon.resize(nExoComponent);
      impactA.resize(nExoComponent);
      impactK.resize(nExoComponent);
      impactX.resize(nExoComponent);
      for (int i = 0; i < nExoComponent; ++i) {
        param.read_var("eIon", impactEIon[i]);
        param.read_var("Acoeff", impactA[i]);
        param.read_var("Kcoeff", impactK[i]);
        param.read_var("Xcoeff", impactX[i]);
      }
    } else if (command == "#CHARGEEXCHANGE") {
      useChargeExchange = true;
      // The neutral component count is inherited from #EXOSPHERE.
      require_exosphere(command);
      param.read_var("nIonSpecies", nCXIonSpecies);
      cxSigma.resize(nExoComponent * nCXIonSpecies);
      for (int iC = 0; iC < nExoComponent; ++iC) {
        for (int iIon = 0; iIon < nCXIonSpecies; ++iIon) {
          param.read_var("sigmaCX", cxSigma[iC * nCXIonSpecies + iIon]);
        }
      }
    } else if (command == "#SHADOWCYLINDER") {
      useShadowCylinder = true;
      param.read_var("solarDirX", solarDir[0]);
      param.read_var("solarDirY", solarDir[1]);
      param.read_var("solarDirZ", solarDir[2]);
      param.read_var("radius", shadowCylinderRadius);
      param.read_var("halfHeight", shadowCylinderHalfHeight);
      amrex::Real norm =
          sqrt(solarDir[0] * solarDir[0] + solarDir[1] * solarDir[1] +
               solarDir[2] * solarDir[2]);
      if (norm > 0.0) {
        solarDir[0] /= norm;
        solarDir[1] /= norm;
        solarDir[2] /= norm;
      }
    } else if (command == "#OPTICALDEPTH") {
      useOpticalDepth = true;
      // The neutral component count is inherited from #EXOSPHERE.
      require_exosphere(command);
      // Per-component EUV absorption cross section [m^2] and, right after,
      // per-component scale height [m].  A negative scale height means "use
      // the local scale height of the #EXOSPHERE profile".
      optCrossSection.assign(nExoComponent, 0.0);
      optScaleHeight.assign(nExoComponent, -1.0);
      for (int i = 0; i < nExoComponent; ++i) {
        param.read_var("crossSection", optCrossSection[i]);
      }
      for (int i = 0; i < nExoComponent; ++i) {
        param.read_var("scaleHeight", optScaleHeight[i]);
      }
      // #OPTICALDEPTH carries its own solar direction [unit vector], kept
      // separate from the #SHADOWCYLINDER one.
      param.read_var("solarDirX", optSolarDir[0]);
      param.read_var("solarDirY", optSolarDir[1]);
      param.read_var("solarDirZ", optSolarDir[2]);
      amrex::Real normOpt = sqrt(optSolarDir[0] * optSolarDir[0] +
                                 optSolarDir[1] * optSolarDir[1] +
                                 optSolarDir[2] * optSolarDir[2]);
      if (normOpt > 0.0) {
        optSolarDir[0] /= normOpt;
        optSolarDir[1] /= normOpt;
        optSolarDir[2] /= normOpt;
      }
      // The remaining parameters are individually optional and default to the
      // BATSRUS values.  The T/F flag is read as a string because
      // read_optional() parses numbers only.
      std::string chapmanFlag;
      if (param.read_optional("chapmanFunction", chapmanFlag)) {
        useChapmanFunction =
            (chapmanFlag == "T" || chapmanFlag == "t" || chapmanFlag == "true");
      }
      if (param.read_optional("minProduction", optMinProduction)) {
        optMinProductionSet = true;
      }
      param.read_optional("tauFloor", optTauFloor);
      param.read_optional("tauCutoff", optTauCutoff);
      param.read_optional("cosSzaFloor", optCosSzaFloor);
      param.read_optional("chapmanTauMax", optChapmanTauMax);
      param.read_optional("chapmanComponent", chapmanComponent);
    } else if (command == "#RECOMBINATION") {
      useRecombination = true;
      int nRecomb;
      param.read_var("nReactions", nRecomb);
      recombIonIndex.resize(nRecomb);
      recombRate0.resize(nRecomb);
      recombTempExp.resize(nRecomb);
      recombRefTemp.resize(nRecomb);
      for (int i = 0; i < nRecomb; ++i) {
        param.read_var("ionSpecies", recombIonIndex[i]);
        param.read_var("rateCoef", recombRate0[i]);
        param.read_var("tempExponent", recombTempExp[i]);
        param.read_var("refTemp", recombRefTemp[i]);
      }
    } else if (command == "#PHOTOCHEMICALBOUNDARY") {
      param.read_var("radius", photochemRadiusSI);
      if (photochemRadiusSI < 0)
        amrex::Abort(printPrefix +
                     "Error: #PHOTOCHEMICALBOUNDARY radius must be >= 0.");
      amrex::Print() << "  Warning: #PHOTOCHEMICALBOUNDARY is deprecated; "
                     << "use #INNERBOUNDARY with the 'photochemical' keyword instead.\n";
      neuAvg.clear();
      neuAvgGrids.clear();
    } else if (command == "#CHEMISTRY") {
      useChemistry = true;
      int nRxns;
      param.read_var("nReactions", nRxns);
      chemReactions.resize(nRxns);
      for (int i = 0; i < nRxns; ++i) {
        param.read_var("reactantIon", chemReactions[i].reactantIon);
        param.read_var("productIon", chemReactions[i].productIon);
        param.read_var("neutralComp", chemReactions[i].neutralComp);
        param.read_var("rateType", chemReactions[i].rateType);
        param.read_var("rateCoef", chemReactions[i].rateCoef);
        param.read_var("tempExp", chemReactions[i].tempExp);
        param.read_var("refTemp", chemReactions[i].refTemp);
      }
    }
  }

  // Validate consistency between exosphere and ionization commands.
  void post_process_param() override {
    // Optical-depth attenuation validation.
    if (useOpticalDepth) {
      if (useShadowCylinder) {
        amrex::Abort(printPrefix + "Error: #OPTICALDEPTH and #SHADOWCYLINDER " +
                     "describe two alternative models of the same physics " +
                     "(EUV attenuation) and cannot be combined. Keep one.");
      }
      if (nExoComponent <= 0) {
        amrex::Abort(printPrefix + "Error: #OPTICALDEPTH requires #EXOSPHERE " +
                     "with nComponent > 0 to be specified before it.");
      }
      const amrex::Real dirNorm = sqrt(optSolarDir[0] * optSolarDir[0] +
                                       optSolarDir[1] * optSolarDir[1] +
                                       optSolarDir[2] * optSolarDir[2]);
      if (dirNorm <= 0.0) {
        amrex::Abort(printPrefix + "Error: #OPTICALDEPTH requires a non-zero " +
                     "solar direction (solarDirX, solarDirY, solarDirZ).");
      }
      if (optCosSzaFloor <= 0.0 || optCosSzaFloor >= 1.0) {
        amrex::Abort(printPrefix + "Error: #OPTICALDEPTH cosSzaFloor must be " +
                     "in (0, 1). Got " + std::to_string(optCosSzaFloor) + ".");
      }
      if (optTauFloor < 0.0) {
        amrex::Abort(printPrefix + "Error: #OPTICALDEPTH tauFloor must be " +
                     ">= 0. Got " + std::to_string(optTauFloor) + ".");
      }
      if (optTauCutoff <= 0.0) {
        amrex::Abort(printPrefix + "Error: #OPTICALDEPTH tauCutoff must be " +
                     "> 0. Got " + std::to_string(optTauCutoff) + ".");
      }
      if (optMinProduction < 0.0 || optMinProduction > 1.0) {
        amrex::Abort(
            printPrefix + "Error: #OPTICALDEPTH minProduction must be " +
            "in [0, 1]. Got " + std::to_string(optMinProduction) + ".");
      }
      const amrex::Real rP = get_rPlanet_SI();
      int nAbsorbing = 0;
      for (int iC = 0; iC < nExoComponent; ++iC) {
        if (optCrossSection[iC] < 0.0) {
          amrex::Abort(printPrefix + "Error: #OPTICALDEPTH crossSection for " +
                       "component " + std::to_string(iC) +
                       " must be >= 0 [m^2]. Got " +
                       std::to_string(optCrossSection[iC]) + ".");
        }
        if (optCrossSection[iC] <= 0.0)
          continue;
        nAbsorbing++;
        if (exosphere_scale_height(rP, iC) <= 0.0) {
          amrex::Abort(printPrefix + "Error: #OPTICALDEPTH component " +
                       std::to_string(iC) + " has a positive crossSection " +
                       "but no usable scale height. Give an explicit " +
                       "scaleHeight, or use an #EXOSPHERE profile that " +
                       "defines one (Exponential with H0 > 0).");
        }
      }
      if (nAbsorbing == 0) {
        amrex::Abort(printPrefix + "Error: #OPTICALDEPTH needs at least one " +
                     "positive crossSection, otherwise the optical depth is " +
                     "always the tauFloor.");
      }
      // BATSRUS uses a lower production floor with the Chapman function.
      if (useChapmanFunction && !optMinProductionSet)
        optMinProduction = 1.0e-6;
      // Resolve the component carrying the Chapman slant column (the one with
      // the largest vertical optical depth, CO2 for Mars).
      if (useChapmanFunction && chapmanComponent < 0) {
        amrex::Real bestTau = -1.0;
        chapmanComponent = 0;
        for (int iC = 0; iC < nExoComponent; ++iC) {
          const amrex::Real h = exosphere_scale_height(rP, iC);
          if (optCrossSection[iC] <= 0.0 || h <= 0.0)
            continue;
          const amrex::Real tau =
              get_exosphere_component_density(rP, iC) * optCrossSection[iC] * h;
          if (tau > bestTau) {
            bestTau = tau;
            chapmanComponent = iC;
          }
        }
      }
    }

    // Recombination validation.
    if (useRecombination) {
      if (nS < 1) {
        amrex::Abort(printPrefix + "Error: #RECOMBINATION requires " +
                     "plasma species. Use #PLASMA to set species.");
      }
      if (!has_electron_species()) {
        amrex::Abort(printPrefix + "Error: #RECOMBINATION requires species 0 " +
                     "to be the electron (negative charge). Reorder " +
                     "#PLASMA so that the electron comes first. The " +
                     "electron may be quasi-neutral " +
                     "(useElectronFluid = false), which is the case for " +
                     "GM-coupled multi-species runs.");
      }
      const int nIonS = nS - 1;
      for (int i = 0; i < static_cast<int>(recombIonIndex.size()); ++i) {
        int iSp = recombIonIndex[i];
        if (iSp < 1 || iSp > nIonS) {
          amrex::Abort(printPrefix + "Error: #RECOMBINATION ionSpecies " +
                       std::to_string(iSp) + " is out of range [1, " +
                       std::to_string(nIonS) + "].");
        }
        if (recombRate0[i] <= 0.0) {
          amrex::Abort(printPrefix + "Error: #RECOMBINATION rateCoef must " +
                       "be positive (got " + std::to_string(recombRate0[i]) +
                       ").");
        }
      }
    }

    // Chemistry validation.
    if (useChemistry) {
      if (usePhotoIonization || useChargeExchange || useRecombination) {
        amrex::Abort(
            printPrefix + "Error: #CHEMISTRY cannot be combined with " +
            "#PHOTOIONIZATION, #CHARGEEXCHANGE, or #RECOMBINATION. " +
            "Use either the unified #CHEMISTRY table or individual commands.");
      }
      if (nS < 2) {
        amrex::Abort(printPrefix + "Error: #CHEMISTRY requires at least " +
                     "2 plasma species (electron + 1 ion).");
      }
      const int nIonS = nS - 1;
      bool needsNeutral = false;
      bool needsElectron = false;
      for (int i = 0; i < static_cast<int>(chemReactions.size()); ++i) {
        const auto& rxn = chemReactions[i];
        if (rxn.reactantIon < 0 || rxn.reactantIon > nIonS) {
          amrex::Abort(printPrefix + "Error: #CHEMISTRY reactantIon " +
                       std::to_string(rxn.reactantIon) + " out of range " +
                       "[0, " + std::to_string(nIonS) + "].");
        }
        if (rxn.productIon < 0 || rxn.productIon > nIonS) {
          amrex::Abort(printPrefix + "Error: #CHEMISTRY productIon " +
                       std::to_string(rxn.productIon) + " out of range " +
                       "[0, " + std::to_string(nIonS) + "].");
        }
        if (rxn.reactantIon == 0 && rxn.productIon == 0) {
          amrex::Abort(printPrefix + "Error: #CHEMISTRY reaction " +
                       std::to_string(i) + " has both reactantIon and " +
                       "productIon = 0 (no-op).");
        }
        if (rxn.rateType == 1 && rxn.neutralComp < 0) {
          amrex::Abort(printPrefix + "Error: #CHEMISTRY photoionization " +
                       "reaction " + std::to_string(i) +
                       " requires a neutral component.");
        }
        if (rxn.rateType == 0 && rxn.neutralComp < 0 && rxn.reactantIon == 0) {
          amrex::Abort(printPrefix + "Error: #CHEMISTRY thermal reaction " +
                       std::to_string(i) + " with no neutral and no " +
                       "reactant ion is invalid.");
        }
        if (rxn.neutralComp >= 0) {
          needsNeutral = true;
          if (nExoComponent > 0 && rxn.neutralComp >= nExoComponent) {
            amrex::Abort(printPrefix + "Error: #CHEMISTRY reaction " +
                         std::to_string(i) + " neutralComp " +
                         std::to_string(rxn.neutralComp) +
                         " >= nExoComponent " + std::to_string(nExoComponent) +
                         ".");
          }
        }
        if (rxn.productIon == 0 || rxn.tempExp != 0.0) {
          needsElectron = true;
        }
      }
      if (needsNeutral && nExoComponent <= 0) {
        amrex::Abort(printPrefix +
                     "Error: #CHEMISTRY reactions with neutralComp >= 0 " +
                     "require #EXOSPHERE with nComponent > 0 to be specified.");
      }
      if (needsElectron && !has_electron_species()) {
        amrex::Abort(printPrefix + "Error: #CHEMISTRY reactions with " +
                     "recombination or temperature dependence require " +
                     "species 0 to be the electron (negative charge). " +
                     "Reorder #PLASMA so that the electron comes first. " +
                     "The electron may be quasi-neutral " +
                     "(useElectronFluid = false), which is the case for " +
                     "GM-coupled multi-species runs.");
      }
    }

    if (exosphereType == "None")
      return;

    if (nS < 1) {
      amrex::Abort(printPrefix + "Error: no plasma species defined. " +
                   "Use #PLASMA to set species.");
    }
    if (QoQi_S[0] >= 0.0) {
      amrex::Abort(
          printPrefix + "Error: species 0 must be the electron " +
          "(negative charge). Got Q/Qi[0] = " + std::to_string(QoQi_S[0]) +
          ". Reorder #PLASMA so the electron is first.");
    }
    // #PHOTOIONIZATION, #ELECTRONIMPACT and #CHARGEEXCHANGE turn component iC
    // into ion species iC + 1. #CHEMISTRY names the product of every reaction
    // explicitly, so any number of components may feed the same ion (e.g. the
    // cold and the hot oxygen corona).
    const bool mapsComponentToIon =
        usePhotoIonization || useElectronImpact || useChargeExchange;
    if (mapsComponentToIon && nExoComponent > nS - 1) {
      amrex::Abort(printPrefix + "Error: #EXOSPHERE nComponent (" +
                   std::to_string(nExoComponent) +
                   ") exceeds the number of ion species (" +
                   std::to_string(nS - 1) +
                   "). Add more ion species in #PLASMA.");
    }

    if (useChargeExchange) {
      if (nCXIonSpecies != nS - 1) {
        amrex::Abort(printPrefix + "Error: #CHARGEEXCHANGE nIonSpecies (" +
                     std::to_string(nCXIonSpecies) +
                     ") != number of ion species (" + std::to_string(nS - 1) +
                     ")");
      }
    }

    for (int iC = 0; iC < nExoComponent; ++iC) {
      if (exoT0[iC] <= 0.0) {
        amrex::Abort(printPrefix + "Error: #EXOSPHERE T0 for component " +
                     std::to_string(iC) + " must be positive (got " +
                     std::to_string(exoT0[iC]) +
                     " K). It is used as the ionization source temperature.");
      }
    }
  }

  // Thermal rate coefficient k(Te) [m^3/s] of a chemistry reaction.
  amrex::Real chem_rate_coef(const ChemistryReaction& rxn,
                             amrex::Real Te_eV) const {
    amrex::Real k_si = rxn.rateCoef * 1e-6; // [cm^3/s] -> [m^3/s]
    if (rxn.tempExp != 0.0 && Te_eV > 0.0) {
      const amrex::Real Te_K = Te_eV * cUnitChargeSI / cBoltzmannSI;
      k_si *= pow(rxn.refTemp / Te_K, rxn.tempExp);
    }
    return k_si;
  }

  // Reaction frequency [1/s] of one reactant ion: the photo rate with the
  // attenuation att, or k(Te) times the neutral density nNeu [m^-3] or, for a
  // reaction without neutral, the electron density neSI [m^-3].
  amrex::Real chem_frequency(const ChemistryReaction& rxn, amrex::Real nNeu,
                             amrex::Real att, amrex::Real neSI,
                             amrex::Real Te_eV) const {
    if (rxn.rateType == 1)
      return rxn.rateCoef * att;
    return chem_rate_coef(rxn, Te_eV) * (rxn.neutralComp >= 0 ? nNeu : neSI);
  }

  //-------------------------------------------------------------------
  // Add the products of a chemistry reaction. sRate [m^-3 s^-1] is the number
  // of reactions per volume and time (already averaged over the dual cell).
  void chem_apply_source(const ChemistryReaction& rxn, amrex::Real sRate,
                         const FluidInterface& other, const amrex::MFIter& mfi,
                         const amrex::IntVect& idx, int iLev, int nIonS,
                         std::vector<amrex::Real>& srcRho,
                         std::vector<amrex::Real>& srcP,
                         std::vector<amrex::Real>& srcRhoUx,
                         std::vector<amrex::Real>& srcRhoUy,
                         std::vector<amrex::Real>& srcRhoUz) const {
    if (rxn.productIon <= 0 || rxn.productIon > nIonS || sRate <= 0)
      return;

    const int iSpProd = rxn.productIon;
    const amrex::Real mass_prod = get_species_mass(iSpProd);
    const amrex::Real srcRho_si = sRate * mass_prod * cProtonMassSI;
    srcRho[iSpProd] += srcRho_si;

    if (rxn.reactantIon > 0 && rxn.reactantIon <= nIonS) {
      // Ion-neutral reaction: the product inherits the velocity and the
      // temperature of the reactant.
      const int iSpReac = rxn.reactantIon;
      srcRhoUx[iSpProd] += srcRho_si * other.get_ux(mfi, idx, iSpReac, iLev);
      srcRhoUy[iSpProd] += srcRho_si * other.get_uy(mfi, idx, iSpReac, iLev);
      srcRhoUz[iSpProd] += srcRho_si * other.get_uz(mfi, idx, iSpReac, iLev);
      // get_p() returns the pressure of species iSpReac itself (it applies
      // the PeRatio share and the multi-species split).
      const amrex::Real nReac =
          species_number_density_si(other, mfi, idx, iLev, iSpReac);
      if (nReac > 0) {
        const amrex::Real p_reac_si = amrex::max(
            0.0, other.get_p(mfi, idx, iSpReac, iLev) / get_Si2NoP());
        srcP[iSpProd] += sRate * p_reac_si / nReac;
      }
    } else {
      // Ionization of a neutral at rest: a new ion-electron pair at the
      // neutral temperature.
      const amrex::Real kT =
          (rxn.neutralComp >= 0 && rxn.neutralComp < nExoComponent)
              ? cBoltzmannSI * exoT0[rxn.neutralComp]
              : 0.0;
      srcP[iSpProd] += sRate * kT;
      if (nS > 0 && get_species_charge(0) < 0) {
        srcRho[0] += sRate * get_species_mass(0) * cProtonMassSI;
        srcP[0] += sRate * kT;
      }
    }
  }

  // Recombination rate coefficient k(Te) [m^3/s] of #RECOMBINATION entry iR.
  amrex::Real recombination_coef(int iR, amrex::Real Te_eV) const {
    amrex::Real k_si = recombRate0[iR] * 1e-6;
    if (Te_eV > 0.0 && recombTempExp[iR] != 0.0) {
      const amrex::Real Te_K = Te_eV * cUnitChargeSI / cBoltzmannSI;
      k_si *= pow(recombRefTemp[iR] / Te_K, recombTempExp[iR]);
    }
    return k_si;
  }

  // Set nodeFluid from plasma-state-dependent ionization processes.
  void set_source(const FluidInterface& other) override {
    timing_func("UserSource::set_source");
    set_node_fluid(other);
    set_node_loss_fluid_to_zero();

    const amrex::Box gbx = convert(Geom(0).Domain(), { AMREX_D_DECL(1, 1, 1) });

    const bool doPhoto = usePhotoIonization;
    const bool doImpact = useElectronImpact;
    const bool doCX = useChargeExchange;
    const bool doRecomb = useRecombination;
    const bool doChem = useChemistry;
    const bool doLoss = use_loss_source();

    const amrex::Real rhoNormPerT = get_Si2NoRho() / get_Si2NoT();
    const amrex::Real pNormPerT = get_Si2NoP() / get_Si2NoT();
    // [m] per FLEKS length unit: 1 for a stand-alone run (the grid is then
    // given in metres) and one planetary radius for a GM-coupled run.
    const int nIonS = nS - 1;
    const int nC = nExoComponent;

    update_neutral_average();

    for (int iLev = 0; iLev < n_lev(); iLev++) {
      if (!nodeFluid[iLev].empty()) {
        for (amrex::MFIter mfi(nodeFluid[iLev]); mfi.isValid(); ++mfi) {
          const amrex::Box& box = mfi.validbox();
          const auto lo = lbound(box);
          const auto hi = ubound(box);

          const amrex::Array4<amrex::Real>& arr = nodeFluid[iLev][mfi].array();
          const auto neu = neuAvg[iLev][mfi].const_array();

          amrex::Array4<amrex::Real> lossArr;
          if (doLoss) {
            lossArr = nodeLossFluid[iLev][mfi].array();
          }

          // Source accumulators in SI units.
          std::vector<amrex::Real> srcRho(nIonS + 1, 0.0);
          std::vector<amrex::Real> srcP(nIonS + 1, 0.0);
          std::vector<amrex::Real> srcRhoUx(nIonS + 1, 0.0);
          std::vector<amrex::Real> srcRhoUy(nIonS + 1, 0.0);
          std::vector<amrex::Real> srcRhoUz(nIonS + 1, 0.0);
          // Loss frequency [1/s] of each species, the part of it that also
          // removes an electron (recombination), and the stiff factor.
          std::vector<amrex::Real> nuLoss(nIonS + 1, 0.0);
          std::vector<amrex::Real> nuRecomb(nIonS + 1, 0.0);
          std::vector<amrex::Real> stiff(nIonS + 1, 1.0);

          for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
              for (int i = lo.x; i <= hi.x; ++i) {
                amrex::IntVect idx = { AMREX_D_DECL(i, j, k) };
                for (int iDim = 0; iDim < nDim; iDim++) {
                  if (Geom(iLev).isPeriodic(iDim)) {
                    idx[iDim] = shift_periodic_index(
                        idx[iDim], gbx.smallEnd(iDim), gbx.bigEnd(iDim));
                  }
                }
                // Neutral exosphere of the node, averaged over its dual cell
                // (see build_neutral_average): ext = <n_c>, extAtt =
                // <n_c * attenuation>, frac = plasma (outside-body) volume
                // fraction. Productions per volume use the dual-cell averages,
                // the frequency of a reactant ion the density in the plasma
                // part, nBar = ext / frac.
                // extProd/extProdAtt: the same above the photochemical
                // boundary, for the ionization that feeds the plasma.
                const amrex::Real frac = neu(i, j, k, 3 * nC);
                auto ext = [&](int iC) { return neu(i, j, k, iC); };
                auto extProdAtt = [&](int iC) { return neu(i, j, k, nC + iC); };
                auto extProd = [&](int iC) {
                  return neu(i, j, k, 2 * nC + iC);
                };
                auto nBar = [&](int iC) {
                  return frac > 0 ? ext(iC) / frac : 0.0;
                };
                // Mean attenuation of the ionizing flux for the photo
                // frequency of an ion (rare; no photochemical cut needed).
                auto attBar = [&](int iC) {
                  return extProd(iC) > 0 ? extProdAtt(iC) / extProd(iC) : 0.0;
                };

                amrex::Real ne = 0, pe = 0, Te_eV = 0;
                bool plasmaFetched = false;
                auto fetch_electron_plasma = [&, &other = other]() {
                  if (plasmaFetched)
                    return;
                  ne = other.get_number_density(mfi, idx, 0, iLev);
                  pe = other.get_p(mfi, idx, 0, iLev);
                  Te_eV = electron_temperature(pe, ne);
                  plasmaFetched = true;
                };

                std::fill(srcRho.begin(), srcRho.end(), 0.0);
                std::fill(srcP.begin(), srcP.end(), 0.0);
                std::fill(srcRhoUx.begin(), srcRhoUx.end(), 0.0);
                std::fill(srcRhoUy.begin(), srcRhoUy.end(), 0.0);
                std::fill(srcRhoUz.begin(), srcRhoUz.end(), 0.0);
                std::fill(nuLoss.begin(), nuLoss.end(), 0.0);
                std::fill(nuRecomb.begin(), nuRecomb.end(), 0.0);

                // ---- Loss frequencies of the ions [1/s] ----
                // They are computed first because the products of a reaction
                // are scaled with the stiff factor of their reactant.
                if (doCX) {
                  for (int iC = 0; iC < nExoComponent; ++iC) {
                    const amrex::Real nn = nBar(iC);
                    if (nn <= 0.0)
                      continue;
                    for (int iSpCX = 1; iSpCX <= nIonS; ++iSpCX)
                      nuLoss[iSpCX] +=
                          nn * cx_sigma_v(other, mfi, idx, iLev, iSpCX, iC);
                  }
                }

                if (doRecomb) {
                  fetch_electron_plasma();
                  const amrex::Real neSI =
                      ne / (get_Si2NoRho() * cProtonMassSI);
                  for (int iR = 0; iR < static_cast<int>(recombIonIndex.size());
                       ++iR) {
                    const int iSp = recombIonIndex[iR];
                    if (iSp < 1 || iSp > nIonS)
                      continue;
                    const amrex::Real nu = recombination_coef(iR, Te_eV) * neSI;
                    nuLoss[iSp] += nu;
                    nuRecomb[iSp] += nu;
                  }
                }

                if (doChem) {
                  fetch_electron_plasma();
                  const amrex::Real neSI =
                      ne / (get_Si2NoRho() * cProtonMassSI);
                  for (const auto& rxn : chemReactions) {
                    if (rxn.reactantIon <= 0 || rxn.reactantIon > nIonS)
                      continue;
                    const int iC = rxn.neutralComp;
                    const amrex::Real rate =
                        chem_frequency(rxn, iC >= 0 ? nBar(iC) : 0.0,
                                       iC >= 0 ? attBar(iC) : 1.0, neSI, Te_eV);
                    if (rate <= 0.0)
                      continue;
                    nuLoss[rxn.reactantIon] += rate;
                    if (rxn.productIon == 0)
                      nuRecomb[rxn.reactantIon] += rate;
                  }
                }

                for (int iSp = 1; iSp <= nIonS; ++iSp)
                  stiff[iSp] = stiff_factor(nuLoss[iSp] * dtSI);

                // ---- Exosphere-based ionization sources ----
                for (int iC = 0; iC < nExoComponent; ++iC) {
                  const int iSp = iC + 1;
                  if (iSp > nIonS)
                    break;

                  const amrex::Real dens_i = ext(iC);
                  if (dens_i <= 0.0)
                    continue;

                  // Ionization: a new ion-electron pair [m^-3 s^-1].
                  amrex::Real S_ion = 0.0;
                  if (doPhoto)
                    S_ion += photoNu0[iC] * extProdAtt(iC);
                  if (doImpact) {
                    fetch_electron_plasma();
                    S_ion +=
                        impact_ionization_rate(ne, Te_eV, iC) * extProd(iC);
                  }

                  // Charge exchange: every ion that loses its charge to a
                  // neutral of component iC leaves a new ion iSp behind.
                  amrex::Real S_cx = 0.0;
                  if (doCX) {
                    for (int iSpCX = 1; iSpCX <= nIonS; ++iSpCX) {
                      const amrex::Real ni = species_number_density_si(
                          other, mfi, idx, iLev, iSpCX);
                      if (ni > 0.0)
                        S_cx += dens_i * ni * stiff[iSpCX] *
                                cx_sigma_v(other, mfi, idx, iLev, iSpCX, iC);
                    }
                  }

                  const amrex::Real S_n = S_ion + S_cx;
                  amrex::Real mass_amu = get_species_mass(iSp);
                  srcRho[iSp] += S_n * mass_amu * cProtonMassSI;
                  srcP[iSp] += S_n * cBoltzmannSI * exoT0[iC];

                  // Co-create neutralizing electrons in full-PIC mode.
                  if (S_ion > 0 && nS > 0 && get_species_charge(0) < 0) {
                    amrex::Real mass_e = get_species_mass(0);
                    srcRho[0] += S_ion * mass_e * cProtonMassSI;
                    srcP[0] += S_ion * cBoltzmannSI * exoT0[iC];
                  }
                }

                // Chemistry sources.
                if (doChem) {
                  fetch_electron_plasma();
                  const amrex::Real neSI =
                      ne / (get_Si2NoRho() * cProtonMassSI);
                  for (const auto& rxn : chemReactions) {
                    const int iC = rxn.neutralComp;
                    amrex::Real sRate = 0.0; // reactions [m^-3 s^-1]
                    if (rxn.reactantIon > 0 && rxn.reactantIon <= nIonS) {
                      // Ion reactant: its frequency in the plasma part of the
                      // dual cell, the stiff fraction that really reacts.
                      const amrex::Real nReac = species_number_density_si(
                          other, mfi, idx, iLev, rxn.reactantIon);
                      sRate = chem_frequency(rxn, iC >= 0 ? nBar(iC) : 0.0,
                                             iC >= 0 ? attBar(iC) : 1.0, neSI,
                                             Te_eV) *
                              stiff[rxn.reactantIon] * nReac * frac;
                    } else if (iC >= 0) {
                      // Ionization of the neutral: photo or by the electrons.
                      sRate = rxn.rateType == 1 ? rxn.rateCoef * extProdAtt(iC)
                                                : chem_rate_coef(rxn, Te_eV) *
                                                      neSI * extProd(iC);
                    }
                    chem_apply_source(rxn, sRate, other, mfi, idx, iLev, nIonS,
                                      srcRho, srcP, srcRhoUx, srcRhoUy,
                                      srcRhoUz);
                  }
                }

                // Write accumulated sources to nodeFluid.
                // set_node_fluid(other) copied the MHD state into nodeFluid,
                // so every slot that may hold a source has to be zeroed
                // first; otherwise a node without a source would re-inject
                // the ambient MHD density as a source.
                for (int iFluid = 0; iFluid < nFluid; iFluid++) {
                  arr(i, j, k, iRho_I[iFluid]) = 0;
                  arr(i, j, k, iUx_I[iFluid]) = 0;
                  arr(i, j, k, iUy_I[iFluid]) = 0;
                  arr(i, j, k, iUz_I[iFluid]) = 0;
                  arr(i, j, k, iP_I[iFluid]) = 0;
                }

                const bool sharedState = shared_ion_state();
                if (sharedState) {
                  // nFluid is 1 for a multi-species MHD state while the ion
                  // species occupy nIon separate density slots.
                  for (int iIon = 0; iIon < nIon; ++iIon)
                    arr(i, j, k, iRho_I[iIon]) = 0;
                }

                // Electron source.  The electron only has a state slot of its
                // own when it is a fluid.  In a GM-coupled quasi-neutral run
                // the electrons are created implicitly with the ions, so
                // nothing is written for species 0 (its pressure is implied
                // by the PeRatio split of the shared pressure slot).
                const int iRhoElectron = node_rho_slot(0);
                if (iRhoElectron >= 0 && srcRho[0] > 0) {
                  arr(i, j, k, iRhoElectron) = srcRho[0] * rhoNormPerT;
                  arr(i, j, k, iP_I[0]) = srcP[0] * pNormPerT;
                  arr(i, j, k, iUx_I[0]) = 0.0;
                  arr(i, j, k, iUy_I[0]) = 0.0;
                  arr(i, j, k, iUz_I[0]) = 0.0;
                }

                bool anySource = false;
                amrex::Real sumRhoUx = 0.0, sumRhoUy = 0.0, sumRhoUz = 0.0;
                amrex::Real sumP = 0.0;
                for (int iSp = 1; iSp <= nIonS; ++iSp) {
                  if (srcRho[iSp] > 0) {
                    anySource = true;
                    const int iRhoSp = node_rho_slot(iSp);
                    if (iRhoSp >= 0)
                      arr(i, j, k, iRhoSp) = srcRho[iSp] * rhoNormPerT;
                    sumRhoUx += srcRhoUx[iSp];
                    sumRhoUy += srcRhoUy[iSp];
                    sumRhoUz += srcRhoUz[iSp];
                    sumP += srcP[iSp];
                    if (!sharedState) {
                      arr(i, j, k, iUx_I[iSp]) = srcRhoUx[iSp] * rhoNormPerT;
                      arr(i, j, k, iUy_I[iSp]) = srcRhoUy[iSp] * rhoNormPerT;
                      arr(i, j, k, iUz_I[iSp]) = srcRhoUz[iSp] * rhoNormPerT;
                      arr(i, j, k, iP_I[iSp]) = srcP[iSp] * pNormPerT;
                    }
                  }
                }

                if (sharedState && anySource) {
                  // All ion species share one momentum and one pressure slot,
                  // so the per-species sources have to be accumulated instead
                  // of overwriting each other.  A single shared slot can only
                  // carry one ion temperature, so the newly created ions end
                  // up with the number-density weighted mean source
                  // temperature.
                  arr(i, j, k, iUx_I[0]) = sumRhoUx * rhoNormPerT;
                  arr(i, j, k, iUy_I[0]) = sumRhoUy * rhoNormPerT;
                  arr(i, j, k, iUz_I[0]) = sumRhoUz * rhoNormPerT;
                  arr(i, j, k, iP_I[0]) =
                      sumP * ion_pressure_factor() * pNormPerT;
                }

                // Electron pressure source.  Only written when the MHD state
                // carries a separate electron pressure variable: iPe == 0 is
                // the BATSRUS sentinel for "no Pe variable" (Pe_ = 1) and
                // would alias the total density slot.
                if (useMhdPe && iPe > 0 && srcP[0] > 0)
                  arr(i, j, k, iPe) = srcP[0] * pNormPerT;

                // ---- Loss frequencies [code units of 1/time] ----
                // apply_loss() multiplies every particle weight by
                // exp(-nu*dt), so the stored value must be a frequency, not
                // a mass-density rate: it does not depend on how much plasma
                // the state `other` happens to hold at this node.
                if (doLoss) {
                  // Electron-ion recombinations per m^3 and s, with the ion
                  // losses of one step folded in by the stiff factor.
                  amrex::Real recombSI = 0.0;
                  for (int iSp = 1; iSp <= nIonS; ++iSp) {
                    lossArr(i, j, k, iSp) = nuLoss[iSp] / get_Si2NoT();
                    if (nuRecomb[iSp] > 0.0) {
                      const amrex::Real ni =
                          species_number_density_si(other, mfi, idx, iLev, iSp);
                      recombSI += ni * nuRecomb[iSp] * stiff[iSp];
                    }
                  }
                  // A kinetic electron species loses as many electrons as
                  // the ions recombine with, i.e. the frequency nu_e with
                  // n_e*(1 - exp(-nu_e*dt)) = recombSI*dt.
                  if (has_electron_species() && recombSI > 0.0) {
                    fetch_electron_plasma();
                    const amrex::Real neSI =
                        ne / (get_Si2NoRho() * cProtonMassSI);
                    amrex::Real nuE = 0.0;
                    if (neSI > 0.0) {
                      if (dtSI > 0.0) {
                        const amrex::Real frac =
                            std::min(recombSI * dtSI / neSI, 1.0 - 1e-12);
                        nuE = -std::log1p(-frac) / dtSI;
                      } else {
                        nuE = recombSI / neSI;
                      }
                    }
                    lossArr(i, j, k, 0) = nuE / get_Si2NoT();
                  }
                }
              } // for k
        }
      }
    }

    if (!grid.is_grid_empty() && useCurrent) {
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        amrex::MultiFab currentMF(nodeFluid[iLev], amrex::make_alias, iJx, 3);
        currentMF.setVal(0, currentMF.nGrow());
      }
    }
  }
};

#endif

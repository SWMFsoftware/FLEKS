#include <cmath>

#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include "DemoConfig.h"

DemoConfig::DemoConfig() {
  amrex::ParmParse pp;
  pp.query("mode", mode);
  pp.query("initial_condition", initialCondition);
  pp.query("n_cell", nCell);
  pp.query("max_grid_size", maxGridSize);
  pp.query("n_iter", nIter);
  pp.query("n_step", nStep);
  pp.query("particles_per_cell", particlesPerCell);
  pp.query("seed", seed);
  pp.query("diagnostic_interval", diagnosticInterval);
  pp.query("gmres_restart", gmresRestart);
  pp.query("gmres_max_iter", gmresMaxIter);
  pp.query("div_e_max_iter", divEMaxIter);
  pp.query("n_div_e_correction", nDivECorrection);
  pp.query("production_ghost", productionGhost);
  pp.query("sort_interval", sortInterval);
  pp.query("dt", dt);
  pp.query("theta", theta);
  pp.query("mass_ratio", massRatio);
  pp.query("perturbation", perturbation);
  pp.query("b0", b0);
  pp.query("solver_tolerance", solverTolerance);
  pp.query("div_e_tolerance", divETolerance);
  pp.query("coef_diff", coefDiff);
  pp.query("correct_div_e", correctDivE);
  pp.query("fused_push", fusedPush);
  pp.query("buffered_response", bufferedResponse);
  pp.query("clustered", clustered);
  pp.query("reference_prefix", referencePrefix);
  if (mode != "fields" && mode != "particles" && mode != "deposit" &&
      mode != "implicit" && mode != "benchmark" && mode != "all")
    amrex::Abort(
        "mode must be fields, particles, deposit, implicit, benchmark, or all");
  if (nCell < 4 || maxGridSize < 2 || nIter < 1 || nStep < 0 ||
      particlesPerCell < 1 || diagnosticInterval < 1 || gmresRestart < 1 ||
      gmresMaxIter < 1 || divEMaxIter < 1 || nDivECorrection < 1 ||
      productionGhost < 1 || sortInterval < 0)
    amrex::Abort("Invalid grid, particle, iteration, or diagnostic count");
  if (!std::isfinite(dt) || dt <= 0 || !std::isfinite(theta) || theta < 0.5 ||
      theta > 1 || !std::isfinite(massRatio) || massRatio <= 0 ||
      !std::isfinite(perturbation) || std::abs(perturbation) > 0.05 ||
      !std::isfinite(b0) || !std::isfinite(solverTolerance) ||
      solverTolerance <= 0 || solverTolerance >= 1)
    amrex::Abort("Require dt>0, 0.5<=theta<=1, mass_ratio>0, "
                 "|perturbation|<=0.05, and 0<solver_tolerance<1");
  if (initialCondition != "uniform" && initialCondition != "langmuir")
    amrex::Abort("initial_condition must be uniform or langmuir");
  if (!std::isfinite(coefDiff) || coefDiff < 0 || coefDiff > 1 ||
      !std::isfinite(divETolerance) || divETolerance <= 0 || divETolerance >= 1)
    amrex::Abort("Require 0<=coef_diff<=1 and 0<div_e_tolerance<1");
}

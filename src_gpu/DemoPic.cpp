#include <algorithm>
#include <cmath>
#include <iomanip>

#include <AMReX_Gpu.H>
#include <AMReX_ParReduce.H>
#include <AMReX_Print.H>

#include "DemoChecks.h"
#include "DemoPic.h"

using namespace amrex;

DemoPic::DemoPic(const DemoConfig& configIn)
    : config(configIn), grid(config), fields(grid, config) {}

int DemoPic::run() {
  Print() << "GPU PIC demo: mode=" << config.mode << ", cells=" << config.nCell
          << "^3, boxes=" << grid.ba.size()
          << ", ranks=" << ParallelDescriptor::NProcs() << "\n";
  if (config.mode == "fields" || config.mode == "all") {
    DemoChecks baseline(grid, config);
    failures += baseline.run();
    failures += fields.run_checks();
  }
  if (config.mode == "particles" || config.mode == "deposit" ||
      config.mode == "all") {
    DemoParticles positive(grid, config, 1, 1, 0);
    failures += positive.run_checks(fields);
    DemoParticles negative(grid, config, -1, 1, 1);
    failures += negative.run_checks(fields);
  }
  if (failures == 0 && (config.mode == "implicit" ||
                        config.mode == "benchmark" || config.mode == "all")) {
    initialize();
    advance();
    if (config.mode == "benchmark")
      benchmark();
  }
  return failures;
}

void DemoPic::initialize() {
  fields.nodeE.setVal(0);
  fields.nodeEth.setVal(0);
  fields.centerB.setVal(0);
  fields.centerB.setVal(config.b0, 2, 1, fields.centerB.nGrow());
  fields.refresh_node_b();
  fields.rhoOld.setVal(0);
  species.clear();
  species.push_back(std::make_unique<DemoParticles>(grid, config, -1, 1, 0));
  species.push_back(
      std::make_unique<DemoParticles>(grid, config, 1, config.massRatio, 1));
  if (config.initialCondition == "uniform")
    species[0]->initialize(0, { config.perturbation, 0, 0 });
  else
    species[0]->initialize(config.perturbation);
  species[1]->initialize();
  // FLEKS stage convention: x^(n+1/2), v^n, E^n and B^n. Start with a
  // half drift from x^0 to x^(1/2); the uniform mode remains exactly neutral.
  for (auto& particles : species)
    particles->drift_positions(0.5 * config.dt);
  deposit_charge(false, config.deposit);
  MultiFab::Copy(fields.rhoOld, fields.rho, 0, 0, 1, 1);
  deposit_moments(config.deposit);
  expectedCount = 2 * Long(config.nCell) * config.nCell * config.nCell *
                  config.particlesPerCell;
  Print() << "STAGES x^(n+1/2),v^n,E^n,B^n -> response -> E^(n+theta) -> "
             "v^(n+1),x^(n+3/2) -> B^(n+1) -> position correction\n"
          << "PHYSICS initial_condition=" << config.initialCondition
          << " theta=" << config.theta << " omega_p_dt="
          << std::sqrt(demo_field::fourPi * (1 + 1 / config.massRatio)) *
                 config.dt
          << " correction=" << config.correctDivE
          << " rho_theta=" << demo_field::rhoTheta << "\n";
}

void DemoPic::deposit_moments(const std::string& method) {
  fields.moments.setVal(0);
  for (auto& particles : species) {
    if (method == "cell")
      particles->deposit_cell(fields.moments);
    else
      particles->deposit(fields.moments);
  }
  fields.moments.SumBoundary(grid.geom.periodicity());
  fields.moments.FillBoundary(grid.geom.periodicity());
}

void DemoPic::deposit_charge(bool matrix, const std::string& method) {
  fields.rho.setVal(0);
  if (matrix)
    fields.centerMM.setVal(0);
  MultiFab* centerMatrix = matrix ? &fields.centerMM : nullptr;
  for (auto& particles : species) {
    if (method == "cell")
      particles->deposit_center_charge_cell(fields.rho, centerMatrix);
    else
      particles->deposit_center_charge(fields.rho, centerMatrix);
  }
  fields.rho.SumBoundary(grid.geom.periodicity());
  fields.rho.FillBoundary(grid.geom.periodicity());
  if (matrix) {
    fields.centerMM.SumBoundary(grid.geom.periodicity());
    fields.centerMM.FillBoundary(grid.geom.periodicity());
  }
}

void DemoPic::accumulate(DemoParticles& particles, const std::string& method) {
  if (method == "cell")
    particles.accumulate_response_cell(fields, config.dt);
  else if (method == "buffered")
    particles.accumulate_response_buffered(fields, config.dt);
  else
    particles.accumulate_response(fields, config.dt);
}

void DemoPic::measure(const std::string& label,
                      const std::function<void()>& operation) {
  BL_PROFILE(label);
  if (!profileActive) {
    operation();
    return;
  }
  Gpu::synchronize();
  const double start = second();
  operation();
  Gpu::synchronize();
  double elapsed = second() - start;
  ParallelDescriptor::ReduceRealMax(elapsed);
  timings[label] += elapsed;
}

int DemoPic::advance() {
  report(0);
  const bool throughput =
      config.mode == "benchmark" && config.benchmarkTiming == "throughput";
  bool windowActive = false;
  double windowStart = 0, diagnosticElapsed = 0;
  for (int step = 1; step <= config.nStep && failures == 0; ++step) {
    const bool measured =
        config.mode == "benchmark" && step > config.warmupSteps;
    profileActive = measured && !throughput;
    // Solver counters cover exactly the measured steps; section timers run
    // only in stage mode with solver_sections=1.
    if (measured && step == config.warmupSteps + 1) {
      fields.electric_solver().reset_stats();
      fields.div_e_solver().reset_stats();
    }
    fields.electric_solver().set_timing(profileActive && config.solverSections);
    fields.div_e_solver().set_timing(profileActive && config.solverSections);
    if (throughput && measured && !windowActive) {
      Gpu::synchronize();
      windowStart = second();
      windowActive = true;
    }
    const auto savedTimings =
        profileActive ? timings : std::map<std::string, double>();
    if (profileActive)
      Gpu::synchronize();
    const double start = second();
    measure("matrix_assembly", [&] {
      fields.clear_response();
      for (auto& particles : species)
        accumulate(*particles, config.response);
    });
    measure("matrix_exchange", [&] { fields.finish_response(); });
    bool solved = false;
    measure("electric_solve", [&] { solved = fields.solve_e(); });
    if (!solved) {
      Print() << "FAIL electric solve: true L2 residual=" << fields.lastResidual
              << "\n";
      ++failures;
      if (profileActive)
        timings = savedTimings;
      break;
    }
    if (measured)
      electricIterations += fields.lastIterations;
    if (!config.fusedPush)
      measure("particle_gather", [&] {
        for (auto& particles : species)
          particles->gather(fields.nodeEth, fields.nodeB);
      });
    measure(config.fusedPush ? "fused_push_migration"
                             : "particle_push_migration",
            [&] {
              for (auto& particles : species)
                particles->push(config.dt, config.fusedPush, &fields.nodeEth,
                                &fields.nodeB);
            });
    measure("magnetic_advance", [&] { fields.advance_b(); });
    if (config.correctDivE) {
      for (int iteration = 0; iteration < config.nDivECorrection; ++iteration) {
        measure("center_deposit",
                [&] { deposit_charge(true, config.deposit); });
        measure("div_e_solve", [&] { solved = fields.solve_div_e(); });
        if (measured) {
          correctionIterations += fields.lastDivIterations;
          correctionWarnings += !solved;
        }
        // Production treats the correction as best effort: LinearSolver
        // stops at its iteration limit and the partial potential is applied.
        if (!solved)
          Print() << "WARN div(E) solve reached div_e_max_iter="
                  << config.divEMaxIter
                  << ", true L2 residual=" << fields.lastDivResidual << "\n";
        measure("position_correction_migration", [&] {
          for (auto& particles : species)
            particles->correct_positions(fields.phi);
        });
      }
    }
    measure("final_charge_deposit",
            [&] { deposit_charge(false, config.deposit); });
    measure("moment_deposit", [&] { deposit_moments(config.deposit); });
    if (config.sortInterval > 0 && step % config.sortInterval == 0)
      measure("particle_sort", [&] {
        for (auto& particles : species)
          particles->sort_by_cell();
      });
    if (profileActive)
      Gpu::synchronize();
    const double evolutionElapsed = second() - start;
    if (step % config.diagnosticInterval == 0 || step == config.nStep ||
        failures) {
      const double diagnosticStart = second();
      report(step);
      if (throughput && measured)
        diagnosticElapsed += second() - diagnosticStart;
    }
    // The Gauss diagnostic above uses (1-rhoTheta)*rhoOld+rhoTheta*rhoNew.
    // Preserve the fresh half-stage charge only after evaluating it.
    const double historyStart = second();
    measure("charge_history",
            [&] { MultiFab::Copy(fields.rhoOld, fields.rho, 0, 0, 1, 1); });
    if (profileActive) {
      Gpu::synchronize();
      double elapsed = evolutionElapsed + second() - historyStart;
      ParallelDescriptor::ReduceRealMax(elapsed);
      timings["timestep"] += elapsed;
    }
    ++completedSteps;
    if (measured)
      ++measuredSteps;
  }
  if (windowActive) {
    Gpu::synchronize();
    // One collective at the end of the entire window; diagnostic time is
    // included in total time, and its maximum is reported separately.
    double elapsed[2] = { second() - windowStart, diagnosticElapsed };
    ParallelDescriptor::ReduceRealMax(elapsed, 2);
    timings["timestep"] = elapsed[0];
    timings["diagnostics_included"] = elapsed[1];
  }
  return failures;
}

void DemoPic::report(int step) {
  Long count = 0;
  Real kinetic = 0, charge = 0, current = 0;
  for (const auto& particles : species) {
    const auto diagnostic = particles->diagnostics();
    count += diagnostic.count;
    kinetic += diagnostic.kineticEnergy;
    charge += diagnostic.charge;
    current += diagnostic.current[0];
  }
  const Real energy = fields.field_energy();
  const Real gauss = fields.gauss_error();
  const Real divB = fields.div_b_error();
  // A point sample retains the sign of the uniform electric oscillation.
  // nodeE.min is globally reduced; for nonuniform runs electric_x is a minimum.
  const Real electric = fields.nodeE.min(0);
  Print() << std::setprecision(17) << "DIAG step=" << step
          << " time=" << step * config.dt << " kinetic=" << kinetic
          << " field=" << energy << " charge=" << charge << " gauss=" << gauss
          << " div_b=" << divB << " residual=" << fields.lastResidual
          << " particles=" << count << " iterations=" << fields.lastIterations
          << " div_residual=" << fields.lastDivResidual
          << " div_iterations=" << fields.lastDivIterations
          << " correction=" << config.correctDivE << " theta=" << config.theta
          << " omega_p_dt="
          << std::sqrt(demo_field::fourPi * (1 + 1 / config.massRatio)) *
                 config.dt
          << " electric_x=" << electric << " current_x=" << current << "\n";
  if (count != expectedCount ||
      !std::isfinite(energy + kinetic + gauss + divB) ||
      std::abs(charge) > 1.e-9) {
    Print() << "FAIL nonfinite diagnostic or periodic particle/charge "
               "conservation\n";
    ++failures;
  }
}

void DemoPic::benchmark() {
  const int measured = measuredSteps;
  Print() << "BENCHMARK warmup_steps=" << config.warmupSteps
          << " measured_steps=" << measured
          << " completed_steps=" << completedSteps
          << " timing=" << config.benchmarkTiming << " diagnostics_included="
          << (config.benchmarkTiming == "throughput")
          << " aborted=" << (failures != 0)
          << " electric_iterations=" << electricIterations
          << " correction_iterations=" << correctionIterations
          << " correction_warnings=" << correctionWarnings
          << " fused=" << config.fusedPush
          << " sort_interval=" << config.sortInterval
          << " response=" << config.response << " deposit=" << config.deposit
          << " clustered=" << config.clustered
          << " pair_symmetry=" << config.pairSymmetry
          << " fold_charge=" << config.foldCharge
          << " cell_chunk_size=" << config.cellChunkSize
          << " gmres_batched=" << config.batchedGmres
          << " div_e_method=" << config.divEMethod
          << " cell_precision=" << config.cellPrecision << "\n";
  for (const auto& entry : timings)
    Print() << "TIMING stage=" << entry.first
            << " seconds=" << entry.second / std::max(1, measured) << "\n";
  // Per-solver costs. Section times are zero unless solver_sections=1 in
  // stage mode.
  for (const auto& entry :
       { std::make_pair("electric", &fields.electric_solver()),
         std::make_pair("div_e", &fields.div_e_solver()) }) {
    const auto& s = entry.second->stats();
    const double perStep = 1. / std::max(1, measured);
    const double perIteration = 1. / std::max(1L, s.iterations);
    Print() << "SOLVER name=" << entry.first << " solves=" << s.solves
            << " iterations=" << s.iterations
            << " host_waits_per_iteration=" << s.hostWaits * perIteration
            << " norm_fallbacks=" << s.normFallbacks
            << " operator_ms_per_step=" << 1e3 * s.operatorSeconds * perStep
            << " orthogonalize_ms_per_step="
            << 1e3 * s.orthogonalizeSeconds * perStep
            << " update_ms_per_step=" << 1e3 * s.updateSeconds * perStep
            << " residual_ms_per_step=" << 1e3 * s.residualSeconds * perStep
            << "\n";
  }
  if (failures == 0 && measured && timings["timestep"] > 0)
    Print() << "THROUGHPUT particles_per_second="
            << expectedCount * measured / timings["timestep"]
            << " cells_per_second="
            << double(config.nCell) * config.nCell * config.nCell * measured /
                   timings["timestep"]
            << "\n";

  if (!config.microbenchmarks || failures)
    return;

  // Compare communication volume for operator-specific and production widths.
  for (int width : { 1, config.productionGhost }) {
    MultiFab test(grid.nodeBa, grid.dm, 3, width);
    test.setVal(1);
    test.FillBoundary(grid.geom.periodicity());
    Gpu::synchronize();
    const double start = second();
    for (int iteration = 0; iteration < config.nIter; ++iteration)
      test.FillBoundary(grid.geom.periodicity());
    Gpu::synchronize();
    double elapsed = (second() - start) / config.nIter;
    ParallelDescriptor::ReduceRealMax(elapsed);
    Print() << "GHOST width=" << width << " seconds=" << elapsed << "\n";
  }

  // Compare equal 243-real response data in the production-style custom
  // FabArray and a component MultiFab. This isolates matrix-vector traffic.
  MultiFab components(grid.nodeBa, grid.dm, 243, 1);
  MultiFab customResult(grid.nodeBa, grid.dm, 3, 0);
  MultiFab componentResult(grid.nodeBa, grid.dm, 3, 0);
  for (auto& particles : species)
    failures += particles->benchmark_layouts(fields.nodeEth, fields.nodeB);

  // Timing includes zeroing the response arrays. The cell variant also
  // includes its counting sort; the buffered variant includes a particle sort,
  // cell-range construction and temporary buffer setup. Each alternative is
  // compared entry by entry with the atomic reference after the exchange.
  auto pack_response = [&] {
    for (MFIter mfi(components); mfi.isValid(); ++mfi) {
      const auto packed = components.array(mfi);
      const auto custom = fields.nodeMM.const_array(mfi);
      ParallelFor(mfi.fabbox(), 243,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int c)
                      noexcept { packed(i, j, k, c) = custom(i, j, k)[c]; });
    }
  };
  for (const std::string method : { "atomic", "cell", "buffered" }) {
    auto assemble = [&] {
      fields.clear_response();
      for (auto& particles : species)
        accumulate(*particles, method);
    };
    assemble();
    Gpu::synchronize();
    const double start = second();
    for (int iteration = 0; iteration < config.nIter; ++iteration)
      assemble();
    Gpu::synchronize();
    double elapsed = (second() - start) / config.nIter;
    ParallelDescriptor::ReduceRealMax(elapsed);
    Print() << "RESPONSE layout=" << method << " seconds=" << elapsed << "\n";
    fields.finish_response();
    if (method == "atomic") {
      pack_response();
      continue;
    }
    const auto packedArrays = components.const_arrays();
    const auto customArrays = fields.nodeMM.const_arrays();
    Real responseError = ParReduce(
        TypeList<ReduceOpMax>{}, TypeList<Real>{}, components, IntVect(0),
        [=] AMREX_GPU_DEVICE(int box, int i, int j, int k) -> GpuTuple<Real> {
          Real error = 0;
          for (int c = 0; c < 243; ++c)
            error = amrex::max(error, std::abs(customArrays[box](i, j, k)[c] -
                                               packedArrays[box](i, j, k, c)));
          return { error };
        });
    ParallelDescriptor::ReduceRealMax(responseError);
    // cell_precision=single sums each cell in float: allow 1e-6 of the
    // largest response entry for the cell method.
    const Real responseTolerance =
        method == "cell" && config.cellPrecision == "single"
            ? Real(1.e-6) * components.norm0(0, 243, IntVect(0))
            : Real(1.e-10);
    Print() << (responseError < responseTolerance ? "PASS" : "FAIL") << " "
            << method << " response matrix error=" << responseError
            << " (tol=" << responseTolerance << ")\n";
    if (!(responseError < responseTolerance))
      ++failures;
  }

  // Particle deposits: time each method for all species, including the
  // counting sort and the boundary exchange.
  MultiFab rhoAtomic(grid.ba, grid.dm, 1, 0), mmAtomic(grid.ba, grid.dm, 27, 0),
      momentsAtomic(grid.nodeBa, grid.dm, 10, 0);
  for (const std::string method : { "atomic", "cell" }) {
    auto time_operation = [&](const std::function<void()>& operation) {
      operation();
      Gpu::synchronize();
      const double start = second();
      for (int iteration = 0; iteration < config.nIter; ++iteration)
        operation();
      Gpu::synchronize();
      double elapsed = (second() - start) / config.nIter;
      ParallelDescriptor::ReduceRealMax(elapsed);
      return elapsed;
    };
    const double chargeOnly =
        time_operation([&] { deposit_charge(false, method); });
    const double chargeMatrix =
        time_operation([&] { deposit_charge(true, method); });
    const double moments = time_operation([&] { deposit_moments(method); });
    Print() << "DEPOSIT layout=" << method << " charge_seconds=" << chargeOnly
            << " charge_matrix_seconds=" << chargeMatrix
            << " moment_seconds=" << moments << "\n";
  }
  // Compare per species: a neutral plasma nearly cancels the summed,
  // charge-weighted fields, which would turn roundoff into large relative
  // errors. Within one species every contribution has the same charge sign.
  // Each error is relative to that field's maximum magnitude.
  auto relative_error = [](MultiFab& reference, const MultiFab& value,
                           int ncomp) {
    const Real scale =
        std::max(reference.norm0(0, ncomp, IntVect(0)), Real(1.e-300));
    MultiFab::Subtract(reference, value, 0, 0, ncomp, 0);
    return reference.norm0(0, ncomp, IntVect(0)) / scale;
  };
  const auto period = grid.geom.periodicity();
  Real rhoError = 0, matrixError = 0, momentError = 0;
  for (auto& particles : species) {
    for (const std::string method : { "atomic", "cell" }) {
      fields.rho.setVal(0);
      fields.centerMM.setVal(0);
      fields.moments.setVal(0);
      if (method == "cell") {
        particles->deposit_center_charge_cell(fields.rho, &fields.centerMM);
        particles->deposit_cell(fields.moments);
      } else {
        particles->deposit_center_charge(fields.rho, &fields.centerMM);
        particles->deposit(fields.moments);
      }
      fields.rho.SumBoundary(period);
      fields.centerMM.SumBoundary(period);
      fields.moments.SumBoundary(period);
      if (method == "atomic") {
        MultiFab::Copy(rhoAtomic, fields.rho, 0, 0, 1, 0);
        MultiFab::Copy(mmAtomic, fields.centerMM, 0, 0, 27, 0);
        MultiFab::Copy(momentsAtomic, fields.moments, 0, 0, 10, 0);
        continue;
      }
      rhoError = std::max(rhoError, relative_error(rhoAtomic, fields.rho, 1));
      matrixError =
          std::max(matrixError, relative_error(mmAtomic, fields.centerMM, 27));
      momentError = std::max(momentError,
                             relative_error(momentsAtomic, fields.moments, 10));
    }
  }
  // Charge and moments stay double precision; the correction matrix uses
  // float per-cell sums with cell_precision=single.
  const Real matrixTolerance =
      config.cellPrecision == "single" ? Real(1.e-6) : Real(1.e-12);
  const bool depositPass = std::max(rhoError, momentError) < 1.e-12 &&
                           matrixError < matrixTolerance;
  Print() << (depositPass ? "PASS" : "FAIL")
          << " cell deposit per-species relative error (rho=" << rhoError
          << " matrix=" << matrixError << " moments=" << momentError
          << ", matrix tol=" << matrixTolerance << ")\n";
  if (!depositPass)
    ++failures;
  for (MFIter mfi(components); mfi.isValid(); ++mfi) {
    const auto packed = components.array(mfi);
    const auto custom = fields.nodeMM.const_array(mfi);
    ParallelFor(mfi.fabbox(), 243,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c)
                    noexcept { packed(i, j, k, c) = custom(i, j, k)[c]; });
  }
  fields.nodeE.FillBoundary(grid.geom.periodicity());
  auto apply = [&](MultiFab& output, bool componentStorage) {
    for (MFIter mfi(output); mfi.isValid(); ++mfi) {
      const auto result = output.array(mfi);
      const auto packed = components.const_array(mfi);
      const auto custom = fields.nodeMM.const_array(mfi);
      const auto input = fields.nodeE.const_array(mfi);
      ParallelFor(mfi.validbox(), 3,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int row) noexcept {
                    Real value = 0;
                    for (int dz = -1; dz <= 1; ++dz)
                      for (int dy = -1; dy <= 1; ++dy)
                        for (int dx = -1; dx <= 1; ++dx)
                          for (int column = 0; column < 3; ++column) {
                            const int index =
                                ((dz + 1) * 9 + (dy + 1) * 3 + dx + 1) * 9 +
                                row * 3 + column;
                            const Real coefficient =
                                componentStorage ? packed(i, j, k, index)
                                                 : custom(i, j, k)[index];
                            value += coefficient *
                                     input(i + dx, j + dy, k + dz, column);
                          }
                    result(i, j, k, row) = value;
                  });
    }
  };
  for (bool componentStorage : { false, true }) {
    auto& result = componentStorage ? componentResult : customResult;
    apply(result, componentStorage);
    Gpu::synchronize();
    const double start = second();
    for (int iteration = 0; iteration < config.nIter; ++iteration)
      apply(result, componentStorage);
    Gpu::synchronize();
    double elapsed = (second() - start) / config.nIter;
    ParallelDescriptor::ReduceRealMax(elapsed);
    Print() << "MATRIX_STORAGE layout="
            << (componentStorage ? "components" : "custom")
            << " seconds=" << elapsed << "\n";
  }
  MultiFab::Subtract(customResult, componentResult, 0, 0, 3, 0);
  const Real error = customResult.norm0(0, 3, IntVect(0));
  Print() << (error < 1.e-12 ? "PASS" : "FAIL")
          << " matrix storage comparison error=" << error << "\n";
  if (!(error < 1.e-12))
    ++failures;
}

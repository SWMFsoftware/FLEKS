#include <algorithm>
#include <cmath>
#include <limits>

#include <AMReX_Gpu.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include "DemoConfig.h"
#include "DemoFields.h"
#include "DemoGrid.h"
#include "DemoParticleKernels.h"
#include "DemoParticles.h"

using namespace amrex;
namespace {
constexpr Real pi = 3.1415926535897932384626433832795;

Real cell_volume(const DemoGrid& grid) {
  const auto dx = grid.geom.CellSizeArray();
  return dx[0] * dx[1] * dx[2];
}
} // namespace

DemoParticles::DemoParticles(const DemoGrid& gridIn, const DemoConfig& configIn,
                             Real chargeIn, Real massIn, int speciesIn)
    : Base(Vector<Geometry>{ gridIn.geom },
           Vector<DistributionMapping>{ gridIn.dm },
           Vector<BoxArray>{ gridIn.ba }, Vector<int>{}),
      grid(gridIn),
      config(configIn),
      charge(chargeIn),
      mass(massIn),
      species(speciesIn) {
  AMREX_ALWAYS_ASSERT(mass > 0 && charge != 0);
}

void DemoParticles::initialize(Real perturbation, GpuArray<Real, 3> drift) {
  clearParticles();
  gatheredValid = false;
  const auto dx = grid.geom.CellSizeArray();
  const auto lo = grid.geom.ProbLoArray();
  const int n = config.nCell;
  const int ppc = config.particlesPerCell;
  const int spec = species;
  const bool clustered = config.clustered;
  const int sampleOffset = (config.seed % ppc + ppc) % ppc;
  const int rank = ParallelDescriptor::MyProc();
  const Real qp = charge * cell_volume(grid) / ppc;
  int lattice = static_cast<int>(std::round(std::cbrt(Real(ppc))));
  const bool cube = lattice * lattice * lattice == ppc;
  const Long speciesOffset = Long(spec) * Long(n) * n * n * ppc;
  MultiFab layout(grid.ba, grid.dm, 1, 0);
  for (MFIter mfi(layout); mfi.isValid(); ++mfi) {
    const auto box = mfi.validbox();
    const auto small = box.smallEnd();
    const auto len = box.length();
    const Long count = box.numPts() * ppc;
    auto& tile = DefineAndReturnParticleTile(0, mfi.index(), 0);
    tile.resize(count);
    auto* particles = tile.GetArrayOfStructs().data();
    ParallelFor(count, [=] AMREX_GPU_DEVICE(Long ip) noexcept {
      const int sample = (ip % ppc + sampleOffset) % ppc;
      Long cell = ip / ppc;
      const int i = small[0] + cell % len[0];
      cell /= len[0];
      const int j = small[1] + cell % len[1];
      const int k = small[2] + cell / len[1];
      auto& p = particles[ip];
      p.id() = speciesOffset + ((Long(k) * n + j) * n + i) * ppc + sample + 1;
      p.cpu() = rank;
      p.idata(0) = spec;
      Real fraction[3];
      if (cube) {
        fraction[0] = (sample % lattice + 0.5) / lattice;
        fraction[1] = ((sample / lattice) % lattice + 0.5) / lattice;
        fraction[2] = (sample / (lattice * lattice) + 0.5) / lattice;
      } else {
        fraction[0] = (sample + 0.5) / ppc;
        fraction[1] = (sample + 0.5) * 0.6180339887498948482;
        fraction[2] = (sample + 0.5) * 0.4142135623730950488;
        fraction[1] -= std::floor(fraction[1]);
        fraction[2] -= std::floor(fraction[2]);
      }
      const int index[3] = { i, j, k };
      for (int d = 0; d < 3; ++d) {
        p.pos(d) = lo[d] + (index[d] + fraction[d]) * dx[d];
        if (clustered)
          p.pos(d) = lo[d] + 0.5 + 0.2 * (p.pos(d) - lo[d] - 0.5);
        p.rdata(d) = drift[d];
      }
      // A velocity perturbation excites the longitudinal plasma mode without
      // changing the initially neutral electron/ion charge distribution.
      p.rdata(0) += perturbation * std::cos(2 * pi * p.pos(0));
      p.rdata(iqp_) = qp;
    });
  }
  Redistribute();
}

void DemoParticles::drift_positions(Real dt) {
  const auto lo = grid.geom.ProbLoArray();
  const auto hi = grid.geom.ProbHiArray();
  for (ParIterType pti(*this, 0); pti.isValid(); ++pti) {
    const int np = pti.numParticles();
    auto* particles = pti.GetArrayOfStructs().data();
    ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
      auto& p = particles[ip];
      for (int d = 0; d < 3; ++d) {
        const Real length = hi[d] - lo[d];
        p.pos(d) += dt * p.rdata(d);
        p.pos(d) -= length * std::floor((p.pos(d) - lo[d]) / length);
      }
    });
  }
  gatheredValid = false;
  Redistribute();
}

void DemoParticles::gather(const MultiFab& e, const MultiFab& b) {
  const auto lo = grid.geom.ProbLoArray();
  const auto invDx = grid.geom.InvCellSizeArray();
  for (ParIterType pti(*this, 0); pti.isValid(); ++pti) {
    const int np = pti.numParticles();
    const auto* particles = pti.GetArrayOfStructs().data();
    auto& values = gathered[{ pti.index(), pti.LocalTileIndex() }];
    values.resize(6 * np);
    auto* data = values.data();
    const auto ea = e.const_array(pti);
    const auto ba = b.const_array(pti);
    ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
      const auto s = demo_particle::shape(particles[ip], lo, invDx);
      Real ep[3], bp[3];
      demo_particle::gather(ea, s, ep);
      demo_particle::gather(ba, s, bp);
      for (int d = 0; d < 3; ++d) {
        data[6 * ip + d] = ep[d];
        data[6 * ip + d + 3] = bp[d];
      }
    });
  }
  gatheredValid = true;
}

void DemoParticles::push(Real dt, bool fused, const MultiFab* e,
                         const MultiFab* b, bool averageDrift) {
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fused ? (e && b) : gatheredValid,
                                   "Call gather before a separate push");
  const auto lo = grid.geom.ProbLoArray();
  const auto hi = grid.geom.ProbHiArray();
  const auto invDx = grid.geom.InvCellSizeArray();
  const GpuArray<Real, 3> length{ hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2] };
  const Real qom = charge / mass;
  for (ParIterType pti(*this, 0); pti.isValid(); ++pti) {
    const int np = pti.numParticles();
    auto* particles = pti.GetArrayOfStructs().data();
    const Real* data =
        fused ? nullptr
              : gathered.at({ pti.index(), pti.LocalTileIndex() }).data();
    const auto ea = fused ? e->const_array(pti) : Array4<Real const>();
    const auto ba = fused ? b->const_array(pti) : Array4<Real const>();
    ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
      auto& p = particles[ip];
      if (p.id() <= 0)
        return;
      Real ep[3], bp[3];
      if (fused) {
        const auto s = demo_particle::shape(p, lo, invDx);
        demo_particle::gather(ea, s, ep);
        demo_particle::gather(ba, s, bp);
      } else {
        for (int d = 0; d < 3; ++d) {
          ep[d] = data[6 * ip + d];
          bp[d] = data[6 * ip + d + 3];
        }
      }
      demo_particle::push(p, ep, bp, qom, dt, lo, length, averageDrift);
    });
  }
  gatheredValid = false;
  Redistribute();
}

void DemoParticles::deposit(MultiFab& moments) {
  AMREX_ALWAYS_ASSERT(moments.nComp() >= 10);
  const auto lo = grid.geom.ProbLoArray();
  const auto invDx = grid.geom.InvCellSizeArray();
  const Real invVol = 1 / cell_volume(grid);
  for (ParConstIterType pti(*this, 0); pti.isValid(); ++pti) {
    const int np = pti.numParticles();
    const auto* particles = pti.GetArrayOfStructs().data();
    const auto out = moments.array(pti);
    ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
      const auto& p = particles[ip];
      if (p.id() <= 0)
        return;
      const auto s = demo_particle::shape(p, lo, invDx);
      const Real u = p.rdata(0), v = p.rdata(1), w = p.rdata(2);
      const Real q = p.rdata(3) * invVol;
      const Real value[10] = { q,         q * u,     q * v,     q * w,
                               q * u * u, q * v * v, q * w * w, q * u * v,
                               q * u * w, q * v * w };
      for (int a = 0; a < 8; ++a) {
        int i, j, k;
        demo_particle::index(s, a, i, j, k);
        for (int c = 0; c < 10; ++c)
          Gpu::Atomic::Add(&out(i, j, k, c), s.weight[a] * value[c]);
      }
    });
  }
}

void DemoParticles::accumulate_response(DemoFields& fields, Real dt) {
  const auto lo = grid.geom.ProbLoArray();
  const auto invDx = grid.geom.InvCellSizeArray();
  const Real qdt2m = charge / mass * dt / 2;
  const Real invVol = 1 / cell_volume(grid);
  for (ParConstIterType pti(*this, 0); pti.isValid(); ++pti) {
    const int np = pti.numParticles();
    const auto* particles = pti.GetArrayOfStructs().data();
    const auto ba = fields.nodeB.const_array(pti);
    const auto current = fields.jHat.array(pti);
    const auto mm = fields.nodeMM.array(pti);
    ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
      const auto& p = particles[ip];
      if (p.id() <= 0)
        return;
      const auto s = demo_particle::shape(p, lo, invDx);
      Real b[3], alpha[9], jhat[3] = { 0, 0, 0 };
      demo_particle::gather(ba, s, b);
      demo_particle::alpha(b, qdt2m, alpha);
      for (int d = 0; d < 3; ++d)
        for (int c = 0; c < 3; ++c)
          jhat[d] += p.rdata(3) * alpha[3 * d + c] * p.rdata(c);
      for (int a = 0; a < 8; ++a) {
        int i, j, k;
        demo_particle::index(s, a, i, j, k);
        for (int d = 0; d < 3; ++d)
          Gpu::Atomic::Add(&current(i, j, k, d), s.weight[a] * jhat[d]);
        // Compute all 64 node pairs directly. Production computes one z
        // half and copies its symmetric geometrical weights afterwards.
        for (int c = 0; c < 8; ++c) {
          int ii, jj, kk;
          demo_particle::index(s, c, ii, jj, kk);
          const int offset =
              (kk - k + 1) * 81 + (jj - j + 1) * 27 + (ii - i + 1) * 9;
          const Real factor =
              s.weight[a] * s.weight[c] * p.rdata(3) * qdt2m * invVol;
          for (int d = 0; d < 9; ++d)
            Gpu::Atomic::Add(&mm(i, j, k)[offset + d], factor * alpha[d]);
        }
      }
    });
  }
}

void DemoParticles::deposit_center_charge(MultiFab& rho,
                                          MultiFab* centerMatrix) {
  const auto lo = grid.geom.ProbLoArray();
  const auto invDx = grid.geom.InvCellSizeArray();
  const Real invVol = 1 / cell_volume(grid);
  const bool makeMatrix = centerMatrix != nullptr;
  for (ParConstIterType pti(*this, 0); pti.isValid(); ++pti) {
    const int np = pti.numParticles();
    const auto* particles = pti.GetArrayOfStructs().data();
    const auto r = rho.array(pti);
    const auto mm = makeMatrix ? centerMatrix->array(pti) : Array4<Real>();
    ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
      const auto& p = particles[ip];
      if (p.id() <= 0)
        return;
      const auto s = demo_particle::shape(p, lo, invDx, 0.5);
      for (int a = 0; a < 8; ++a) {
        int i, j, k;
        demo_particle::index(s, a, i, j, k);
        Gpu::Atomic::Add(&r(i, j, k), p.rdata(3) * invVol * s.weight[a]);
        if (makeMatrix)
          for (int c = 0; c < 8; ++c) {
            int ii, jj, kk;
            demo_particle::index(s, c, ii, jj, kk);
            Real dot = 0;
            for (int d = 0; d < 3; ++d)
              dot += s.derivative[a][d] * s.derivative[c][d];
            const int offset =
                (ii - i + 1) * 9 + (jj - j + 1) * 3 + (kk - k + 1);
            Gpu::Atomic::Add(&mm(i, j, k, offset),
                             std::abs(p.rdata(3)) * invVol * dot);
          }
      }
    });
  }
}

void DemoParticles::correct_positions(const MultiFab& phi) {
  const auto lo = grid.geom.ProbLoArray();
  const auto hi = grid.geom.ProbHiArray();
  const auto invDx = grid.geom.InvCellSizeArray();
  const auto dx = grid.geom.CellSizeArray();
  const Real sign = charge > 0 ? 1 : -1;
  for (ParIterType pti(*this, 0); pti.isValid(); ++pti) {
    const int np = pti.numParticles();
    auto* particles = pti.GetArrayOfStructs().data();
    const auto potential = phi.const_array(pti);
    ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
      auto& p = particles[ip];
      if (p.id() <= 0)
        return;
      const auto s = demo_particle::shape(p, lo, invDx, 0.5);
      Real shift[3] = { 0, 0, 0 };
      for (int a = 0; a < 8; ++a) {
        int i, j, k;
        demo_particle::index(s, a, i, j, k);
        for (int d = 0; d < 3; ++d)
          shift[d] += 4 * pi * sign * potential(i, j, k) * s.derivative[a][d];
      }
      Real maximum = 0, norm = 0;
      for (int d = 0; d < 3; ++d) {
        maximum = amrex::max(maximum, std::abs(shift[d] * invDx[d]));
        norm += shift[d] * shift[d];
      }
      const Real factor = maximum > 0.1 ? 0.1 * dx[0] / std::sqrt(norm) : 1;
      for (int d = 0; d < 3; ++d) {
        const Real length = hi[d] - lo[d];
        p.pos(d) += factor * shift[d];
        p.pos(d) -= length * std::floor((p.pos(d) - lo[d]) / length);
      }
    });
  }
  gatheredValid = false;
  Redistribute();
}

void DemoParticles::sort_by_cell() {
  SortParticlesByCell();
  gatheredValid = false;
}

DemoParticleDiagnostics DemoParticles::diagnostics() const {
  ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum,
            ReduceOpSum, ReduceOpSum>
      op;
  ReduceData<Real, Real, Real, Real, Real, Real, Real> data(op);
  const Real massOverCharge = mass / charge;
  for (ParConstIterType pti(*this, 0); pti.isValid(); ++pti) {
    const int np = pti.numParticles();
    const auto* particles = pti.GetArrayOfStructs().data();
    op.eval(np, data,
            [=] AMREX_GPU_DEVICE(
                int ip) -> GpuTuple<Real, Real, Real, Real, Real, Real, Real> {
              const auto& p = particles[ip];
              if (p.id() <= 0)
                return { 0, 0, 0, 0, 0, 0, 0 };
              const Real q = p.rdata(3);
              const Real u = p.rdata(0), v = p.rdata(1), w = p.rdata(2);
              return { 1,
                       q,
                       0.5 * massOverCharge * q * (u * u + v * v + w * w),
                       q * u,
                       q * v,
                       q * w,
                       Real(p.id()) };
            });
  }
  const auto result = data.value(op);
  Real values[7] = { get<0>(result), get<1>(result), get<2>(result),
                     get<3>(result), get<4>(result), get<5>(result),
                     get<6>(result) };
  ParallelDescriptor::ReduceRealSum(values, 7);
  DemoParticleDiagnostics out;
  out.count = static_cast<Long>(values[0]);
  out.charge = values[1];
  out.kineticEnergy = values[2];
  for (int d = 0; d < 3; ++d)
    out.current[d] = values[d + 3];
  out.idSum = values[6];
  return out;
}

std::vector<DemoParticles::ParticleType> DemoParticles::snapshots(
    bool global) const {
  std::vector<ParticleType> out;
  for (ParConstIterType pti(*this, 0); pti.isValid(); ++pti) {
    const auto& aos = pti.GetArrayOfStructs();
    const auto oldSize = out.size();
    out.resize(oldSize + aos.numParticles());
    Gpu::dtoh_memcpy_async(out.data() + oldSize, aos.data(),
                           aos.numParticles() * sizeof(ParticleType));
    Gpu::streamSynchronize();
  }
  Gpu::streamSynchronize();
#ifdef AMREX_USE_MPI
  if (global) {
    const int nproc = ParallelDescriptor::NProcs();
    std::vector<int> counts(nproc), displacements(nproc);
    int bytes = out.size() * sizeof(ParticleType);
    MPI_Allgather(&bytes, 1, MPI_INT, counts.data(), 1, MPI_INT,
                  ParallelDescriptor::Communicator());
    int total = 0;
    for (int rank = 0; rank < nproc; ++rank) {
      displacements[rank] = total;
      total += counts[rank];
    }
    std::vector<ParticleType> all(total / sizeof(ParticleType));
    MPI_Allgatherv(out.data(), bytes, MPI_BYTE, all.data(), counts.data(),
                   displacements.data(), MPI_BYTE,
                   ParallelDescriptor::Communicator());
    out.swap(all);
  }
#else
  amrex::ignore_unused(global);
#endif
  std::sort(out.begin(), out.end(),
            [](const ParticleType& a, const ParticleType& b) {
              return a.id() < b.id();
            });
  return out;
}

void DemoParticles::set_particles(
    const std::vector<ParticleType>& globalParticles) {
  clearParticles();
  gatheredValid = false;
  if (ParallelDescriptor::IOProcessor()) {
    MultiFab layout(grid.ba, grid.dm, 1, 0);
    for (MFIter mfi(layout); mfi.isValid(); ++mfi) {
      auto& tile = DefineAndReturnParticleTile(0, mfi.index(), 0);
      tile.resize(globalParticles.size());
      Gpu::htod_memcpy_async(tile.GetArrayOfStructs().data(),
                             globalParticles.data(),
                             globalParticles.size() * sizeof(ParticleType));
      break;
    }
  }
  Gpu::streamSynchronize();
  Redistribute();
}

namespace {
// These reference routines intentionally do not call the device CIC or Boris
// helpers: their indexing, arithmetic order, and rotation implementation
// differ.
struct ReferenceShape {
  int index[8][3];
  Real weight[8], gradient[8][3];
};
ReferenceShape reference_shape(const DemoParticles::ParticleType& p, int n,
                               bool center = false) {
  ReferenceShape out;
  int cell[3];
  Real f[3];
  for (int d = 0; d < 3; ++d) {
    const Real scaled = p.pos(d) * n - (center ? 0.5 : 0);
    cell[d] = int(std::floor(scaled));
    f[d] = scaled - cell[d];
  }
  int a = 0;
  for (int iz = 0; iz < 2; ++iz)
    for (int iy = 0; iy < 2; ++iy)
      for (int ix = 0; ix < 2; ++ix, ++a) {
        const int corner[3] = { ix, iy, iz };
        const Real w[3] = { ix ? f[0] : 1 - f[0], iy ? f[1] : 1 - f[1],
                            iz ? f[2] : 1 - f[2] };
        out.weight[a] = w[0] * w[1] * w[2];
        for (int d = 0; d < 3; ++d) {
          out.index[a][d] = cell[d] + corner[d];
          out.gradient[a][d] =
              (corner[d] ? n : -n) * w[(d + 1) % 3] * w[(d + 2) % 3];
        }
      }
  return out;
}
std::size_t reference_offset(int i, int j, int k, int n, int nc, int c = 0) {
  i = (i % n + n) % n;
  j = (j % n + n) % n;
  k = (k % n + n) % n;
  return ((std::size_t(k) * n + j) * n + i) * nc + c;
}
Real compare_reference(const MultiFab& mf, const std::vector<Real>& ref,
                       int n) {
  Real error = 0;
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const auto& fab = mf[mfi];
    FArrayBox host(fab.box(), mf.nComp(), The_Pinned_Arena());
    Gpu::dtoh_memcpy_async(host.dataPtr(), fab.dataPtr(), fab.nBytes());
    Gpu::streamSynchronize();
    const auto a = host.const_array();
    amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
      for (int c = 0; c < mf.nComp(); ++c)
        error = std::max(
            error, std::abs(a(i, j, k, c) -
                            ref[reference_offset(i, j, k, n, mf.nComp(), c)]));
    });
  }
  ParallelDescriptor::ReduceRealMax(error);
  return error;
}
Real compare_matrix_reference(const DemoNodeMMFab& mf,
                              const std::vector<Real>& ref, int n) {
  Real error = 0;
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const auto& fab = mf[mfi];
    std::vector<DemoRealMM> host(fab.box().numPts());
    Gpu::dtoh_memcpy_async(host.data(), fab.dataPtr(),
                           fab.box().numPts() * sizeof(DemoRealMM));
    Gpu::streamSynchronize();
    const auto box = fab.box();
    const auto small = box.smallEnd();
    const auto len = box.length();
    amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
      const std::size_t at =
          ((std::size_t(k - small[2]) * len[1] + j - small[1]) * len[0] + i -
           small[0]);
      for (int c = 0; c < 243; ++c)
        error = std::max(
            error,
            std::abs(host[at][c] - ref[reference_offset(i, j, k, n, 243, c)]));
    });
  }
  ParallelDescriptor::ReduceRealMax(error);
  return error;
}
void reference_boris(DemoParticles::ParticleType& p, const Real e[3],
                     const Real b[3], Real qom, Real dt) {
  Real minus[3], prime[3], t[3], s[3], velocity[3];
  Real t2 = 0;
  for (int d = 0; d < 3; ++d) {
    minus[d] = p.rdata(d) + qom * dt * e[d] / 2;
    t[d] = qom * dt * b[d] / 2;
    t2 += t[d] * t[d];
  }
  for (int d = 0; d < 3; ++d) {
    prime[d] = minus[d] + minus[(d + 1) % 3] * t[(d + 2) % 3] -
               minus[(d + 2) % 3] * t[(d + 1) % 3];
    s[d] = 2 * t[d] / (1 + t2);
  }
  for (int d = 0; d < 3; ++d) {
    velocity[d] = minus[d] + prime[(d + 1) % 3] * s[(d + 2) % 3] -
                  prime[(d + 2) % 3] * s[(d + 1) % 3] + qom * dt * e[d] / 2;
    p.pos(d) += dt * velocity[d];
    p.pos(d) -= std::floor(p.pos(d));
    p.rdata(d) = velocity[d];
  }
}
} // namespace

int DemoParticles::run_checks(DemoFields& fields) {
  int failures = 0;
  auto report = [&](const char* name, Real error, Real tolerance) {
    const bool pass = std::isfinite(error) && error <= tolerance;
    failures += !pass;
    Print() << (pass ? "  PASS  " : "  FAIL  ") << name << ": error = " << error
            << " (tol = " << tolerance << ")\n";
  };
  const Real tol =
      std::max(Real(3.e-11), Real(200) * std::numeric_limits<Real>::epsilon());
  initialize();
  const auto initial = diagnostics();
  const Long expectedCount = Long(config.nCell) * config.nCell * config.nCell *
                             config.particlesPerCell;
  report("deterministic particle initialization count",
         std::abs(Real(initial.count - expectedCount)), 0);
  report("particle weight normalization", std::abs(initial.charge - charge),
         tol);

  const Real xyz[7][3] = { { 0, 0, 0 },
                           { .5, .5, .5 },
                           { .5 / config.nCell, .5 / config.nCell,
                             .5 / config.nCell },
                           { .98, .02, .51 },
                           { .123, .456, .789 },
                           { 1 - 1.e-9, .5, .25 },
                           { .5, .25, 1 - 1.e-9 } };
  std::vector<ParticleType> samples(7);
  for (int ip = 0; ip < 7; ++ip) {
    auto& p = samples[ip];
    p.id() = ip + 1;
    p.cpu() = 0;
    p.idata(0) = species;
    for (int d = 0; d < 3; ++d) {
      p.pos(d) = xyz[ip][d];
      p.rdata(d) = .07 * (ip + 1) * (d == 1 ? -1 : 1) + .03 * d;
    }
    p.rdata(3) = charge / 7;
  }
  auto compare_particles = [&](const std::vector<ParticleType>& expected) {
    const auto actual = snapshots(true);
    if (actual.size() != expected.size())
      return Real(1.e30);
    Real error = 0;
    for (std::size_t ip = 0; ip < actual.size(); ++ip) {
      if (actual[ip].id() != expected[ip].id() ||
          actual[ip].idata(0) != expected[ip].idata(0))
        return Real(1.e30);
      for (int d = 0; d < 3; ++d) {
        Real delta = std::abs(actual[ip].pos(d) - expected[ip].pos(d));
        delta = std::min(delta, 1 - delta);
        error = std::max(error, delta);
        error = std::max(error,
                         std::abs(actual[ip].rdata(d) - expected[ip].rdata(d)));
      }
      error = std::max(error,
                       std::abs(actual[ip].rdata(3) - expected[ip].rdata(3)));
    }
    return error;
  };

  fields.nodeE.setVal(0);
  fields.nodeB.setVal(0);
  set_particles(samples);
  auto expected = samples;
  for (auto& p : expected)
    for (int d = 0; d < 3; ++d) {
      p.pos(d) += 2.17 * p.rdata(d);
      p.pos(d) -= std::floor(p.pos(d));
    }
  gather(fields.nodeE, fields.nodeB);
  push(2.17);
  report("free streaming across boxes/ranks/periodic seams",
         compare_particles(expected), tol);
  const auto migrated = diagnostics();
  report("migration preserves count and IDs",
         std::abs(Real(migrated.count) - 7) + std::abs(migrated.idSum - 28), 0);

  fields.nodeE.setVal(.2, 0, 1, fields.nodeE.nGrow());
  set_particles(samples);
  expected = samples;
  const Real uniformE[3] = { .2, 0, 0 }, zero[3] = { 0, 0, 0 };
  for (auto& p : expected)
    reference_boris(p, uniformE, zero, charge / mass, .071);
  gather(fields.nodeE, fields.nodeB);
  push(.071);
  report("uniform E acceleration", compare_particles(expected), tol);

  fields.nodeE.setVal(0);
  fields.nodeB.setVal(.4, 2, 1, fields.nodeB.nGrow());
  set_particles(samples);
  expected = samples;
  const Real uniformB[3] = { 0, 0, .4 };
  for (auto& p : expected)
    reference_boris(p, zero, uniformB, charge / mass, .071);
  const auto beforeGyro = diagnostics();
  gather(fields.nodeE, fields.nodeB);
  push(.071);
  report("Boris uniform B independent cross-product reference",
         compare_particles(expected), tol);
  report("uniform B kinetic energy",
         std::abs(diagnostics().kineticEnergy - beforeGyro.kineticEnergy), tol);
  Real gyroError[2] = { 0, 0 };
  for (int refinement = 0; refinement < 2; ++refinement) {
    set_particles(samples);
    const int steps = 10 * (refinement + 1);
    const Real gyroDt = .4 / steps;
    for (int step = 0; step < steps; ++step) {
      gather(fields.nodeE, fields.nodeB);
      push(gyroDt);
    }
    const auto actual = snapshots(true);
    const Real angle = charge / mass * .4 * .4;
    for (std::size_t ip = 0; ip < actual.size(); ++ip) {
      const Real vx = samples[ip].rdata(0) * std::cos(angle) +
                      samples[ip].rdata(1) * std::sin(angle);
      const Real vy = samples[ip].rdata(1) * std::cos(angle) -
                      samples[ip].rdata(0) * std::sin(angle);
      gyroError[refinement] = std::max(
          gyroError[refinement], std::max(std::abs(actual[ip].rdata(0) - vx),
                                          std::abs(actual[ip].rdata(1) - vy)));
    }
  }
  report("uniform B rotation converges under dt refinement",
         std::max(Real(0), gyroError[1] - .3 * gyroError[0]), tol);

  // Nonuniform periodic nodal field tests gather independently of a constant
  // field. Every reference interpolation evaluates the sampled grid values.
  const int n = config.nCell;
  for (MFIter mfi(fields.nodeE); mfi.isValid(); ++mfi) {
    const auto ea = fields.nodeE.array(mfi);
    ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      ea(i, j, k, 0) =
          std::sin(2 * pi * i / n) + .25 * std::cos(2 * pi * j / n);
      ea(i, j, k, 1) = .1 * std::sin(2 * pi * k / n);
      ea(i, j, k, 2) = 0;
    });
  }
  fields.nodeE.FillBoundary(grid.geom.periodicity());
  set_particles(samples);
  expected = samples;
  for (auto& p : expected) {
    const auto s = reference_shape(p, n);
    Real ep[3] = { 0, 0, 0 };
    for (int a = 0; a < 8; ++a) {
      const int i = s.index[a][0], j = s.index[a][1], k = s.index[a][2];
      ep[0] += s.weight[a] *
               (std::sin(2 * pi * i / n) + .25 * std::cos(2 * pi * j / n));
      ep[1] += s.weight[a] * .1 * std::sin(2 * pi * k / n);
    }
    reference_boris(p, ep, uniformB, charge / mass, .013);
  }
  gather(fields.nodeE, fields.nodeB);
  push(.013);
  report("CIC gather plus push independent serial reference",
         compare_particles(expected), tol);
  const auto separate = snapshots(true);
  set_particles(samples);
  push(.013, true, &fields.nodeE, &fields.nodeB);
  report("fused gather/push matches separate kernels",
         compare_particles(separate), tol);

  set_particles(samples);
  std::vector<Real> moments(std::size_t(n) * n * n * 10, 0);
  const Real invVol = 1 / cell_volume(grid);
  for (const auto& p : samples) {
    const auto s = reference_shape(p, n);
    const Real q = p.rdata(3) * invVol, u = p.rdata(0), v = p.rdata(1),
               w = p.rdata(2);
    const Real value[10] = { q,         q * u,     q * v,     q * w,
                             q * u * u, q * v * v, q * w * w, q * u * v,
                             q * u * w, q * v * w };
    for (int a = 0; a < 8; ++a)
      for (int c = 0; c < 10; ++c)
        moments[reference_offset(s.index[a][0], s.index[a][1], s.index[a][2], n,
                                 10, c)] += s.weight[a] * value[c];
  }
  fields.moments.setVal(0);
  deposit(fields.moments);
  fields.moments.SumBoundary(grid.geom.periodicity());
  fields.moments.FillBoundary(grid.geom.periodicity());
  report("atomic ten-moment CIC scatter serial reference",
         compare_reference(fields.moments, moments, n) / invVol, tol);

  // Full nodal response reference for seven particles, including collisions
  // on the same nodes and deposits at periodic corners.
  fields.clear_response();
  accumulate_response(fields, .071);
  fields.finish_response();
  std::vector<Real> current(std::size_t(n) * n * n * 3, 0),
      matrix(std::size_t(n) * n * n * 243, 0);
  const Real t = charge / mass * .071 * .4 / 2, den = 1 + t * t;
  const Real aMatrix[9] = {
    1 / den, t / den, 0, -t / den, 1 / den, 0, 0, 0, 1
  };
  for (const auto& p : samples) {
    const auto s = reference_shape(p, n);
    for (int a = 0; a < 8; ++a) {
      for (int d = 0; d < 3; ++d)
        for (int c = 0; c < 3; ++c)
          current[reference_offset(s.index[a][0], s.index[a][1], s.index[a][2],
                                   n, 3, d)] += s.weight[a] * invVol *
                                                p.rdata(3) *
                                                aMatrix[3 * d + c] * p.rdata(c);
      for (int b = 0; b < 8; ++b) {
        const int di = s.index[b][0] - s.index[a][0] + 1,
                  dj = s.index[b][1] - s.index[a][1] + 1,
                  dk = s.index[b][2] - s.index[a][2] + 1;
        for (int c = 0; c < 9; ++c)
          matrix[reference_offset(s.index[a][0], s.index[a][1], s.index[a][2],
                                  n, 243, dk * 81 + dj * 27 + di * 9 + c)] +=
              s.weight[a] * s.weight[b] * invVol * p.rdata(3) * charge / mass *
              .071 / 2 * aMatrix[c];
      }
    }
  }
  report("jHat scatter independent serial reference",
         compare_reference(fields.jHat, current, n) / invVol, tol);
  report("243-component mass matrix serial reference",
         compare_matrix_reference(fields.nodeMM, matrix, n) / invVol, tol);
  fields.clear_response();
  accumulate_response_buffered(fields, .071);
  fields.finish_response();
  report("cell-buffered jHat serial reference",
         compare_reference(fields.jHat, current, n) / invVol, tol);
  report("cell-buffered mass matrix serial reference",
         compare_matrix_reference(fields.nodeMM, matrix, n) / invVol, tol);

  fields.rho.setVal(0);
  fields.centerMM.setVal(0);
  deposit_center_charge(fields.rho, &fields.centerMM);
  fields.rho.SumBoundary(grid.geom.periodicity());
  fields.centerMM.SumBoundary(grid.geom.periodicity());
  fields.rho.FillBoundary(grid.geom.periodicity());
  fields.centerMM.FillBoundary(grid.geom.periodicity());
  std::vector<Real> rho(std::size_t(n) * n * n, 0),
      centerMatrix(std::size_t(n) * n * n * 27, 0);
  for (const auto& p : samples) {
    const auto s = reference_shape(p, n, true);
    for (int a = 0; a < 8; ++a) {
      rho[reference_offset(s.index[a][0], s.index[a][1], s.index[a][2], n,
                           1)] += p.rdata(3) * invVol * s.weight[a];
      for (int b = 0; b < 8; ++b) {
        Real dot = 0;
        for (int d = 0; d < 3; ++d)
          dot += s.gradient[a][d] * s.gradient[b][d];
        const int di = s.index[b][0] - s.index[a][0] + 1,
                  dj = s.index[b][1] - s.index[a][1] + 1,
                  dk = s.index[b][2] - s.index[a][2] + 1;
        centerMatrix[reference_offset(
            s.index[a][0], s.index[a][1], s.index[a][2], n, 27,
            di * 9 + dj * 3 + dk)] += std::abs(p.rdata(3)) * invVol * dot;
      }
    }
  }
  report("cell charge CIC serial reference",
         compare_reference(fields.rho, rho, n) / invVol, tol);
  report("divE position-response matrix serial reference",
         compare_reference(fields.centerMM, centerMatrix, n) / (invVol * n * n),
         tol);

  // A spatially constant potential cannot change particle positions.
  fields.phi.setVal(.005);
  correct_positions(fields.phi);
  report("constant divE potential leaves particles unchanged",
         compare_particles(samples), tol);
  initialize();
  return failures;
}

void DemoParticles::accumulate_response_buffered(DemoFields& fields, Real dt) {
  // Sorting and range construction are intentionally inside this operation:
  // the benchmark must account for the cost needed to form cell buffers.
  sort_by_cell();
  const auto lo = grid.geom.ProbLoArray();
  const auto invDx = grid.geom.InvCellSizeArray();
  const Real qdt2m = charge / mass * dt / 2, invVol = 1 / cell_volume(grid);
  for (ParConstIterType pti(*this, 0); pti.isValid(); ++pti) {
    const int np = pti.numParticles();
    if (!np)
      continue;
    const auto box = pti.validbox();
    const auto small = box.smallEnd();
    const auto length = box.length();
    const int ncells = box.numPts();
    Gpu::DeviceVector<int> first(ncells), last(ncells);
    int* begin = first.data();
    int* end = last.data();
    ParallelFor(ncells, [=] AMREX_GPU_DEVICE(int cell) {
      begin[cell] = np;
      end[cell] = 0;
    });
    const auto* particles = pti.GetArrayOfStructs().data();
    ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) {
      const auto s = demo_particle::shape(particles[ip], lo, invDx);
      const int cell =
          ((s.cell[2] - small[2]) * length[1] + s.cell[1] - small[1]) *
              length[0] +
          s.cell[0] - small[0];
      Gpu::Atomic::Min(&begin[cell], ip);
      Gpu::Atomic::Max(&end[cell], ip + 1);
    });
    const auto ba = fields.nodeB.const_array(pti);
    const auto current = fields.jHat.array(pti);
    const auto mm = fields.nodeMM.array(pti);
    ParallelFor(ncells, [=] AMREX_GPU_DEVICE(int cell) {
      if (begin[cell] >= end[cell])
        return;
      Real jbuffer[24] = { 0 };
      Real mbuffer[576] = { 0 };
      for (int ip = begin[cell]; ip < end[cell]; ++ip) {
        const auto& p = particles[ip];
        const auto s = demo_particle::shape(p, lo, invDx);
        Real b[3], alpha[9], jhat[3] = { 0, 0, 0 };
        demo_particle::gather(ba, s, b);
        demo_particle::alpha(b, qdt2m, alpha);
        for (int d = 0; d < 3; ++d)
          for (int c = 0; c < 3; ++c)
            jhat[d] += p.rdata(3) * alpha[3 * d + c] * p.rdata(c);
        for (int a = 0; a < 8; ++a) {
          for (int d = 0; d < 3; ++d)
            jbuffer[3 * a + d] += s.weight[a] * jhat[d];
          for (int c = 0; c < 8; ++c) {
            const Real factor =
                s.weight[a] * s.weight[c] * p.rdata(3) * qdt2m * invVol;
            for (int d = 0; d < 9; ++d)
              mbuffer[(a * 8 + c) * 9 + d] += factor * alpha[d];
          }
        }
      }
      const auto s = demo_particle::shape(particles[begin[cell]], lo, invDx);
      for (int a = 0; a < 8; ++a) {
        int i, j, k;
        demo_particle::index(s, a, i, j, k);
        for (int d = 0; d < 3; ++d)
          Gpu::Atomic::Add(&current(i, j, k, d), jbuffer[3 * a + d]);
        for (int c = 0; c < 8; ++c) {
          int ii, jj, kk;
          demo_particle::index(s, c, ii, jj, kk);
          const int offset =
              (kk - k + 1) * 81 + (jj - j + 1) * 27 + (ii - i + 1) * 9;
          for (int d = 0; d < 9; ++d)
            Gpu::Atomic::Add(&mm(i, j, k)[offset + d],
                             mbuffer[(a * 8 + c) * 9 + d]);
        }
      }
    });
    // Range storage must remain alive until its kernel has finished.
    Gpu::streamSynchronize();
  }
}

namespace {
struct SoAScratchParticle {
  Real* data;
  int size, index;
  AMREX_GPU_DEVICE Real& pos(int d) { return data[d * size + index]; }
  AMREX_GPU_DEVICE const Real& pos(int d) const {
    return data[d * size + index];
  }
  AMREX_GPU_DEVICE Real& rdata(int d) { return data[(3 + d) * size + index]; }
};
struct LayoutScratch {
  Gpu::DeviceVector<DemoParticles::ParticleType> aos;
  Gpu::DeviceVector<Real> soa;
};
} // namespace

int DemoParticles::benchmark_layouts(const MultiFab& e, const MultiFab& b) {
  std::map<std::pair<int, int>, LayoutScratch> scratch;
  const auto lo = grid.geom.ProbLoArray(), hi = grid.geom.ProbHiArray();
  const auto invDx = grid.geom.InvCellSizeArray();
  const GpuArray<Real, 3> length{ hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2] };
  const Real qom = charge / mass, dt = config.dt;
  double conversion = amrex::second();
  for (ParConstIterType pti(*this, 0); pti.isValid(); ++pti) {
    const int np = pti.numParticles();
    auto& tile = scratch[{ pti.index(), pti.LocalTileIndex() }];
    tile.aos.resize(np);
    tile.soa.resize(7 * np);
    const auto* source = pti.GetArrayOfStructs().data();
    auto* soa = tile.soa.data();
    Gpu::dtod_memcpy_async(tile.aos.data(), source, np * sizeof(ParticleType));
    ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) {
      for (int d = 0; d < 3; ++d)
        soa[d * np + ip] = source[ip].pos(d);
      for (int d = 0; d < 4; ++d)
        soa[(3 + d) * np + ip] = source[ip].rdata(d);
    });
  }
  Gpu::synchronize();
  conversion = amrex::second() - conversion;
  ParallelDescriptor::ReduceRealMax(conversion);
  double times[2];
  for (int layout = 0; layout < 2; ++layout) {
    // One warmup is included on BOTH paths and excluded from timing. Frozen
    // positions keep the stencil local so this isolates memory layout, not MPI.
    for (int iter = -1; iter < config.nIter; ++iter) {
      if (iter == 0) {
        Gpu::synchronize();
        ParallelDescriptor::Barrier();
        times[layout] = amrex::second();
      }
      for (ParConstIterType pti(*this, 0); pti.isValid(); ++pti) {
        const int np = pti.numParticles();
        auto& tile = scratch.at({ pti.index(), pti.LocalTileIndex() });
        auto* aos = tile.aos.data();
        auto* soa = tile.soa.data();
        const auto ea = e.const_array(pti), ba = b.const_array(pti);
        if (layout == 0) {
          ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) {
            auto& p = aos[ip];
            Real position[3], ep[3], bp[3];
            for (int d = 0; d < 3; ++d)
              position[d] = p.pos(d);
            const auto s = demo_particle::shape(p, lo, invDx);
            demo_particle::gather(ea, s, ep);
            demo_particle::gather(ba, s, bp);
            demo_particle::push(p, ep, bp, qom, dt, lo, length, false);
            for (int d = 0; d < 3; ++d)
              p.pos(d) = position[d];
          });
        } else {
          ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) {
            SoAScratchParticle p{ soa, np, ip };
            Real position[3], ep[3], bp[3];
            for (int d = 0; d < 3; ++d)
              position[d] = p.pos(d);
            const auto s = demo_particle::shape(p, lo, invDx);
            demo_particle::gather(ea, s, ep);
            demo_particle::gather(ba, s, bp);
            demo_particle::push(p, ep, bp, qom, dt, lo, length, false);
            for (int d = 0; d < 3; ++d)
              p.pos(d) = position[d];
          });
        }
      }
    }
    Gpu::synchronize();
    times[layout] = (amrex::second() - times[layout]) / config.nIter;
    ParallelDescriptor::ReduceRealMax(times[layout]);
  }
  Real error = 0;
  for (auto& entry : scratch) {
    const int np = entry.second.aos.size();
    std::vector<ParticleType> aos(np);
    std::vector<Real> soa(7 * np);
    Gpu::dtoh_memcpy_async(aos.data(), entry.second.aos.data(),
                           np * sizeof(ParticleType));
    Gpu::dtoh_memcpy_async(soa.data(), entry.second.soa.data(),
                           soa.size() * sizeof(Real));
    Gpu::streamSynchronize();
    for (int ip = 0; ip < np; ++ip)
      for (int d = 0; d < 3; ++d)
        error = std::max(error,
                         std::abs(aos[ip].rdata(d) - soa[(3 + d) * np + ip]));
  }
  ParallelDescriptor::ReduceRealMax(error);
  Print() << "Particle layout microbenchmark (frozen-position gather/push, no "
             "migration):\n"
          << "  AoS " << times[0] * 1.e3 << " ms, component SoA scratch "
          << times[1] * 1.e3 << " ms\n"
          << "  transpose/setup " << conversion * 1.e3
          << " ms (excluded), max velocity difference " << error << "\n";
  return (!std::isfinite(error) || error > 1.e-10) ? 1 : 0;
}

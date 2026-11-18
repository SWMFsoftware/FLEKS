#ifndef _TESTPARTICLES_H_
#define _TESTPARTICLES_H_

#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "Particles.h"

enum class TPOutputGroup : int {
  B = 0,
  E,
  Drifts,
  GradB,
  VExB,
  Mu,
  Pitch,
  Energy,
  Gamma
};

struct TPOutputSpec {
  int recordSize = 7;
  bool needGradB = false;
  std::vector<TPOutputGroup> groups;
  std::vector<std::string> colNames;

  // Parse a #TPOUTPUT variable string (presets P/PB/PBE/PBED/PBEG and/or
  // composable tokens b, e, drifts, gradb, vexb, mu, pitch, energy, gamma).
  // When empty, falls back to the compile-time ptRecordSize setting.
  static TPOutputSpec parse(std::string_view spec);
};

class TestParticles : public PTParticles {
public:
  struct TrajEntry {
    int cpu = 0;
    int id = 0;
    amrex::Real qp = 0.0;
    bool merged = false;
    std::vector<amrex::Real> data;
  };

private:
  static constexpr int iRegionBoundary_ = 1;
  static constexpr int iRegionUniform_ = 2;
  static constexpr int iRegionSideXp_ = 3;
  static constexpr int iRegionUser_ = 4;

public:
  TestParticles(Grid* gridIn, FluidInterface* const fluidIn,
                TimeCtr* const tcIn, const int speciesID,
                const amrex::Real charge, const amrex::Real mass, int id = 0,
                const ParticlesInfo& picInfo = ParticlesInfo{});

  ~TestParticles() = default;

  // Push and optionally record particles on level iLev without redistributing.
  // nodeBMF is the evolved field and nodeB0MF (optional) is the frozen
  // intrinsic field; nodeEMF and nodeE0MF (optional) behave analogously.
  void move_and_save_particles(int iLev, const amrex::MultiFab& nodeEMF,
                               const amrex::MultiFab* nodeE0MF,
                               const amrex::MultiFab& nodeBMF,
                               const amrex::MultiFab* nodeB0MF, amrex::Real dt,
                               amrex::Real dtNext, amrex::Real tNow,
                               bool doSave,
                               const amrex::MultiFab* nodeJacBMF = nullptr);

  void move_and_save_charged_particles(
      int iLev, const amrex::MultiFab& nodeEMF, const amrex::MultiFab* nodeE0MF,
      const amrex::MultiFab& nodeBMF, const amrex::MultiFab* nodeB0MF,
      amrex::Real dt, amrex::Real dtNext, amrex::Real tNow, bool doSave,
      const amrex::MultiFab* nodeJacBMF = nullptr);

  void move_and_save_neutrals(int iLev, amrex::Real dt, amrex::Real tNow,
                              bool doSave);

  void read_test_particle_list(const amrex::Vector<std::string>& listFiles);

  void add_test_particles_from_fluid(
      const amrex::Vector<Vel>& tpStates = amrex::Vector<Vel>());

  void add_test_particles_from_pic(PicParticles* pts);

  void reset_record_counter();

  bool write_particles(int cycle);

  void set_output_spec(const TPOutputSpec& spec) {
    outputSpec = spec;
    rebuild_col_scales();
  }

  const TPOutputSpec& get_output_spec() const { return outputSpec; }

  void set_IO_units(amrex::Real no2outLIn, amrex::Real no2outVIn,
                    amrex::Real no2outMIn, amrex::Real no2outBIn,
                    amrex::Real no2outEIn, amrex::Real no2siVIn,
                    amrex::Real no2siBIn, const std::string& ioUnitIn);

  void set_interval(amrex::IntVect in) { nIntervalCell = in; }

  void set_particle_region(std::string& sRegion,
                           amrex::Vector<std::shared_ptr<Shape> >& shapes) {
    if (!sRegion.empty()) {
      // TODO: std::string::starts_with (C++20)
      if (sRegion[0] == '+' || sRegion[0] == '-') {
        iPartRegion = iRegionUser_;
        tpRegions = Regions(shapes, sRegion);
      } else if (sRegion == "boundary") {
        iPartRegion = iRegionBoundary_;
      } else if (sRegion == "uniform") {
        iPartRegion = iRegionUniform_;
      } else if (sRegion == "sideX+") {
        iPartRegion = iRegionSideXp_;
      } else {
        amrex::Abort("Error:Unknown test particle region!");
      }
    }
  }

  template <typename T> void gather_accumulate_and_scatter(T& local, T& ahead) {
    int nProc = amrex::ParallelDescriptor::NProcs();

    // The following two vectors are only useful on the root processor. They are
    // allocated on all processors to avoid memory leak since they are passed to
    // functions like Gather()
    amrex::Vector<T> perProc, accumulated;
    perProc.resize(nProc, 0);
    accumulated.resize(nProc, 0);

    amrex::ParallelDescriptor::Gather(
        &local, 1, &perProc[0], 1,
        amrex::ParallelDescriptor::IOProcessorNumber());

    if (amrex::ParallelDescriptor::IOProcessor()) {
      for (int i = 1; i < accumulated.size(); ++i) {
        accumulated[i] = accumulated[i - 1] + perProc[i - 1];
      }
    }

    amrex::ParallelDescriptor::Scatter(
        &ahead, 1, &accumulated[0], 1,
        amrex::ParallelDescriptor::IOProcessorNumber());
  }

  void update_initial_particle_number() {
    nInitPart = TotalNumberOfParticles(false, false);
  }

  unsigned long int init_particle_number() const { return nInitPart; }
  void init_particle_number(unsigned long int in) { nInitPart = in; }

  void set_dt_save(amrex::Real dt, amrex::Real tNow = 0.0) {
    dtSave = dt;
    if (dtSave > 0.0) {
      if (tNow > 0.0) {
        long int nPass = std::floor(tNow / dtSave);
        tNextSave = (nPass + 1) * dtSave;
      } else {
        tNextSave = 0.0;
      }
    }
  }

  amrex::Real get_dt_save() const { return dtSave; }
  void set_t_next_save(amrex::Real t) { tNextSave = t; }
  amrex::Real get_t_next_save() const { return tNextSave; }
  void advance_next_save() {
    if (dtSave > 0.0) {
      tNextSave += dtSave;
    }
  }

  bool is_time_to_record(amrex::Real tNow, long int cycle, int dnSave) const {
    if (dtSave > 0.0) {
      return (tNow >= tNextSave - 1.0e-9 * dtSave);
    }
    return (dnSave > 0 && cycle % dnSave == 0);
  }

private:
  struct SlotStage {
    struct Meta {
      int cpu = 0;
      int id = 0;
      amrex::Real qp = 0.0;
    };
    std::vector<Meta> metas;
    std::vector<amrex::Real> values;
  };

  std::string outputDir;

  std::string printPrefix;
  std::string gridName;
  int gridID;

  amrex::IntVect nIntervalCell = { AMREX_D_DECL(1, 1, 1) };

  int iPartRegion = iRegionBoundary_;

  amrex::Real no2outL = 1.0, no2outV = 1.0, no2outM = 1.0, no2outB = 1.0,
              no2outE = 1.0;
  amrex::Real no2siV = 1.0, no2siB = 1.0;
  amrex::Real energyScale = 1.0, muScale = 1.0;
  std::string ioUnit = "planet";

  TPOutputSpec outputSpec = TPOutputSpec::parse("");
  std::vector<amrex::Real> colScales;

  unsigned long int nInitPart;

  std::vector<PID> vIDs;

  Regions tpRegions;

  amrex::Real dtSave = -1.0;
  amrex::Real tNextSave = 0.0;

  // Rank-local trajectory buffer keyed by packed (cpu, id). Decoupled from
  // AMReX Particle structs so particle containers remain compact and segments
  // of exited particles are preserved until write_particles() flushes them.
  std::unordered_map<std::uint64_t, TrajEntry> trajMap_;

  void rebuild_col_scales();
  void commit_stages(const std::vector<SlotStage>& stages);
  void coalesce_trajectories(std::vector<TrajEntry>& sortedEntries);
  void pack_particles_and_records(const std::vector<TrajEntry>& entries,
                                  char* recordBuff, char* listBuff,
                                  unsigned long long int shift) const;

  void interpolate_record_trajectory(
      amrex::Real xp, amrex::Real yp, amrex::Real zp, amrex::Real up,
      amrex::Real vp, amrex::Real wp, amrex::Real unp1, amrex::Real vnp1,
      amrex::Real wnp1, amrex::Real dtStep, amrex::Real dtElse,
      amrex::Real tNowSI, amrex::Real& tRec, amrex::Real& xRec,
      amrex::Real& yRec, amrex::Real& zRec, amrex::Real& uRec,
      amrex::Real& vRec, amrex::Real& wRec) const;
};

#endif

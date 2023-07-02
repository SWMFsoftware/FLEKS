#ifndef _TESTPARTICLES_H_
#define _TESTPARTICLES_H_

#include "BitArray.h"
#include "Particles.h"

class TestParticles : public Particles<nPTPartReal, nPTPartInt> {
public:
  static const int iTPt_ = 0;
  static const int iTPx_ = 1;
  static const int iTPy_ = 2;
  static const int iTPz_ = 3;
  static const int iTPu_ = 4;
  static const int iTPv_ = 5;
  static const int iTPw_ = 6;
  static const int iTPBx_ = 7;
  static const int iTPBy_ = 8;
  static const int iTPBz_ = 9;
  static const int iTPEx_ = 10;
  static const int iTPEy_ = 11;
  static const int iTPEz_ = 12;

private:
  static const int iRegionBoundary_ = 1;
  static const int iRegionUniform_ = 2;

public:
  TestParticles(Grid* gridIn, FluidInterface* const fluidIn,
                TimeCtr* const tcIn, const int speciesID,
                const amrex::Real charge, const amrex::Real mass, int id = 0);

  ~TestParticles() = default;

  inline int record_var_index(int iPart, int iVar = 0) {
    return nPicPartReal + ptRecordSize * iPart + iVar;
  }

  void move_and_save_particles(const amrex::MultiFab& nodeEMF,
                               const amrex::MultiFab& nodeBMF, amrex::Real dt,
                               amrex::Real dtNext, amrex::Real tNow,
                               bool doSave);

  void move_and_save_charged_particles(const amrex::MultiFab& nodeEMF,
                                       const amrex::MultiFab& nodeBMF,
                                       amrex::Real dt, amrex::Real dtNext,
                                       amrex::Real tNow, bool doSave);

  void move_and_save_neutrals(amrex::Real dt, amrex::Real tNow, bool doSave);

  void read_test_particle_list(const amrex::Vector<std::string>& listFiles);

  void add_test_particles_from_fluid(
      amrex::Vector<Vel> tpStates = amrex::Vector<Vel>());

  void add_test_particles_from_pic(Particles<>* pts);

  void reset_record_counter();

  bool write_particles(int cycle);

  unsigned long long int loop_particles(
      std::string action = "count_record_size", char* buff = nullptr,
      unsigned long long int sizeLimit = 0, unsigned long long int shift = 0);

  void print_record_buffer(char* buffer, unsigned long long int nBuffer);

  void set_IO_units(amrex::Real no2outLIn, amrex::Real no2outVIn,
                    amrex::Real no2outMIn, amrex::Real no2outBIn,
                    amrex::Real no2outEIn) {
    no2outL = no2outLIn;
    no2outV = no2outVIn;
    no2outM = no2outMIn;
    no2outB = no2outBIn;
    no2outE = no2outEIn;
  }

  void set_interval(amrex::IntVect in) { nIntervalCell = in; };

  void set_particle_region(std::string& sRegion) {
    if (!sRegion.empty()) {
      if (sRegion.find("boundary") != std::string::npos) {
        iPartRegion = iRegionBoundary_;
      } else if (sRegion.find("uniform") != std::string::npos) {
        iPartRegion = iRegionUniform_;
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
      for (int i = 1; i < accumulated.size(); i++) {
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

private:
  std::string outputDir;

  std::string printPrefix;
  std::string gridName;
  int gridID;

  amrex::IntVect nIntervalCell = { AMREX_D_DECL(1, 1, 1) };

  int iPartRegion = iRegionBoundary_;

  amrex::Real no2outL, no2outV, no2outM, no2outB, no2outE;

  unsigned long int nInitPart;

  std::vector<PID> vIDs;
};

#endif

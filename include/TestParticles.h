#ifndef _TESTPARTICLES_H_
#define _TESTPARTICLES_H_

#include "Particles.h"

class TestParticles : public Particles<nPTPartReal, nPTPartInt> {
public:
  static const int iRecordt_ = 0;
  static const int iRecordx_ = 1;
  static const int iRecordy_ = 2;
  static const int iRecordz_ = 3;
  static const int iRecordu_ = 4;
  static const int iRecordv_ = 5;
  static const int iRecordw_ = 6;

public:
  TestParticles(const amrex::BoxArray& regionBAIn, const amrex::Geometry& geom,
                const amrex::DistributionMapping& dm, const amrex::BoxArray& ba,
                FluidInterface* const fluidIn, TimeCtr* const tcIn,
                const int speciesID, const amrex::Real charge,
                const amrex::Real mass, int domainIDIn = 0);

  ~TestParticles() = default;

  inline int record_var_index(int iPart, int iVar = 0) {
    return nPicPartReal + ptRecordSize * iPart + iVar;
  }

  void move_and_save_particles(const amrex::MultiFab& nodeEMF,
                               const amrex::MultiFab& nodeBMF, amrex::Real dt,
                               amrex::Real dtNext, amrex::Real tNow);

  void add_test_particles(const amrex::iMultiFab& cellStatus);

  void reset_record_counter();

  bool write_particles(int cycle);

  unsigned long long int loop_particles(
      std::string action = "count_record_size", char* buff = nullptr,
      unsigned long long int sizeLimit = 0, unsigned long long int shift = 0);

  void print_record_buffer(char* buffer, unsigned long long int nBuffer);

  void set_IO_units(amrex::Real no2outLIn, amrex::Real no2outVIn,
                    amrex::Real no2outMIn) {
    no2outL = no2outLIn;
    no2outV = no2outVIn;
    no2outM = no2outMIn;
  }

  template <typename T> void gather_accumulate_and_scatter(T& local, T& ahead) {
    using namespace amrex;

    int nProc = ParallelDescriptor::NProcs();

    // The following two vectors are only useful on the root processor. They are
    // allocated on all processors to avoid memory leak since they are passed to
    // functions like Gather()
    Vector<T> perProc, accumulated;
    perProc.resize(nProc, 0);
    accumulated.resize(nProc, 0);

    ParallelDescriptor::Gather(&local, 1, &perProc[0], 1,
                               ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor()) {
      for (int i = 1; i < accumulated.size(); i++) {
        accumulated[i] = accumulated[i - 1] + perProc[i - 1];
      }
    }

    ParallelDescriptor::Scatter(&ahead, 1, &accumulated[0], 1,
                                ParallelDescriptor::IOProcessorNumber());
  }

  void update_initial_particle_number() {
    nInitPart = TotalNumberOfParticles(false, false);
  }

  unsigned long int init_particle_number() const { return nInitPart; }
  void init_particle_number(unsigned long int in) { nInitPart = in; }

private:
  int iStep;
  std::string outputDir;

  std::string printPrefix;
  std::string domainName;
  int domainID;

  amrex::Real no2outL, no2outV, no2outM;

  unsigned long int nInitPart;
};

#endif
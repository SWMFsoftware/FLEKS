#ifndef _TIMECTR_H_
#define _TIMECTR_H_

#include <math.h>
#include <string>

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include "PlotWriter.h"

class TimeCtr;

class EventCtr {
private:
  TimeCtr *tc;

  // In SI unit
  amrex::Real dtEvent;
  amrex::Real lastTime; 
  amrex::Real nextTime;
  bool useDt;
  int nCount;

  int dnEvent;
  int nLast; 
  int nNext;
  bool useDn;

public:
  EventCtr(TimeCtr *tcIn, const amrex::Real dtIn = -1, const int dnIn = -1);

  void init(const amrex::Real dtIn, const int dnIn = -1);

  bool is_time_to(bool doForce=false);
};

class PlotCtr : public EventCtr {
public:
  PlotWriter writer;

public:
  PlotCtr(TimeCtr *tcIn, const int idIn = 0, const amrex::Real dtIn = -1,
          const int dnIn = -1, const std::string plotStringIN = "",
          const double dxIn = 1, const std::string plotVarIn = "",
          const std::array<double, 3> plotMinIn_D = { { 1, 1, 1 } },
          const std::array<double, 3> plotMaxIn_D = { { -1, -1, -1 } },
          const int nSpeciesIn = 2)
      : EventCtr(tcIn, dtIn, dnIn),
        writer(idIn, plotStringIN, dxIn, plotVarIn, plotMinIn_D, plotMaxIn_D,
               nSpeciesIn) {}
};

class TimeCtr {
private:
  amrex::Real timeNowSI;
  amrex::Real dtSI;
  amrex::Real si2no;
  amrex::Real no2si;
  long int cycle;

  // public member variables.
public:
  amrex::Vector<PlotCtr> plots;

  // public methods
public:
  TimeCtr() : timeNowSI(0), dtSI(0), si2no(1), no2si(1), cycle(0) {}

  void set_si2no(const amrex::Real si2noIn) {
    si2no = si2noIn;
    no2si = 1. / si2no;
  }
  void set_no2si(const amrex::Real no2siIn) {
    no2si = no2siIn;
    si2no = 1. / no2si;
  }

  void set_cycle(long int cycleIn) { cycle = cycleIn; }
  long int get_cycle() const { return cycle; }

  void set_time(const amrex::Real timeIn) { timeNowSI = timeIn * no2si; }
  void set_time_si(const amrex::Real timeIn) { timeNowSI = timeIn; }
  amrex::Real get_time() const { return timeNowSI * si2no; }
  amrex::Real get_time_si() const { return timeNowSI; }

  void set_dt(const amrex::Real dtIn) { dtSI = dtIn * no2si; }
  void set_dt_si(const amrex::Real dtIn) { dtSI = dtIn; }
  amrex::Real get_dt() const { return dtSI * si2no; }
  amrex::Real get_dt_si() const { return dtSI; }

  void update() {
    timeNowSI += dtSI;
    cycle++;
  }

  void write_plots(bool doForceWrite=false);
};

#endif

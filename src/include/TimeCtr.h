#ifndef _TIMECTR_H_
#define _TIMECTR_H_

#include <math.h>
#include <string>

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

class TimeCtr;

class EventCtr {
private:
  TimeCtr *tc;

  amrex::Real dtEvent;
  amrex::Real nextTime;
  bool useDt;
  int nCount;

  int dnEvent;
  int nNext;
  bool useDn;

public:
  EventCtr(TimeCtr *tcIn) { tc = tcIn; }

  void init(const amrex::Real time, const amrex::Real dtIn,
            const int nCurrent = -1, int dnIn = -1) {
    useDt = dtIn > 0;
    useDn = dnIn > 0;

    if (useDt && useDn) {
      amrex::Abort("Error: The event can not be controlled by dt and dn at the "
                   "same time");
    }

    if (useDt) {
      dtEvent = dtIn;
      nCount = floor(time / dtEvent) + 1;
      nextTime = nCount * dtEvent;
    }

    if (useDn) {
      dnEvent = dnIn;
      nNext = (floor((amrex::Real)nCurrent / dnEvent) + 1) * dnEvent;
    }
  }

  bool do_reach_time(const amrex::Real timeIn, const int nCycle) {
    bool doReach = false;
    if (useDt && timeIn >= nextTime) {
      nCount++;
      nextTime = nCount * dtEvent;
      doReach = true;
    }

    if (useDn && nCycle >= nNext) {
      nNext += dnEvent;
      doReach = true;
    }
    return doReach;
  }
};

class PlotCtr: public EventCtr{
public:
  PlotCtr(TimeCtr *tcIn):EventCtr(tcIn){}

private:
  std::string varStr; 
  
};


class TimeCtr {
private:
  amrex::Real timeNow;
  amrex::Real dt;
  amrex::Real si2no;
  amrex::Real no2si;
  long int cycle;

  // public member variables.
public:
  amrex::Vector<PlotCtr> plots;
  
  // public methods
public:
  TimeCtr() : timeNow(0), dt(0), si2no(1), no2si(1), cycle(0) {}

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

  void set_time(const amrex::Real timeIn) { timeNow = timeIn; }
  void set_time_si(const amrex::Real timeIn) { timeNow = timeIn * si2no; }
  amrex::Real get_time() const { return timeNow; }
  amrex::Real get_time_si() const { return timeNow * no2si; }

  void set_dt(const amrex::Real dtIn) { dt = dtIn; }
  void set_dt_si(const amrex::Real dtIn) { dt = dtIn * si2no; }
  amrex::Real get_dt() const { return dt; }
  amrex::Real get_dt_si() const { return dt * no2si; }

  void update() {
    timeNow += dt;
    cycle++;
  }
};

#endif

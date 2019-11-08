#include "TimeCtr.h"

EventCtr::EventCtr(TimeCtr *tcIn, const amrex::Real dtIn, const int dnIn) {
  tc = tcIn;
  init(dtIn, dnIn);
}

void EventCtr::init(const amrex::Real dtIn, const int dnIn) {
  useDt = dtIn > 0;
  useDn = dnIn > 0;

  if (useDt && useDn) {
    amrex::Abort("Error: The event can not be controlled by dt and dn at the "
                 "same time");
  }

  if (useDt) {
    dtEvent = dtIn;
    nCount = floor(tc->get_time_si() / dtEvent) + 1;
    nextTime = nCount * dtEvent;
  }

  if (useDn) {
    dnEvent = dnIn;
    nNext = (floor((amrex::Real)(tc->get_cycle()) / dnEvent) + 1) * dnEvent;
  }
}

bool EventCtr::is_time_to() {

  bool isTime = false;

  // Q: Why use 'while' instead of 'if'?
  // A: Make sure 'nextTime' would be larger than the curren time. 
  while (useDt && (tc->get_time_si() >= nextTime - dtEvent * 1e-6)) {
    nCount++;
    nextTime = nCount * dtEvent;
    isTime = true;
  }

  if (useDn && tc->get_cycle() >= nNext) {
    nNext += dnEvent;
    isTime = true;
  }
  return isTime;
}

#ifndef _TIMER_H_
#define _TIMER_H_

#include <cstring>
#include <string>

#include <AMReX_BLProfiler.H>

#include "Timing_c.h"

inline void timing_start(const char* name) {
  if (name) {
    size_t nameLen = std::strlen(name);
    timing_start_c(&nameLen, const_cast<char*>(name));
  }
}

inline void timing_stop(const char* name) {
  if (name) {
    size_t nameLen = std::strlen(name);
    timing_stop_c(&nameLen, const_cast<char*>(name));
  }
}

class Timer {
private:
  std::string nameStr;
  const char* nameC = nullptr;
  bool isTiming = true;

public:
  explicit Timer(const char* nameIn) : nameC(nameIn) { timing_start(nameC); }

  explicit Timer(const std::string& nameIn) : nameStr(nameIn) {
    nameC = nameStr.c_str();
    timing_start(nameC);
  }

  ~Timer() { stop(); }

  void stop() {
    if (isTiming) {
      timing_stop(nameC);
      isTiming = false;
    }
  }
};

#define timing_func(name)                                                      \
  Timer funcTimer(name);                                                       \
  BL_PROFILE(name);

#endif

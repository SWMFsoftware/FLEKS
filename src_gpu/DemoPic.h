#ifndef DEMO_PIC_H_
#define DEMO_PIC_H_

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "DemoConfig.h"
#include "DemoFields.h"
#include "DemoGrid.h"
#include "DemoParticles.h"

class DemoPic {
public:
  explicit DemoPic(const DemoConfig& configIn);
  int run();
  void initialize();
  void deposit_moments();
  void deposit_charge(bool matrix);
  void report(int step);
  int advance();
  void benchmark();
  void accumulate(DemoParticles& particles, const std::string& method);
  void measure(const std::string& label,
               const std::function<void()>& operation);

private:
  const DemoConfig& config;
  DemoGrid grid;
  DemoFields fields;
  std::vector<std::unique_ptr<DemoParticles> > species;
  std::map<std::string, double> timings;
  bool profileActive = false;
  amrex::Long expectedCount = 0;
  int failures = 0;
};

#endif

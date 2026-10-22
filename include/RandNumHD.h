#ifndef _RANDNUMHD_H_
#define _RANDNUMHD_H_

#include <cstdint>

#include <AMReX_GpuQualifiers.H>

// Park-Miller minimal standard generator with Schrage's factorization.
//
// This is a host-device copy of RandNum (share/Library/src/RandNum.h): the
// same algorithm with the same arithmetic, so the CPU and the GPU draw
// bitwise-identical sequences from the same seed. FLEKS seeds one stream per
// cell (Particles::set_random_seed), so a kernel that handles one cell per
// thread reproduces the CPU results exactly. Keep the arithmetic below in
// sync with RandNum; tests/gpu/randnum_test.cpp checks the two agree.
class RandNumHD {
private:
  // The constants of RandNum, renamed: RandNum.h defines its names as macros.
  static constexpr std::int32_t kIA = 16807;
  static constexpr std::int32_t kIM = 2147483647;
  static constexpr double kAM = 1.0 / kIM;
  static constexpr std::int32_t kIQ = 127773;
  static constexpr double kInvIQ = 1.0 / kIQ;
  static constexpr std::int32_t kIR = 2836;
  static constexpr std::int32_t kMASK = 123459876;

  std::int32_t idum = 1;

public:
  AMREX_GPU_HOST_DEVICE RandNumHD() = default;

  AMREX_GPU_HOST_DEVICE explicit RandNumHD(long seed) { set_seed(seed); }

  AMREX_GPU_HOST_DEVICE void set_seed(long in) {
    idum = in;

    if ((*this)() == 0) {
      idum = in % kMASK;
    } else {
      idum = in;
    }
  }

  // Uniform random number in (0, 1).
  AMREX_GPU_HOST_DEVICE double operator()() {
    std::int32_t k;
    double ans;

    idum ^= kMASK;
    k = (idum)*kInvIQ;
    idum = kIA * (idum - k * kIQ) - kIR * k;
    if (idum < 0)
      idum += kIM;
    ans = kAM * idum;
    idum ^= kMASK;
    return ans;
  }
};

#endif

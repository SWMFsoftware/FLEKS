// Check that RandNumHD (include/RandNumHD.h), the host-device random number
// generator of FLEKS, draws exactly the same numbers as RandNum
// (share/Library/src/RandNum.h), the generator FLEKS used before the GPU port.
// The particle initial conditions and sources depend on it bit for bit.
//
// Build and run from the FLEKS root (the AMReX headers only provide the
// AMREX_GPU_HOST_DEVICE qualifier):
//   g++ -std=c++17 -Iinclude -Ishare/Library/src -Iutil/AMREX/InstallDir/include
//       tests/gpu/randnum_test.cpp -o /tmp/randnum_test && /tmp/randnum_test

#include <cstdio>
#include <cstring>

#include "RandNum.h"
#include "RandNumHD.h"

int main() {
  // Seeds of the form used by Particles::calc_random_seed, including large
  // values that overflow 32 bits, plus small and special ones.
  const long seeds[] = { 0,          1,           2,          7,
                         123459876,  2147483647L, 2147483648L, 4294967296L,
                         -5,         987654321012L, 3 * 7 * 64 * 1234567L,
                         (5 + 3) * 7 * 787L * 1000003L };
  long nDraw = 0;
  for (long seed : seeds) {
    RandNum ref;
    RandNumHD hd;
    ref.set_seed(seed);
    hd.set_seed(seed);
    for (int i = 0; i < 100000; ++i) {
      const double a = ref();
      const double b = hd();
      if (std::memcmp(&a, &b, sizeof(double)) != 0) {
        std::printf("FAILED: seed %ld draw %d: %.17g != %.17g\n", seed, i, a,
                    b);
        return 1;
      }
      ++nDraw;
    }
  }
  std::printf("PASSED: %ld draws are bitwise identical\n", nDraw);
  return 0;
}

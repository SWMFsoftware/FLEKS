#ifndef _PARTICLEPUSHER_H_
#define _PARTICLEPUSHER_H_

#include <cmath>
#include <string>
#include <string_view>

#include <AMReX.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include "Utility.h"

/**
 * @brief Supported Lorentz velocity integrators for charged particles.
 */
enum class PusherType : int { Boris = 0, Vay = 1, HigueraCary = 2 };

/**
 * @brief Parse a user-facing pusher name into PusherType.
 */
inline PusherType parse_pusher_type(std::string_view name) {
  const std::string token = normalize_string_token(std::string(name));
  if (token == "boris")
    return PusherType::Boris;
  if (token == "vay")
    return PusherType::Vay;
  if (token == "higueracary" || token == "hc")
    return PusherType::HigueraCary;
  amrex::Abort("Error: unknown particle pusher '" + std::string(name) +
               "'. Expected 'boris', 'vay', or 'higueracary'.");
  return PusherType::Boris;
}

/**
 * @brief Canonical string name for a PusherType.
 */
inline const char* pusher_type_name(PusherType type) {
  switch (type) {
    case PusherType::Boris:
      return "boris";
    case PusherType::Vay:
      return "vay";
    case PusherType::HigueraCary:
      return "higueracary";
  }
  return "boris";
}

namespace detail {

// Closed-form intermediate Lorentz factor shared by Vay (2008, Eq. 12) and
// Higuera-Cary (2017, Eq. 20) in normalized c = 1 units.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real solve_relativistic_gamma(
    const amrex::Real (&u)[3], const amrex::Real (&tau)[3]) {
  const amrex::Real u2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
  const amrex::Real tau2 = tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2];
  const amrex::Real uStar = u[0] * tau[0] + u[1] * tau[1] + u[2] * tau[2];
  const amrex::Real sigma = 1.0 + u2 - tau2;
  const amrex::Real disc = sigma * sigma + 4.0 * (tau2 + uStar * uStar);
  return std::sqrt(0.5 * (sigma + std::sqrt(disc)));
}

} // namespace detail

/**
 * @brief Advance a charged particle's 3-velocity over one step in (E, B).
 *
 * Uses normalized c = 1 (uNormSI). For non-relativistic motion, Boris, Vay,
 * and Higuera-Cary are algebraically identical and share the fast Boris kernel.
 * Returns false if isRelativistic is true and |vIn| >= 1 (superluminal).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool push_lorentz_velocity(
    PusherType pusher, bool isRelativistic, const amrex::Real (&vIn)[3],
    const amrex::Real (&ep)[3], const amrex::Real (&bp)[3], amrex::Real qdto2m,
    amrex::Real (&vOut)[3]) {
  if (!isRelativistic) {
    const amrex::Real ut[3] = { vIn[0] + qdto2m * ep[0],
                                vIn[1] + qdto2m * ep[1],
                                vIn[2] + qdto2m * ep[2] };
    const amrex::Real omega[3] = { qdto2m * bp[0], qdto2m * bp[1],
                                   qdto2m * bp[2] };
    boris_push_nonrelativistic(vIn, ut, omega, vOut);
    return true;
  }

  const amrex::Real v2 = vIn[0] * vIn[0] + vIn[1] * vIn[1] + vIn[2] * vIn[2];
  if (!(v2 < amrex::Real(1.0))) {
    return false;
  }

  const amrex::Real gammaN =
      amrex::Real(1.0) / std::sqrt(amrex::Real(1.0) - v2);
  const amrex::Real uN[3] = { gammaN * vIn[0], gammaN * vIn[1],
                              gammaN * vIn[2] };
  const amrex::Real eps[3] = { qdto2m * ep[0], qdto2m * ep[1], qdto2m * ep[2] };
  const amrex::Real tau[3] = { qdto2m * bp[0], qdto2m * bp[1], qdto2m * bp[2] };

  amrex::Real uNp1[3] = { 0.0, 0.0, 0.0 };

  if (pusher == PusherType::Vay) {
    // Vay (2008), Phys. Plasmas 15, 056701, Eqs. (9)-(12).
    const amrex::Real uPrime[3] = {
      uN[0] + amrex::Real(2.0) * eps[0] + (vIn[1] * tau[2] - vIn[2] * tau[1]),
      uN[1] + amrex::Real(2.0) * eps[1] + (vIn[2] * tau[0] - vIn[0] * tau[2]),
      uN[2] + amrex::Real(2.0) * eps[2] + (vIn[0] * tau[1] - vIn[1] * tau[0])
    };
    const amrex::Real invGammaNp1 =
        amrex::Real(1.0) / detail::solve_relativistic_gamma(uPrime, tau);
    const amrex::Real t[3] = { tau[0] * invGammaNp1, tau[1] * invGammaNp1,
                               tau[2] * invGammaNp1 };
    const amrex::Real s = amrex::Real(1.0) / (amrex::Real(1.0) + t[0] * t[0] +
                                              t[1] * t[1] + t[2] * t[2]);
    const amrex::Real uPrimeDotT =
        uPrime[0] * t[0] + uPrime[1] * t[1] + uPrime[2] * t[2];
    uNp1[0] = s * (uPrime[0] + uPrimeDotT * t[0] +
                   (uPrime[1] * t[2] - uPrime[2] * t[1]));
    uNp1[1] = s * (uPrime[1] + uPrimeDotT * t[1] +
                   (uPrime[2] * t[0] - uPrime[0] * t[2]));
    uNp1[2] = s * (uPrime[2] + uPrimeDotT * t[2] +
                   (uPrime[0] * t[1] - uPrime[1] * t[0]));
  } else {
    const amrex::Real uMinus[3] = { uN[0] + eps[0], uN[1] + eps[1],
                                    uN[2] + eps[2] };
    amrex::Real invGamma = 0.0;
    if (pusher == PusherType::HigueraCary) {
      // Higuera & Cary (2017), Phys. Plasmas 24, 052104.
      invGamma =
          amrex::Real(1.0) / detail::solve_relativistic_gamma(uMinus, tau);
    } else {
      invGamma = amrex::Real(1.0) /
                 std::sqrt(amrex::Real(1.0) + uMinus[0] * uMinus[0] +
                           uMinus[1] * uMinus[1] + uMinus[2] * uMinus[2]);
    }
    const amrex::Real t[3] = { tau[0] * invGamma, tau[1] * invGamma,
                               tau[2] * invGamma };
    boris_push_nonrelativistic(uN, uMinus, t, uNp1);
  }

  const amrex::Real invGammaOut =
      amrex::Real(1.0) / std::sqrt(amrex::Real(1.0) + uNp1[0] * uNp1[0] +
                                   uNp1[1] * uNp1[1] + uNp1[2] * uNp1[2]);
  vOut[0] = uNp1[0] * invGammaOut;
  vOut[1] = uNp1[1] * invGammaOut;
  vOut[2] = uNp1[2] * invGammaOut;
  return true;
}

#endif

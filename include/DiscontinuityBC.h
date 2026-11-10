#ifndef FLEKS_DISCONTINUITYBC_H
#define FLEKS_DISCONTINUITYBC_H

#include <string>
#include <vector>

#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include "FluidInterface.h"
#include "ReadParam.h"

enum class DiscontinuityType { TD = 0, RD = 1, General = 2 };

struct DiscontinuitySpeciesState {
  amrex::Real nDens = 0.0; // code units
  amrex::Real ux = 0.0;    // code units
  amrex::Real uy = 0.0;
  amrex::Real uz = 0.0;
  amrex::Real vth = 0.0;  // 1-D thermal std sqrt(kT/m) in code units
  amrex::Real mass = 1.0; // m / m_p
  amrex::Real charge = 1.0;
};

// Runtime manager for a planar solar wind discontinuity (#DISCONTINUITY)
// propagating into the domain across 'inflow' boundaries.
class DiscontinuityBoundaryManager {
public:
  bool active = false;
  bool isCodeUnits = false;

  std::string typeStr = "td";
  DiscontinuityType type = DiscontinuityType::TD;
  int rdRotationSign = 1; // +1: shortest arc (+pi if 180 deg), -1: reverse arc

  // Raw parameters read from #DISCONTINUITY
  amrex::Real normal[3] = { 1.0, 0.0, 0.0 };
  amrex::Real width = 0.0;    // code length units
  amrex::Real tStartSI = 0.0; // SI seconds [s]
  amrex::Real tStart = 0.0;   // code time units
  amrex::Real vnInput = -1.0; // optional override [km/s]; <= 0 means auto
  amrex::Real vn = 0.0;       // normal propagation speed in code units

  std::vector<amrex::Real> rho2In;         // [1/cc]
  std::vector<amrex::Real> ux2In;          // [km/s]
  std::vector<amrex::Real> uy2In;          // [km/s]
  std::vector<amrex::Real> uz2In;          // [km/s]
  std::vector<amrex::Real> t2In;           // [K]
  amrex::Real b2In[3] = { 0.0, 0.0, 0.0 }; // [T]

  // Converted code-unit states for State 1 (pre-discontinuity) and State 2
  std::vector<DiscontinuitySpeciesState> state1;
  std::vector<DiscontinuitySpeciesState> state2;
  amrex::Real b1[3] = { 0.0, 0.0, 0.0 };
  amrex::Real b2[3] = { 0.0, 0.0, 0.0 };
  amrex::Real uBulk1[3] = { 0.0, 0.0, 0.0 };
  amrex::Real uBulk2[3] = { 0.0, 0.0, 0.0 };

  // First-contact domain projection dMin = min_{x in Domain}(n . x)
  amrex::Real dMin = 0.0;

  // Precomputed RD basis and rotation parameters in plane perpendicular to n
  amrex::Real bn = 0.0;
  amrex::Real e1[3] = { 0.0, 1.0, 0.0 };
  amrex::Real e2[3] = { 0.0, 0.0, 1.0 };
  amrex::Real bt1 = 0.0;
  amrex::Real bt2 = 0.0;
  amrex::Real theta1 = 0.0;
  amrex::Real dTheta = 0.0;
  amrex::Real walenSign = 1.0;
  amrex::Real invSqrtFourPiRho1 = 0.0;

  // Hybrid fluid electron parameters for TD pressure balance
  bool isHybridPIC = false;
  bool useHybridElectronP = false;
  amrex::Real te0Code = 0.0;
  amrex::Real ne0Code = 1.0;
  amrex::Real gammaE = 1.0;

  void read_param(ReadParam& param);

  void convert_and_validate(const FluidInterface* fi,
                            const amrex::Geometry& geom, int nDim,
                            bool isFake2D, bool useHybridPIC,
                            amrex::Real electronTemperatureCode,
                            amrex::Real electronDensity0Code,
                            amrex::Real electronGammaIn);

  // Signed phase distance s(x, t) ahead of the discontinuity center plane.
  amrex::Real phase_distance(const amrex::Real* pos, amrex::Real t) const;

  // Smooth transition weight alpha(x, t) in [0, 1] (0 = State 1, 1 = State 2).
  amrex::Real alpha_at(const amrex::Real* pos, amrex::Real t) const;

  // Evaluate magnetic field B(x, t) in code units.
  void eval_B(const amrex::Real* pos, amrex::Real t, amrex::Real bOut[3]) const;

  // Evaluate bulk plasma velocity u_bulk(x, t) in code units.
  void eval_u_bulk(const amrex::Real* pos, amrex::Real t,
                   amrex::Real uOut[3]) const;

  // Evaluate convective electric field E(x, t) = -u_bulk(x, t) x B(x, t).
  void eval_E(const amrex::Real* pos, amrex::Real t, amrex::Real eOut[3]) const;

  // Evaluate per-species plasma moments (nDens, u[3], vth) in code units.
  void eval_species(int iS, const amrex::Real* pos, amrex::Real t,
                    amrex::Real& nOut, amrex::Real uOut[3],
                    amrex::Real& vthOut) const;
};

#endif // FLEKS_DISCONTINUITYBC_H

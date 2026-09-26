#include <Constants.h>

#include <AMReX_BC_TYPES.H>
#include <AMReX_PhysBCFunct.H>

#include "BC.h"

using namespace amrex;

// Input spellings accepted by ParticleBC::parse().
const std::vector<bc_detail::Entry>& ParticleBC::table() {
  using bc_detail::Legacy;
  static const std::vector<bc_detail::Entry> tbl = {
    { "periodic", ParticleBC::periodic, Legacy::none },
    { "coupled", ParticleBC::coupled, Legacy::none },
    { "outflow", ParticleBC::outflow, Legacy::none },
    { "open", ParticleBC::outflow, Legacy::silent },
    { "vacuum", ParticleBC::vacuum, Legacy::none },
    { "vacume", ParticleBC::vacuum, Legacy::silent },
    { "reflect", ParticleBC::reflect, Legacy::none },
    { "absorb", ParticleBC::absorb, Legacy::none },
    { "inflow", ParticleBC::inflow, Legacy::none },
    // Field-domain spellings the old shared enum silently accepted.
    { "conducting", ParticleBC::reflect, Legacy::deprecated },
    { "fixed", ParticleBC::inflow, Legacy::deprecated },
    { "wave", ParticleBC::outflow, Legacy::deprecated },
  };
  return tbl;
}

// Input spellings accepted by FieldBC::parse().
const std::vector<bc_detail::Entry>& FieldBC::table() {
  using bc_detail::Legacy;
  static const std::vector<bc_detail::Entry> tbl = {
    { "periodic", FieldBC::periodic, Legacy::none },
    { "coupled", FieldBC::coupled, Legacy::none },
    { "outflow", FieldBC::outflow, Legacy::none },
    { "open", FieldBC::outflow, Legacy::silent },
    { "vacuum", FieldBC::vacuum, Legacy::none },
    { "vacume", FieldBC::vacuum, Legacy::silent },
    { "conducting", FieldBC::conducting, Legacy::none },
    { "absorb", FieldBC::absorb, Legacy::none },
    { "inflow", FieldBC::inflow, Legacy::none },
    { "fixed", FieldBC::fixed, Legacy::none },
    { "wave", FieldBC::wave, Legacy::none },
    // Particle-domain spellings the old shared enum silently accepted.
    { "reflect", FieldBC::conducting, Legacy::deprecated },
  };
  return tbl;
}

// Input spellings accepted by BodyFieldBC::parse().
const std::vector<bc_detail::Entry>& BodyFieldBC::table() {
  using bc_detail::Legacy;
  static const std::vector<bc_detail::Entry> tbl = {
    { "linetied", BodyFieldBC::linetied, Legacy::none },
    { "conducting", BodyFieldBC::conducting, Legacy::none },
    { "insulating", BodyFieldBC::insulating, Legacy::none },
  };
  return tbl;
}

amrex::Vector<amrex::BCRec> FieldBC::create_bcrec(
    const BoxBC<FieldBC::Type>& physBC, Quantity qty, int nComp,
    const amrex::Geometry& geom) {
  amrex::Vector<amrex::BCRec> bcr(nComp);

  for (int c = 0; c < nComp; ++c) {
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      for (int side = 0; side < 2; ++side) {
        int mathType = amrex::BCType::bogus;

        if (geom.isPeriodic(d)) {
          mathType = amrex::BCType::int_dir;
        } else {
          const auto pType = static_cast<FieldBC::Type>(physBC.face(d, side));
          switch (pType) {
            case FieldBC::periodic:
              mathType = amrex::BCType::int_dir;
              break;
            case FieldBC::outflow:
            case FieldBC::vacuum:
            case FieldBC::inflow:
            case FieldBC::absorb:
            case FieldBC::wave:
              mathType = amrex::BCType::foextrap;
              break;
            case FieldBC::conducting:
              if (qty == Quantity::Magnetic) {
                mathType = (c == d) ? amrex::BCType::reflect_odd
                                    : amrex::BCType::reflect_even;
              } else if (qty == Quantity::Electric) {
                mathType = (c == d) ? amrex::BCType::reflect_even
                                    : amrex::BCType::reflect_odd;
              } else {
                mathType = amrex::BCType::reflect_even;
              }
              break;
            case FieldBC::coupled:
            case FieldBC::fixed:
              mathType = amrex::BCType::ext_dir;
              break;
            default:
              mathType = amrex::BCType::foextrap;
              break;
          }
        }

        if (side == 0) {
          bcr[c].setLo(d, mathType);
        } else {
          bcr[c].setHi(d, mathType);
        }
      }
    }
  }

  return bcr;
}

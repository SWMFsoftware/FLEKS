#ifndef _BC_H_
#define _BC_H_

#include <map>
#include <string>
#include <vector>

#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_iMultiFab.H>

void apply_float_boundary(const amrex::iMultiFab &status, amrex::MultiFab &mf,
                          const amrex::Geometry &gm, const int iStart,
                          const int nComp, const int nshift = 0);

//==========================================================================
/// Storage shared by the two boundary-condition domains below.
///
/// FLEKS has two physically distinct boundary-condition domains:
///
///  * `ParticleBC` -- kinetic boundary conditions, defined **per species**.
///  * `FieldBC`    -- electromagnetic boundary conditions, defined once for
///                    the whole box.  A single face type is read by both the
///                    B operator and the E operator; the difference between
///                    B and E is carried by the operators themselves (their
///                    `isB` / `iField` argument), not by the input.
///
/// See Doc/Boundary_Conditions.md for the full design.
//==========================================================================
template <typename EnumT> struct BoxBC {
  amrex::IntVect lo, hi;

  BoxBC()
      : lo(static_cast<int>(EnumT::coupled)),
        hi(static_cast<int>(EnumT::coupled)) {}

  /// Boundary type on face `side` (0 = lo, 1 = hi) of dimension `d`.
  int face(const int d, const int side) const {
    return side == 0 ? lo[d] : hi[d];
  }

  void set(const int d, const int side, const EnumT type) {
    if (side == 0)
      lo[d] = static_cast<int>(type);
    else
      hi[d] = static_cast<int>(type);
  }
};

//==========================================================================
/// Implementation details shared by the two parse tables.
//==========================================================================
namespace bc_detail {

/// How an accepted-but-non-canonical input string is treated.
enum class Legacy {
  none,       ///< canonical spelling
  silent,     ///< documented alias; accepted without complaint
  deprecated  ///< accepted, but reported through parse()
};

/// One accepted input string.  `type` is an int so a single table layout can
/// serve both enum domains.
struct Entry {
  const char *name;
  int type;
  Legacy legacy;
};

/// Canonical (non-legacy) spelling of `type`.
inline const char *canonical_name(const std::vector<Entry> &table,
                                  const int type) {
  for (const Entry &e : table) {
    if (e.type == type && e.legacy == Legacy::none)
      return e.name;
  }
  return "?";
}

/// Look up `str` in `table`.
///
/// Returns false when `str` is unknown, or when it is a deprecated spelling
/// and `strict` is true; `valid` then holds the accepted spellings.  On a
/// legacy hit, `mapped` receives the canonical type and `warn` a
/// human-readable note (empty for `Legacy::silent`).
inline bool lookup(const std::string &str, const std::vector<Entry> &table,
                   const bool strict, int &mapped, std::string &warn,
                   std::string &valid) {
  valid.clear();
  warn.clear();

  for (const Entry &e : table) {
    if (e.legacy != Legacy::none)
      continue;
    if (!valid.empty())
      valid += ", ";
    valid += e.name;
  }

  for (const Entry &e : table) {
    if (str != e.name)
      continue;

    if (e.legacy == Legacy::none) {
      mapped = e.type;
      return true;
    }
    // Documented aliases stay accepted even in strict mode; only deprecated
    // spellings are rejected.
    if (e.legacy == Legacy::deprecated && strict)
      return false;

    mapped = e.type;
    if (e.legacy == Legacy::deprecated) {
      warn = "the boundary type '" + str + "' is deprecated; use '" +
             canonical_name(table, e.type) + "' instead.";
    }
    return true;
  }

  return false;
}

} // namespace bc_detail

//==========================================================================
/// Kinetic (particle) boundary conditions, one instance per species.
//==========================================================================
namespace ParticleBC {

enum Type {
  unset = -1, ///< internal sentinel, never set from input
  periodic = 0,
  coupled = 1, ///< default: state supplied by the MHD/fluid interface
  outflow = 2,
  vacuum = 3,
  reflect = 4,
  absorb = 5,
  inflow = 6,
  thermal = 7 ///< reserved: delete + re-emit at the wall temperature
};

/// Accepted input spellings (defined in src/BC.cpp).
const std::vector<bc_detail::Entry> &table();

/// True for a canonical type (i.e. one the user may end up with).
inline bool is_valid(const Type t) {
  for (const bc_detail::Entry &e : table()) {
    if (e.type == t && e.legacy == bc_detail::Legacy::none)
      return true;
  }
  return false;
}

/// Canonical input spelling of `t`, for diagnostics.
inline const char *to_string(const Type t) {
  return bc_detail::canonical_name(table(), t);
}

/// Parse an input string.  Aborts on an unknown string.  A non-empty `warn`
/// means a deprecated spelling was mapped; the caller decides whether to
/// print it (and should de-duplicate).
inline Type parse(const std::string &str, const bool strict,
                  std::string &warn) {
  int mapped = unset;
  std::string valid;
  if (!bc_detail::lookup(str, table(), strict, mapped, warn, valid)) {
    amrex::Abort("Error: unrecognized particle boundary type '" + str +
                 "'. Accepted values: " + valid);
  }
  return static_cast<Type>(mapped);
}

} // namespace ParticleBC

//==========================================================================
/// Electromagnetic field boundary conditions, one instance for the box.
//==========================================================================
namespace FieldBC {

enum Type {
  unset = -1, ///< internal sentinel, never set from input
  periodic = 0,
  coupled = 1, ///< default: state supplied by the MHD/fluid interface
  outflow = 2,
  vacuum = 3,
  conducting = 4,
  symmetry = 5,
  absorb = 6,
  inflow = 7,
  fixed = 8, ///< deprecated: currently equivalent to `inflow`
  wave = 9
};

/// Accepted input spellings (defined in src/BC.cpp).
const std::vector<bc_detail::Entry> &table();

/// True for a canonical type (i.e. one the user may end up with).
inline bool is_valid(const Type t) {
  for (const bc_detail::Entry &e : table()) {
    if (e.type == t && e.legacy == bc_detail::Legacy::none)
      return true;
  }
  return false;
}

/// Canonical input spelling of `t`, for diagnostics.
inline const char *to_string(const Type t) {
  return bc_detail::canonical_name(table(), t);
}

/// Parse an input string.  Aborts on an unknown string.  A non-empty `warn`
/// means a deprecated spelling was mapped; the caller decides whether to
/// print it (and should de-duplicate).
inline Type parse(const std::string &str, const bool strict,
                  std::string &warn) {
  int mapped = unset;
  std::string valid;
  if (!bc_detail::lookup(str, table(), strict, mapped, warn, valid)) {
    amrex::Abort("Error: unrecognized field boundary type '" + str +
                 "'. Accepted values: " + valid);
  }
  return static_cast<Type>(mapped);
}

} // namespace FieldBC

#endif

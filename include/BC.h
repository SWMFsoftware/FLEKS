#ifndef _BC_H_
#define _BC_H_

#include <string>
#include <vector>

#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Geometry.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Vector.H>
#include <AMReX_iMultiFab.H>

//==========================================================
template <typename EnumT> struct BoxBC {
  amrex::IntVect lo, hi;

  BoxBC()
      : lo(static_cast<int>(EnumT::coupled)),
        hi(static_cast<int>(EnumT::coupled)) {}

  // side: 0 = lo, 1 = hi
  int face(const int d, const int side) const {
    return side == 0 ? lo[d] : hi[d];
  }

  void set(const int d, const int side, const EnumT type) {
    if (side == 0)
      lo[d] = static_cast<int>(type);
    else
      hi[d] = static_cast<int>(type);
  }

  int loFace(const int d) const { return lo[d]; }
  int hiFace(const int d) const { return hi[d]; }
  void setLo(const int d, const EnumT type) { lo[d] = static_cast<int>(type); }
  void setHi(const int d, const EnumT type) { hi[d] = static_cast<int>(type); }

  bool has(const EnumT type) const {
    for (int d = 0; d < amrex::SpaceDim; ++d) {
      if (lo[d] == static_cast<int>(type) || hi[d] == static_cast<int>(type))
        return true;
    }
    return false;
  }
};

//==========================================================
namespace bc_detail {

enum class Legacy { none, silent, deprecated };

struct Entry {
  const char *name;
  int type;
  Legacy legacy;
};

inline const char *canonical_name(const std::vector<Entry> &table,
                                  const int type) {
  for (const Entry &e : table) {
    if (e.type == type && e.legacy == Legacy::none)
      return e.name;
  }
  return "?";
}

inline bool lookup(const std::string &str, const std::vector<Entry> &table,
                   int &mapped, std::string &valid) {
  valid.clear();

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

    if (e.legacy == Legacy::none || e.legacy == Legacy::silent) {
      mapped = e.type;
      return true;
    }
    if (e.legacy == Legacy::deprecated)
      return false;
  }

  return false;
}

} // namespace bc_detail

//==========================================================
namespace ParticleBC {

enum Type {
  unset = -1,
  periodic = 0,
  coupled = 1,
  outflow = 2,
  vacuum = 3,
  reflect = 4,
  absorb = 5,
  inflow = 6
};

const std::vector<bc_detail::Entry> &table();

inline bool is_valid(const Type t) {
  for (const bc_detail::Entry &e : table()) {
    if (e.type == t && e.legacy == bc_detail::Legacy::none)
      return true;
  }
  return false;
}

inline const char *to_string(const Type t) {
  return bc_detail::canonical_name(table(), t);
}

inline Type parse(const std::string &str) {
  int mapped = unset;
  std::string valid;
  if (!bc_detail::lookup(str, table(), mapped, valid)) {
    amrex::Abort("Error: unrecognized particle boundary type '" + str +
                 "'. Accepted values: " + valid);
  }
  return static_cast<Type>(mapped);
}

} // namespace ParticleBC

//==========================================================
namespace FieldBC {

enum Type {
  unset = -1,
  periodic = 0,
  coupled = 1,
  outflow = 2,
  vacuum = 3,
  conducting = 4,
  absorb = 5,
  inflow = 6,
  fixed = 7,
  wave = 8
};

enum class Quantity { Magnetic, Electric, Scalar };

// Translate physical FieldBC boundary conditions on domain faces into
// mathematical AMReX BCRec records (amrex::BCType) for each component of a
// MultiFab.
amrex::Vector<amrex::BCRec> create_bcrec(const BoxBC<FieldBC::Type> &physBC,
                                         Quantity qty, int nComp,
                                         const amrex::Geometry &geom);

const std::vector<bc_detail::Entry> &table();

inline bool is_valid(const Type t) {
  for (const bc_detail::Entry &e : table()) {
    if (e.type == t && e.legacy == bc_detail::Legacy::none)
      return true;
  }
  return false;
}

inline const char *to_string(const Type t) {
  return bc_detail::canonical_name(table(), t);
}

inline Type parse(const std::string &str) {
  int mapped = unset;
  std::string valid;
  if (!bc_detail::lookup(str, table(), mapped, valid)) {
    amrex::Abort("Error: unrecognized field boundary type '" + str +
                 "'. Accepted values: " + valid);
  }
  return static_cast<Type>(mapped);
}

} // namespace FieldBC

//==========================================================
namespace BodyFieldBC {

// Electromagnetic condition on the surface of the inner body (see the #BODY
// and #BODYBOUNDARY commands). Unlike FieldBC, these are applied on the
// nodes/cells of a body that sits inside the domain, and the surface normal
// is the radial direction from the body center, not a domain-face axis.
enum Type {
  unset = -1,
  // Perfectly absorbing / line-tied: E vanishes on the body nodes (the body
  // nodes are not unknowns of the implicit E solve) and B is frozen at its
  // initial value inside, because no plasma and no E are left there.
  linetied = 0,
  // Perfect conductor: the tangential electric field and the radial magnetic
  // field vanish, E_t = 0 and B_r = 0 (the radial E remains an unknown and
  // the tangential B carries the surface current).
  conducting = 1,
  // Perfect insulator: no conduction and no macroscopic surface shielding
  // current, so the magnetic field passes through undistorted
  // (B_inside = B_outside), n x (E_out - E_in) = 0 and the jump of the normal
  // D is the accumulated surface charge. Discretely: no field constraint.
  insulating = 2
};

const std::vector<bc_detail::Entry> &table();

inline bool is_valid(const Type t) {
  for (const bc_detail::Entry &e : table()) {
    if (e.type == t && e.legacy == bc_detail::Legacy::none)
      return true;
  }
  return false;
}

inline const char *to_string(const Type t) {
  return bc_detail::canonical_name(table(), t);
}

inline Type parse(const std::string &str) {
  int mapped = unset;
  std::string valid;
  if (!bc_detail::lookup(str, table(), mapped, valid)) {
    amrex::Abort("Error: unrecognized body field boundary type '" + str +
                 "'. Accepted values: " + valid);
  }
  return static_cast<Type>(mapped);
}

} // namespace BodyFieldBC

// Boundary geometry metadata helper for domain walls
struct BoundaryBounds {
  amrex::Dim3 domLo;
  amrex::Dim3 domHi;
  bool isNode[3] = { false, false, false };
  int loBnd[3] = { 0, 0, 0 };
  int hiBnd[3] = { 0, 0, 0 };
  int bcLo[3] = { 0, 0, 0 };
  int bcHi[3] = { 0, 0, 0 };

  BoundaryBounds() = default;
  BoundaryBounds(const amrex::Geometry &geom, amrex::IndexType ixType,
                 const BoxBC<FieldBC::Type> *bc = nullptr) {
    domLo = geom.Domain().smallEnd().dim3();
    domHi = geom.Domain().bigEnd().dim3();
    const int *dLo = geom.Domain().smallEnd().getVect();
    const int *dHi = geom.Domain().bigEnd().getVect();
    for (int d = 0; d < 3; ++d) {
      if (d < amrex::SpaceDim) {
        isNode[d] = (ixType[d] == amrex::IndexType::NODE);
        loBnd[d] = dLo[d];
        hiBnd[d] = isNode[d] ? (dHi[d] + 1) : dHi[d];
        if (bc) {
          bcLo[d] = bc->face(d, 0);
          bcHi[d] = bc->face(d, 1);
        }
      }
    }
  }
};

#endif

#ifndef _SHAPE_H_
#define _SHAPE_H_

#include <string>

#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_RealBox.H>

// Abstract base class for shapes. All shapes should be defined in normalized
// units
class Shape {
public:
  virtual ~Shape();
  virtual bool is_inside(const amrex::Real* xyz) const = 0;

  std::string get_name() const { return name; };

private:
  std::string name;
};

//=========== Box ===========//
class BoxShape : public Shape {
public:
  BoxShape(const amrex::Real* lo, const amrex::Real* hi) {
    rb.setLo(lo);
    rb.setHi(hi);
  }

  BoxShape(amrex::RealBox& in) { rb = in; };

  bool is_inside(const amrex::Real* xyz) const override {
    return rb.contains(xyz);
  };

private:
  amrex::RealBox rb;
};

//========== Sphere =========//
class Sphere : public Shape {
public:
  Sphere(const amrex::Real* c, amrex::Real r) {
    radius = r;

    for (int i = 0; i < 3; ++i)
      center[i] = c[i];
  }

  bool is_inside(const amrex::Real* xyz) const override {
    amrex::Real l =
        sqrt(pow(xyz[0] - center[0], 2) + pow(xyz[1] - center[1], 2) +
             pow(xyz[2] - center[2], 2));
    return l < radius;
  };

private:
  amrex::Real radius;
  amrex::Real center[3];
};

//========= Shell ==========//
class Shell : public Shape {
  Shell(const amrex::Real* c, amrex::Real rInner, amrex::Real rOuter) {
    radiusInner = rInner;
    radiusOuter = rOuter;

    for (int i = 0; i < 3; ++i)
      center[i] = c[i];
  }

  bool is_inside(const amrex::Real* xyz) const override {
    amrex::Real l =
        sqrt(pow(xyz[0] - center[0], 2) + pow(xyz[1] - center[1], 2) +
             pow(xyz[2] - center[2], 2));
    return l < radiusOuter && l > radiusInner;
  };

private:
  amrex::Real radiusInner, radiusOuter;
  amrex::Real center[3];
};

#endif
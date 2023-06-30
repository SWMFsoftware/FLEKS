#ifndef _SHAPE_H_
#define _SHAPE_H_

#include <string>

#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_RealBox.H>

#include "Constants.h"

// Abstract base class for shapes. All shapes should be defined in normalized
// units
class Shape {
public:
  virtual ~Shape(){};
  virtual bool is_inside(const amrex::Real* xyz) const = 0;

  std::string get_name() const { return name; };

protected:
  std::string name;
};

//=========== Box ===========//
class BoxShape : public Shape {
public:
  BoxShape(const std::string& nameIn, const amrex::Real* lo,
           const amrex::Real* hi) {
    name = nameIn;
    rb.setLo(lo);
    rb.setHi(hi);
  }

  BoxShape(const std::string& nameIn, amrex::RealBox& in) {
    name = nameIn;
    rb = in;
  };

  bool is_inside(const amrex::Real* xyz) const override {
    return rb.contains(xyz);
  };

private:
  amrex::RealBox rb;
};

//========== Sphere =========//
class Sphere : public Shape {
public:
  Sphere(const std::string& nameIn, const amrex::Real* c, amrex::Real r) {
    name = nameIn;

    radius = r;
    for (int i = 0; i < nDim; ++i)
      center[i] = c[i];
  }

  bool is_inside(const amrex::Real* xyz) const override {
    amrex::Real l2 = 0;

    for (int i = 0; i < nDim; ++i)
      l2 += pow(xyz[i] - center[i], 2);
    
    return l2 < pow(radius, 2);
  };

private:
  amrex::Real radius;
  amrex::Real center[nDim];
};

//========= Shell ==========//
class Shell : public Shape {
public:
  Shell(const std::string& nameIn, const amrex::Real* c, amrex::Real rInner,
        amrex::Real rOuter) {
    name = nameIn;

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
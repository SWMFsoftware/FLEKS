#ifndef _BC_H_
#define _BC_H_

#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_iMultiFab.H>
#include <map>

void apply_float_boundary(const amrex::iMultiFab& status, amrex::MultiFab& mf,
                          const amrex::Geometry& gm, const int iStart,
                          const int nComp, const int nshift = 0);

class BC {
public:
  static const int unset = -1;
  static const int periodic = 0;
  static const int coupled = 1;
  static const int outflow = 2;
  static const int vacume = 3;

  amrex::IntVect lo;
  amrex::IntVect hi;

  BC() : lo(coupled), hi(coupled) {}

  int num_type(const std::string& str) {
    if (str == "periodic")
      return periodic;
    else if (str == "coupled")
      return coupled;
    else if (str == "outflow")
      return outflow;
    else if (str == "vacume")
      return vacume;
    else
      return unset;
  }
};

#endif

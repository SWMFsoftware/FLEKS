#ifndef _BC_H_
#define _BC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>

void apply_float_boundary(amrex::MultiFab& mf, const amrex::Geometry& geom,
                          const int iStart, const int nComp,                          
                          const int nshift = 0);
#endif
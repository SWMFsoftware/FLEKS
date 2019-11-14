#ifndef _BC_H_
#define _BC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>



void apply_zero_boundary(amrex::MultiFab& mf, const amrex::Geometry& geom);
void zero_boundary_cpu(amrex::Box const& bx,
                       amrex::Array4<amrex::Real> const& arr, const int iStart,
                       const int nComp, amrex::GeometryData const& geom,
                       const amrex::Real time, const amrex::BCRec* bcr,
                       const int bcomp, const int orig_comp);
// void apply_float_boundary(amrex::MultiFab& mf, const amrex::Geometry& geom);
// void apply_external_boundary();

#endif
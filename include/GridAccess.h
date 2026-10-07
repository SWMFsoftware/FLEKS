#ifndef _GRID_ACCESS_H_
#define _GRID_ACCESS_H_

#include "Grid.h"

// Non-owning query facade over the one simulation Grid that Domain owns.
//
// Pic, FluidInterface, and ParticleTracker all need the same mesh queries.  The
// definitions live here once instead of being repeated in every consumer.  The
// reference is stored once in the base, so a consumer must not declare a second
// Grid member.
//
// GridAccess exposes no mesh mutation, no ownership, and no virtual hooks: the
// hierarchy is changed only through Domain's orchestration path.
class GridAccess {
protected:
  Grid& grid;

  explicit GridAccess(Grid& gridIn) : grid(gridIn) {}
  ~GridAccess() = default;

public:
  GridAccess(const GridAccess&) = delete;
  GridAccess& operator=(const GridAccess&) = delete;
  GridAccess(GridAccess&&) = delete;
  GridAccess& operator=(GridAccess&&) = delete;

  int n_lev() const { return grid.n_lev(); }
  int n_lev_max() const { return grid.n_lev_max(); }

  // The queries below are inherited from amrex::AmrCore/AmrMesh and keep their
  // original library names on purpose; see the "AMReX API exception" note in
  // doc/Coding_standards.md. Grid-owned queries elsewhere in this facade keep
  // the FLEKS snake_case spelling (cell_status, node_box_array, ...).
  int finestLevel() const { return grid.finestLevel(); }
  const amrex::Geometry& Geom(int iLev) const { return grid.Geom(iLev); }
  const amrex::Vector<amrex::Geometry>& Geom() const { return grid.Geom(); }
  // No-argument overload; returns all refinement ratios by value. AmrMesh
  // returns a const vector reference here; the copy is kept so this facade's
  // established return semantics do not change.
  amrex::Vector<amrex::IntVect> refRatio() const { return grid.refRatio(); }
  // Level-specific overload; returns the IntVect by value, matching AmrMesh.
  amrex::IntVect refRatio(int iLev) const { return grid.refRatio(iLev); }
  int find_mpi_rank_from_coord(const amrex::RealVect& xyz) const {
    return grid.find_mpi_rank_from_coord(xyz);
  }
  const amrex::BoxArray& boxArray(int iLev) const {
    return grid.boxArray(iLev);
  }
  const amrex::BoxArray& active_region_ref() const {
    return grid.active_region_ref();
  }
  const amrex::BoxArray& node_box_array(int iLev) const {
    return grid.node_box_array(iLev);
  }
  const amrex::Vector<amrex::BoxArray>& node_box_arrays() const {
    return grid.node_box_arrays();
  }
  const amrex::DistributionMapping& DistributionMap(int iLev) const {
    return grid.DistributionMap(iLev);
  }
  const amrex::Vector<amrex::iMultiFab>& cell_status() const {
    return grid.cell_status();
  }
  const amrex::iMultiFab& cell_status(int iLev) const {
    return grid.cell_status(iLev);
  }
  const amrex::iMultiFab& node_status(int iLev) const {
    return grid.node_status(iLev);
  }
  bool is_grid_empty() const { return grid.is_grid_empty(); }
  bool is_new_grid() const { return grid.is_new_grid(); }
  bool is_fake_2d() const { return grid.is_fake_2d(); }
  int get_n_ghost() const { return grid.get_n_ghost(); }
  const amrex::AmrInfo& get_amr_info() const { return grid.get_amr_info(); }
  const RefineRegions* get_refine_regions() const {
    return grid.get_refine_regions();
  }
  amrex::BoxArray get_base_grid() const { return grid.get_base_grid(); }
  std::string lev_string(int iLev) const { return grid.lev_string(iLev); }
  bool is_inside_domain(const amrex::Real* loc) const {
    return grid.is_inside_domain(loc);
  }
  const amrex::Vector<amrex::RealBox>& domain_range() const {
    return grid.domain_range();
  }
  bool use_body() const { return grid.use_body(); }
  amrex::Real get_body_radius() const { return grid.get_body_radius(); }
  const amrex::Real* get_body_center() const { return grid.get_body_center(); }
  int get_dim() const { return grid.get_dim(); }
  const amrex::Vector<amrex::MultiFab>& get_cost() const {
    return grid.get_cost();
  }
  amrex::Real get_cell_volume(int iLev) const {
    return grid.get_cell_volume(iLev);
  }
  bool is_inside_body(const amrex::Real* loc) const {
    return grid.is_inside_body(loc);
  }
  int get_finest_lev(const amrex::RealVect& xyz) const {
    return grid.get_finest_lev(xyz);
  }
};

#endif

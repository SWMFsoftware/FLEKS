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
  int finestLevel() const { return grid.finestLevel(); }
  const amrex::Geometry& Geom(int iLev) const { return grid.Geom(iLev); }
  const amrex::BoxArray& boxArray(int iLev) const {
    return grid.box_array(iLev);
  }
  const amrex::BoxArray& node_box_array(int iLev) const {
    return grid.node_box_array(iLev);
  }
  const amrex::Vector<amrex::BoxArray>& node_box_arrays() const {
    return grid.node_box_arrays();
  }
  const amrex::DistributionMapping& DistributionMap(int iLev) const {
    return grid.get_dmap(iLev);
  }
  const amrex::iMultiFab& cell_status(int iLev) const {
    return grid.cell_status(iLev);
  }
  const amrex::iMultiFab& node_status(int iLev) const {
    return grid.node_status(iLev);
  }
  bool is_grid_empty() const { return grid.is_grid_empty(); }
  bool is_new_grid() const { return grid.is_new_grid(); }
  int get_n_ghost() const { return grid.get_n_ghost(); }
  const amrex::AmrInfo& get_amr_info() const { return grid.get_amr_info(); }
  const RefineRegions* get_refine_regions() const {
    return grid.get_refine_regions();
  }
  amrex::BoxArray get_base_grid() const { return grid.get_base_grid(); }
  std::string lev_string(int iLev) const { return grid.lev_string(iLev); }
  int get_finest_lev(const amrex::RealVect& xyz) const {
    return grid.get_finest_lev(xyz);
  }
};

#endif

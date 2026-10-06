#ifndef _MESH_CHANGE_H_
#define _MESH_CHANGE_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Vector.H>

enum class MeshChangeReason { None, Init, Topology, LoadBalance, Restart };

struct MeshChange {
  MeshChangeReason reason = MeshChangeReason::None;
  int oldFinestLevel = 0;
  amrex::Vector<amrex::BoxArray> oldCellGrids;
  amrex::Vector<amrex::BoxArray> oldNodeGrids;
  amrex::Vector<amrex::DistributionMapping> oldDmaps;

  bool is_topology_change() const {
    return reason == MeshChangeReason::Topology ||
           reason == MeshChangeReason::Init ||
           reason == MeshChangeReason::Restart;
  }

  bool is_load_balance() const {
    return reason == MeshChangeReason::LoadBalance;
  }
};

#endif

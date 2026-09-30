#ifndef _REFINEREGIONS_H_
#define _REFINEREGIONS_H_

#include <AMReX_Vector.H>

#include "Regions.h"

// Refinement regions and whether their definitions changed since the last
// regrid. Call mark_modified() after changing a contained Regions object.
class RefineRegions : public amrex::Vector<Regions> {
public:
  using amrex::Vector<Regions>::Vector;

  bool is_modified() const { return modified; }
  void mark_modified() { modified = true; }
  void clear_modified() { modified = false; }

private:
  bool modified = false;
};

#endif

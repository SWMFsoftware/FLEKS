#ifndef _DOMAINGRID_H_
#define _DOMAINGRID_H_

#include <memory>
#include <set>
#include <sstream>
#include <string>

#include <AMReX_AmrMesh.H>
#include <AMReX_BCRec.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include "Constants.h"
#include "FleksDistributionMap.h"
#include "GridInfo.h"
#include "RefineRegions.h"

class DomainGrid {

protected:
  const int nGst = 2;

  bool isFake2D = false;

  amrex::Vector<int> nCell = { 1, 1, 1 };

  amrex::IntVect maxBlockSize;
  amrex::IntVect periodicity;
  amrex::IntVect centerBoxLo;
  amrex::IntVect centerBoxHi;
  amrex::Box centerBox;
  amrex::RealBox domainRange;

  amrex::Real lNormSI = 0;
  amrex::Real uNormSI = 0;
  amrex::Real mNormSI = 0;
  int scalingFactor = 1;

  const int coord = 0; // Cartesian grid
  amrex::Geometry gm;

  amrex::AmrInfo amrInfo;

  GridInfo gridInfo;

  int iGrid = 1;
  int iDecomp = 1;

  int gridID;
  std::string printPrefix;
  std::string gridName;
  int refinementRatio = 2;

  // Shapes persist across parameter sessions. A selector may start using a
  // previously defined shape, but a name never acquires new geometry.
  amrex::Vector<std::shared_ptr<Shape> > shapeList;
  std::set<std::string> shapeNames;
  amrex::Vector<std::string> refineRegionsStr;
  RefineRegions refineRegions;

  void add_shape(const std::shared_ptr<Shape>& shape) {
    const std::string name = shape->get_name();
    if (!shapeNames.insert(name).second)
      amrex::Abort("Duplicate #REGION name: " + name);
    shapeList.push_back(shape);
  }

  bool set_refine_region(int iLev, const std::string& selector) {
    if (iLev < 0 || iLev >= static_cast<int>(refineRegionsStr.size()) - 1)
      amrex::Abort(
          "Invalid refinement level " + std::to_string(iLev) +
          " in #REFINEREGION: max allowed level is " +
          std::to_string(static_cast<int>(refineRegionsStr.size()) - 2));

    std::stringstream input(selector);
    std::string token, normalized;
    bool hasNone = false;
    int tokenCount = 0;
    while (input >> token) {
      tokenCount++;
      if (token == "none")
        hasNone = true;
      if (!normalized.empty())
        normalized += ' ';
      normalized += token;
    }
    if (hasNone && tokenCount > 1)
      amrex::Abort("Cannot combine 'none' with other regions in #REFINEREGION");

    // An empty selector tags no cells; "none" is the input spelling for it.
    if (normalized == "none")
      normalized.clear();

    // Repeating a selector does not request another AMR rebuild.
    if (refineRegionsStr[iLev] == normalized)
      return false;
    refineRegionsStr[iLev] = normalized;
    return true;
  }

  void rebuild_refine_regions() {
    for (int i = 0; i < static_cast<int>(refineRegionsStr.size()); ++i) {
      std::stringstream input(refineRegionsStr[i]);
      std::string token;
      while (input >> token) {
        if (token.size() < 2 || (token[0] != '+' && token[0] != '-'))
          amrex::Abort("Invalid shape prefix in #REFINEREGION: '" + token +
                       "' (must start with '+' or '-')");
        if (shapeNames.count(token.substr(1)) == 0)
          amrex::Abort("Unknown shape in #REFINEREGION: '" + token.substr(1) +
                       "'");
      }
      if (i < static_cast<int>(refineRegions.size()))
        refineRegions[i].define(shapeList, refineRegionsStr[i]);
    }
  }

  // "This threshold value, which defaults to 0.7 (or 70%), is used to ensure
  // that grids do not contain too large a fraction of un-tagged cells." - AMReX
  // online docs
  amrex::Real gridEfficiency = 0.7;

  // If the grid has not been initialized or the grid changed due to AMR,
  // isNewGrid is true.
  bool isNewGrid = true;

  bool doSplitLevs = false;

public:
  DomainGrid() {
    for (int i = 0; i < nDim; ++i) {
      periodicity[i] = 0;
      maxBlockSize[i] = 8;
    }
  }
  ~DomainGrid() = default;

  int get_iGrid() const { return iGrid; }
  int get_iDecomp() const { return iDecomp; }
  int get_dim() const { return (isFake2D || nDim == 2) ? 2 : nDim; }
  void set_periodicity(const int iDir, const bool isPeriodic) {
    periodicity[iDir] = (isPeriodic ? 1 : 0);
  }
};
#endif

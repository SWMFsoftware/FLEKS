#ifndef _DOMAINGRID_H_
#define _DOMAINGRID_H_

#include <map>
#include <memory>
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

  amrex::Vector<std::shared_ptr<Shape> > shapeList;
  std::map<std::string, std::string> shapeSignatures;
  amrex::Vector<std::string> refineRegionsStr;
  RefineRegions refineRegions;

  bool upsert_shape(const std::shared_ptr<Shape>& shape,
                    const std::string& signature) {
    const std::string name = shape->get_name();
    const auto signatureIt = shapeSignatures.find(name);
    if (signatureIt != shapeSignatures.end() &&
        signatureIt->second == signature)
      return false;

    for (auto& existing : shapeList) {
      if (existing->get_name() == name) {
        existing = shape;
        shapeSignatures[name] = signature;
        return true;
      }
    }

    shapeList.push_back(shape);
    shapeSignatures[name] = signature;
    return true;
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

    if (normalized == "none")
      normalized.clear();

    if (refineRegionsStr[iLev] == normalized)
      return false;
    refineRegionsStr[iLev] = normalized;
    return true;
  }

  bool is_shape_used_for_refinement(const std::string& name) const {
    for (const auto& selector : refineRegionsStr) {
      std::stringstream input(selector);
      std::string token;
      while (input >> token) {
        if (token.size() > 1 && token.substr(1) == name)
          return true;
      }
    }
    return false;
  }

  void rebuild_refine_regions() {
    for (int i = 0; i < static_cast<int>(refineRegionsStr.size()); ++i) {
      std::stringstream input(refineRegionsStr[i]);
      std::string token;
      while (input >> token) {
        if (token.size() < 2 || (token[0] != '+' && token[0] != '-'))
          amrex::Abort("Invalid shape prefix in #REFINEREGION: '" + token +
                       "' (must start with '+' or '-')");
        if (shapeSignatures.count(token.substr(1)) == 0)
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
  void set_periodicity(const int iDir, const bool isPeriodic) {
    periodicity[iDir] = (isPeriodic ? 1 : 0);
  }
};
#endif

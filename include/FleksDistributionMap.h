#ifndef _FLEKSDISTRIBUTIONMAP_H_
#define _FLEKSDISTRIBUTIONMAP_H_

#include <AMReX_DistributionMapping.H>

class FleksDistributionMap : public amrex::DistributionMapping {

  // This method is the same as DistributionMapping::makeKnapSack(), except that
  // 'sort' is provided as an argument below.
  static DistributionMapping make_knapsack_for_fleks(
      const amrex::MultiFab& weight, amrex::Real& eff,
      int nmax = std::numeric_limits<int>::max(), bool sort = false);

  static DistributionMapping make_knapsack_for_fleks(
      const amrex::MultiFab& weight, const amrex::Vector<int>& ord,
      amrex::Real& eff, int nmax = std::numeric_limits<int>::max(),
      bool sort = false);
};

#endif // _FLEKSDISTRIBUTIONMAP_H_
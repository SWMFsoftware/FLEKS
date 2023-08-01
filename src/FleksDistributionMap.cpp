
#include <AMReX_MultiFab.H>

#include "FleksDistributionMap.h"

using namespace amrex;

namespace {
Vector<Long> gather_weights_fleks(const MultiFab& weight) {
  BL_PROFILE("gather_weights_fleks");

  LayoutData<Real> costld(weight.boxArray(), weight.DistributionMap());
  for (MFIter mfi(weight); mfi.isValid(); ++mfi) {
    costld[mfi] = weight[mfi].sum<RunOn::Device>(mfi.validbox(), 0);
  }
  Vector<Real> rcost(weight.size());
  ParallelDescriptor::GatherLayoutDataToVector(
      costld, rcost, ParallelContext::IOProcessorNumberSub());
  ParallelDescriptor::Bcast(rcost.data(), rcost.size(),
                            ParallelContext::IOProcessorNumberSub());
  Real wmax = *std::max_element(rcost.begin(), rcost.end());
  Real scale = (wmax == 0) ? 1.e9_rt : 1.e9_rt / wmax;
  Vector<Long> lcost(rcost.size());
  for (int i = 0; i < rcost.size(); ++i) {
    lcost[i] = static_cast<Long>(rcost[i] * scale) + 1L;
  }
  return lcost;
}
} // namespace

DistributionMapping FleksDistributionMap::make_knapsack_for_fleks(
    const MultiFab& weight, Real& eff, int nmax, bool sort) {
  BL_PROFILE("make_knapsack_for_fleks_v1");

  Vector<Long> cost = gather_weights_fleks(weight);

  int nprocs = ParallelContext::NProcsSub();
  DistributionMapping r;
  r.KnapSackProcessorMap(cost, nprocs, &eff, sort, nmax);
  return r;
}

DistributionMapping FleksDistributionMap::make_knapsack_for_fleks(
    const MultiFab& weight, const Vector<int>& ord, Real& eff, int nmax,
    bool sort) {
  BL_PROFILE("make_knapsack_for_fleks_v2");

  if (ord.size() != ParallelDescriptor::NProcs()) {
    amrex::Abort("ord.size()!=ParallelDescriptor::NProcs()");
  }

  DistributionMapping dm = make_knapsack_for_fleks(weight, eff, nmax, sort);

  Vector<int> pmapNew;
  pmapNew.resize(dm.size());

  const Vector<int>& pmap = dm.ProcessorMap();

  for (long i = 0; i < pmap.size(); ++i) {
    pmapNew[i] = ord[pmap[i]];
  }

  return DistributionMapping(pmapNew);
}
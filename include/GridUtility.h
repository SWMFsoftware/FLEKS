#ifndef _GRIDUTILITY_H_
#define _GRIDUTILITY_H_

#include <AMReX_DistributionMapping.H>
#include <AMReX_FabArray.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PhysBCFunct.H>

#include "Bit.h"
#include "Constants.h"
#include "GridInfo.h"
#include "UInterp.h"
#include "Utility.h"

void curl_center_to_node(const amrex::MultiFab& centerMF,
                         amrex::MultiFab& nodeMF, const amrex::Real* invDx);

void curl_node_to_center(const amrex::MultiFab& nodeMF,
                         amrex::MultiFab& centerMF, const amrex::Real* invDx);

void lap_node_to_node(const amrex::MultiFab& srcMF, amrex::MultiFab& dstMF,
                      const amrex::DistributionMapping dm,
                      const amrex::Geometry& gm);

void grad_node_to_center(const amrex::MultiFab& nodeMF,
                         amrex::MultiFab& centerMF, const amrex::Real* invDx);

void grad_center_to_node(const amrex::MultiFab& centerMF,
                         amrex::MultiFab& nodeMF, const amrex::Real* invDx);

void div_center_to_node(const amrex::MultiFab& centerMF,
                        amrex::MultiFab& nodeMF, const amrex::Real* invDx);

void div_node_to_center(const amrex::MultiFab& nodeMF,
                        amrex::MultiFab& centerMF, const amrex::Real* invDx);

void div_center_to_center(const amrex::MultiFab& srcMF, amrex::MultiFab& dstMF,
                          const amrex::Real* invDx);

void average_center_to_node(const amrex::MultiFab& centerMF,
                            amrex::MultiFab& nodeMF);

void print_MultiFab(const amrex::iMultiFab& data, std::string tag,
                    int nshift = 0);

void print_MultiFab(const amrex::MultiFab& data, std::string tag,
                    const int iVarStart, const int iVarEnd, int nshift = 0);

void print_MultiFab(const amrex::MultiFab& data, std::string tag,
                    amrex::Geometry& gm, int nshift = 0);

template <class FAB>
void print_fab(const amrex::FabArray<FAB>& mf, std::string tag,
               const int iStart, const int nComp, int nshift = 0) {
  amrex::AllPrint() << "-----" << tag << " begin-----" << std::endl;
  amrex::Real sum = 0;
  amrex::Real sum2 = 0;

  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const FAB& fab = mf[mfi];
    const auto& box = mfi.validbox();
    const auto& data = fab.array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    amrex::AllPrint() << "------ box = " << box << std::endl;
    for (int i = lo.x - nshift; i <= hi.x + nshift; ++i)
      for (int j = lo.y - nshift; j <= hi.y + nshift; ++j)
        // for (int k = lo.z; k <= hi.z; ++k)
        for (int iVar = iStart; iVar < iStart + nComp; iVar++) {
          int k = 0;
          amrex::AllPrint() << " i = " << i << " j = " << j << " k = " << k
                            << " iVar = " << iVar
                            << " data = " << data(i, j, k, iVar) << std::endl;
          sum += data(i, j, k, iVar);
          sum2 += pow(data(i, j, k, iVar), 2);
        }
  }
  amrex::AllPrint() << "sum = " << sum << " sum2 = " << sqrt(sum2)
                    << " on proc = " << amrex::ParallelDescriptor::MyProc()
                    << std::endl;
  amrex::AllPrint() << "-----" << tag << " end-----" << std::endl;
}

inline int get_local_node_or_cell_number(const amrex::MultiFab& MF) {
  int nTotal = 0;
  if (!MF.empty())
    for (amrex::MFIter mfi(MF); mfi.isValid(); ++mfi) {
      const amrex::Box& box = mfi.validbox();
      const auto lo = lbound(box);
      const auto hi = ubound(box);
      nTotal += (hi.x - lo.x + 1) * (hi.y - lo.y + 1) * (hi.z - lo.z + 1);
    }
  return nTotal;
}

inline void find_node_index(const amrex::RealVect& xyz,
                            const amrex::Real* const plo,
                            const amrex::Real* const invDx,
                            amrex::IntVect& loIdx, amrex::RealVect& dShift) {
  for (int i = 0; i < nDim; ++i) {
    dShift[i] = (xyz[i] - plo[i]) * invDx[i];
    loIdx[i] = fastfloor(dShift[i]);
    dShift[i] = dShift[i] - loIdx[i];
  }
}

inline void find_cell_index(const amrex::RealVect& xyz,
                            const amrex::Real* const plo,
                            const amrex::Real* const invDx,
                            amrex::IntVect& loIdx, amrex::RealVect& dShift) {
  for (int i = 0; i < nDim; ++i) {
    // plo is the corner location => -0.5
    dShift[i] = (xyz[i] - plo[i]) * invDx[i] - 0.5;
    loIdx[i] = fastfloor(dShift[i]);
    dShift[i] = dShift[i] - loIdx[i];
  }
}

inline amrex::Real get_value_at_loc(const amrex::MultiFab& mf,
                                    const amrex::MFIter& mfi,
                                    const amrex::Geometry& gm,
                                    const amrex::RealVect xyz, const int iVar) {
  amrex::IntVect loIdx;
  amrex::RealVect dx;
  find_node_index(xyz, gm.ProbLo(), gm.InvCellSize(), loIdx, dx);

  amrex::Real coef[2][2][2];
  {
    amrex::Real xi[2];
    amrex::Real eta[2];
    amrex::Real zeta[2];
    xi[0] = dx[0];
    eta[0] = dx[1];
    zeta[0] = nDim > 2 ? dx[2] : 1;
    xi[1] = 1 - xi[0];
    eta[1] = 1 - eta[0];
    zeta[1] = nDim > 2 ? 1 - zeta[0] : 1;

    amrex::Real multi[2][2];
    multi[0][0] = xi[0] * eta[0];
    multi[0][1] = xi[0] * eta[1];
    multi[1][0] = xi[1] * eta[0];
    multi[1][1] = xi[1] * eta[1];

    // coef[k][j][i]: This is so wired. But is may be faster since it matches
    // the AMREX multifab data ordering.
    if (nDim == 2) {
      coef[0][0][0] = multi[1][1];
      coef[0][0][1] = multi[0][1];
      coef[0][1][0] = multi[1][0];
      coef[0][1][1] = multi[0][0];
    } else {
      coef[0][0][0] = multi[1][1] * zeta[1];
      coef[1][0][0] = multi[1][1] * zeta[0];
      coef[0][1][0] = multi[1][0] * zeta[1];
      coef[1][1][0] = multi[1][0] * zeta[0];
      coef[0][0][1] = multi[0][1] * zeta[1];
      coef[1][0][1] = multi[0][1] * zeta[0];
      coef[0][1][1] = multi[0][0] * zeta[1];
      coef[1][1][1] = multi[0][0] * zeta[0];
    }
  }

  const auto& arr = mf[mfi].array();
  amrex::Real val = 0;

  amrex::Box box = amrex::Box(amrex::IntVect(0), amrex::IntVect(1));

  amrex::ParallelFor(box, [&](int ii, int jj, int kk) noexcept {
    int iIdx = loIdx[ix_] + ii;
    int jIdx = loIdx[iy_] + jj;
    int kIdx = nDim > 2 ? loIdx[iz_] + kk : 0;
    val += arr(iIdx, jIdx, kIdx, iVar) * coef[kk][jj][ii];
  });

  return val;
}

inline amrex::Real get_value_at_loc(const amrex::MultiFab& mf,
                                    const amrex::Geometry& gm,
                                    const amrex::RealVect xyz, const int iVar) {
  auto idx = gm.CellIndex(xyz.begin());

  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    // Cell box
    const amrex::Box& bx =
        amrex::convert(mfi.validbox(), { AMREX_D_DECL(0, 0, 0) });
    if (bx.contains(idx))
      return get_value_at_loc(mf, mfi, gm, xyz, iVar);
  }

  amrex::AllPrint() << "xyz = " << xyz << std::endl;

  amrex::Abort("Error: can not find this point!");
  return -1; // To suppress compiler warnings.
}

inline void add_to_mf(const amrex::Real& val, amrex::MultiFab& mf,
                      const amrex::MFIter& mfi, const amrex::Geometry& gm,
                      const amrex::RealVect xyz, const int iVar) {
  const auto plo = gm.ProbLo();

  const auto invDx = gm.InvCellSize();

  int loIdx[3];
  amrex::Real dx[3] = { 0, 0, 0 };
  for (int i = 0; i < nDim; ++i) {
    dx[i] = (xyz[i] - plo[i]) * invDx[i];
    loIdx[i] = fastfloor(dx[i]);
    dx[i] = dx[i] - loIdx[i];
  }

  amrex::Real coef[2][2][2];
  {
    amrex::Real xi[2];
    amrex::Real eta[2];
    amrex::Real zeta[2];
    xi[0] = dx[0];
    eta[0] = dx[1];
    zeta[0] = dx[2];
    xi[1] = 1 - xi[0];
    eta[1] = 1 - eta[0];
    zeta[1] = 1 - zeta[0];

    amrex::Real multi[2][2];
    multi[0][0] = xi[0] * eta[0];
    multi[0][1] = xi[0] * eta[1];
    multi[1][0] = xi[1] * eta[0];
    multi[1][1] = xi[1] * eta[1];

    // coef[k][j][i]
    coef[0][0][0] = multi[1][1] * zeta[1];
    coef[1][0][0] = multi[1][1] * zeta[0];
    coef[0][1][0] = multi[1][0] * zeta[1];
    coef[1][1][0] = multi[1][0] * zeta[0];
    coef[0][0][1] = multi[0][1] * zeta[1];
    coef[1][0][1] = multi[0][1] * zeta[0];
    coef[0][1][1] = multi[0][0] * zeta[1];
    coef[1][1][1] = multi[0][0] * zeta[0];
  }

  const amrex::Array4<amrex::Real>& arr = mf[mfi].array();
  for (int kk = 0; kk < 2; ++kk) {
    const int kIdx = loIdx[iz_] + kk;
    for (int jj = 0; jj < 2; ++jj) {
      const int jIdx = loIdx[iy_] + jj;
      for (int ii = 0; ii < 2; ++ii) {
        arr(loIdx[ix_] + ii, jIdx, kIdx, iVar) += val * coef[kk][jj][ii];
      }
    }
  }

  return;
}

template <class FAB> class PhysBCFunctNoOpFab {
public:
  void operator()(amrex::FabArray<FAB>& /*mf*/, int /*dcomp*/, int /*ncomp*/,
                  amrex::IntVect const& /*nghost*/, amrex::Real /*time*/,
                  int /*bccomp*/) {}
};

template <class FAB>
void distribute_FabArray(amrex::FabArray<FAB>& fa, amrex::BoxArray baNew,
                         const amrex::DistributionMapping& dm, int nComp,
                         int nGst, bool doCopy = true) {
  // Assume 'dm' is the new dm.
  amrex::FabArray<FAB> tmp;

  if (!baNew.empty()) {
    tmp.define(baNew, dm, nComp, nGst);
    tmp.setVal(-7777);
    if (doCopy && !fa.empty()) {
      tmp.ParallelCopy(fa, 0, 0, fa.nComp(), nGst, nGst);
      // tmp.Redistribute(fa, 0, 0, tmp.nComp(), tmp.nGrowVect());
    } else {
      tmp.setVal(0);
    }
  }

  fa = std::move(tmp);
}

template <class FAB>
void distribute_FabArray(amrex::FabArray<FAB>& fa, amrex::BoxArray baNew,
                         const amrex::DistributionMapping& dm,
                         bool doCopy = true) {

  distribute_FabArray(fa, baNew, dm, fa.nComp(), fa.nGrow(), doCopy);
}

// Interpolate from coarse lev to fine lev
template <class FAB, class Interp>
void interp_from_coarse_to_fine(amrex::FabArray<FAB>& coarse,
                                amrex::FabArray<FAB>& fine, const int iStart,
                                const int nComp, const amrex::IntVect ratio,
                                const amrex::Geometry& cgeom,
                                const amrex::Geometry& fgeom, Interp* mapper,
                                int nGst = 0) {
  BL_PROFILE("interp_from_coarse_to_fine");
  amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);

  int lo_bc[3] = { amrex::BCType::int_dir, amrex::BCType::int_dir,
                   amrex::BCType::int_dir };
  int hi_bc[3] = { amrex::BCType::int_dir, amrex::BCType::int_dir,
                   amrex::BCType::int_dir };

  amrex::Vector<amrex::BCRec> bcs(nComp, amrex::BCRec(lo_bc, hi_bc));

  PhysBCFunctNoOpFab<FAB> cphysbc, fphysbc;

  // See definition in
  // AMREX/InstallDir/include/AMReX_FillPatchUtil_I.H
  amrex::InterpFromCoarseLevel(f, amrex::IntVect(nGst), 0.0, c, 0, 0, nComp,
                               cgeom, fgeom, cphysbc, 0, fphysbc, 0, ratio,
                               mapper, bcs, 0);
}

// Sum from fine level to coarse level for nodes at the interface of two levels.
template <class FAB>
void sum_fine_to_coarse_lev_bny_node(amrex::FabArray<FAB>& coarse,
                                     amrex::FabArray<FAB>& fine,
                                     const int iStart, const int nComp,
                                     const amrex::IntVect ratio) {
  BL_PROFILE("sum_fine_to_coarse_lev_bny_node");
  amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);

  amrex::FabArray<FAB> ctmp(c.boxArray(), c.DistributionMap(), nComp, 0);
  ctmp.setVal(0.0);

  // // f and c/ctmp have different distribution maps. Inside
  // average_down_nodal(), parallel copy is called to copy coarsen(f) to ctmp.
  amrex::average_down_nodal(f, ctmp, ratio);

  // This is the Add() declared in AMReX_FabArray.H, and it works for any FAB.
  amrex::Add(c, ctmp, 0, 0, nComp, 0);
}

// Sum from coarse level to fine level for nodes at the boundary of two levels.

template <class FAB>
void fill_lev_bny_from_value(amrex::FabArray<FAB>& dst,
                             const amrex::iMultiFab& fstatus,
                             amrex::Real value) {

  for (amrex::MFIter mfi(dst); mfi.isValid(); ++mfi) {
    FAB& fab = dst[mfi];
    const auto& box = mfi.fabbox();
    const auto& data = fab.array();

    const auto& statusArr = fstatus[mfi].array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int iVar = 0; iVar < dst.nComp(); iVar++)
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (bit::is_lev_boundary(statusArr(i, j, k))) {
              data(i, j, k, iVar) = value;
            }
          }
  }
}

template <class FAB>
void fill_lev_from_value(amrex::FabArray<FAB>& dst, amrex::Real value,
                         int startvar = 0, int stopvar = -1) {
  if (stopvar == -1) {
    stopvar = dst.nComp() - 1;
  }
  for (amrex::MFIter mfi(dst); mfi.isValid(); ++mfi) {
    FAB& fab = dst[mfi];
    const auto& box = mfi.fabbox();
    const auto& data = fab.array();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int iVar = startvar; iVar <= stopvar; iVar++)
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            data(i, j, k, iVar) = value;
          }
  }
}

template <class FAB>
void fill_fine_lev_bny_from_coarse(amrex::FabArray<FAB>& coarse,
                                   amrex::FabArray<FAB>& fine, const int iStart,
                                   const int nComp, const amrex::IntVect ratio,
                                   const amrex::Geometry& cgeom,
                                   const amrex::Geometry& fgeom,
                                   const amrex::iMultiFab& fstatus,
                                   amrex::Interpolater& mapper) {
  BL_PROFILE("fill_fine_lev_bny_from_coarse");

  amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);

  amrex::FabArray<FAB> ftmp(f.boxArray(), f.DistributionMap(), nComp,
                            fine.nGrow());
  ftmp.setVal(0.0);

  interp_from_coarse_to_fine(c, ftmp, 0, nComp, ratio, cgeom, fgeom, &mapper,
                             f.nGrow());

  for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
    FAB& fab = f[mfi];
    const auto& box = mfi.fabbox();
    const auto& data = fab.array();

    const auto& statusArr = fstatus[mfi].array();
    const auto& tmp = ftmp[mfi].array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int iVar = 0; iVar < f.nComp(); iVar++)
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (bit::is_lev_boundary(statusArr(i, j, k))) {
              data(i, j, k, iVar) = tmp(i, j, k, iVar);
            }
          }
  }
}

template <class FAB>
void fill_fine_lev_edge_from_coarse(
    amrex::FabArray<FAB>& coarse, amrex::FabArray<FAB>& fine, const int iStart,
    const int nComp, const amrex::IntVect ratio, const amrex::Geometry& cgeom,
    const amrex::Geometry& fgeom, const amrex::iMultiFab& fstatus,
    amrex::Interpolater& mapper) {
  BL_PROFILE("fill_fine_lev_bny_from_coarse");

  amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);

  amrex::FabArray<FAB> ftmp(f.boxArray(), f.DistributionMap(), nComp,
                            fine.nGrow());
  ftmp.setVal(0.0);

  interp_from_coarse_to_fine(c, ftmp, 0, nComp, ratio, cgeom, fgeom, &mapper,
                             f.nGrow());

  for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
    FAB& fab = f[mfi];
    const auto& box = mfi.fabbox();
    const auto& data = fab.array();

    const auto& statusArr = fstatus[mfi].array();
    const auto& tmp = ftmp[mfi].array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int iVar = 0; iVar < f.nComp(); iVar++)
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (bit::is_lev_edge(statusArr(i, j, k))) {
              data(i, j, k, iVar) = tmp(i, j, k, iVar);
            }
          }
  }
}

template <class FAB>
void fill_fine_lev_from_coarse(amrex::FabArray<FAB>& coarse,
                               amrex::FabArray<FAB>& fine, const int iStart,
                               const int nComp, const amrex::IntVect ratio,
                               const amrex::Geometry& cgeom,
                               const amrex::Geometry& fgeom,
                               amrex::Interpolater& mapper) {
  BL_PROFILE("fill_fine_lev_from_coarse");

  amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);

  amrex::FabArray<FAB> ftmp(f.boxArray(), f.DistributionMap(), nComp,
                            fine.nGrow());
  ftmp.setVal(0.0);

  interp_from_coarse_to_fine(c, ftmp, 0, nComp, ratio, cgeom, fgeom, &mapper,
                             f.nGrow());

  for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
    FAB& fab = f[mfi];
    const auto& box = mfi.fabbox();
    const auto& data = fab.array();

    const auto& tmp = ftmp[mfi].array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int iVar = 0; iVar < f.nComp(); iVar++)
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {

            data(i, j, k, iVar) = tmp(i, j, k, iVar);
          }
  }
  fine.FillBoundary();
}

// Sum from coarse level to fine level for nodes at the boundary of two levels.
template <class FAB>
void sum_coarse_to_fine_lev_bny_node(
    amrex::FabArray<FAB>& coarse, amrex::FabArray<FAB>& fine, const int iStart,
    const int nComp, const amrex::IntVect ratio, const amrex::Geometry& cgeom,
    const amrex::Geometry& fgeom, const amrex::iMultiFab& fstatus) {
  BL_PROFILE("sum_coarse_to_fine_lev_bny_node");

  amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);

  amrex::FabArray<FAB> ftmp(f.boxArray(), f.DistributionMap(), nComp, 0);
  ftmp.setVal(0.0);
  amrex::UNodeBilinear<typename FAB::value_type> mapper;
  interp_from_coarse_to_fine(c, ftmp, 0, nComp, ratio, cgeom, fgeom, &mapper);

  for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
    FAB& fab = f[mfi];
    const auto& box = mfi.validbox();
    const auto& data = fab.array();

    const auto& statusArr = fstatus[mfi].array();
    const auto& tmp = ftmp[mfi].array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int iVar = 0; iVar < f.nComp(); iVar++)
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (bit::is_lev_edge(statusArr(i, j, k))) {
              data(i, j, k, iVar) = tmp(i, j, k, iVar);
            }
          }
  }
}

// Sum from coarse level to fine level for domain boundary edge nodes
template <class FAB>
void interp_from_coarse_to_fine_for_domain_edge(
    amrex::FabArray<FAB>& coarse, amrex::FabArray<FAB>& fine, const int iStart,
    const int nComp, const amrex::IntVect ratio, const amrex::Geometry& cgeom,
    const amrex::Geometry& fgeom, const amrex::iMultiFab& fstatus) {
  BL_PROFILE("interp_from_coarse_to_fine_for_domain_edge");

  amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);

  amrex::FabArray<FAB> ftmp(f.boxArray(), f.DistributionMap(), nComp, 0);
  ftmp.setVal(0.0);

  interp_from_coarse_to_fine(c, ftmp, 0, nComp, ratio, cgeom, fgeom,
                             &amrex::node_bilinear_interp);

  for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
    FAB& fab = f[mfi];
    const auto& box = mfi.validbox();
    const auto& data = fab.array();

    const auto& statusArr = fstatus[mfi].array();
    const auto& tmp = ftmp[mfi].array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int iVar = 0; iVar < f.nComp(); iVar++)
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (bit::is_domain_edge(statusArr(i, j, k))) {
              data(i, j, k, iVar) = tmp(i, j, k, iVar);
            }
          }
  }
}

template <class FAB>
void sum_two_lev_interface_node(amrex::FabArray<FAB>& coarse,
                                amrex::FabArray<FAB>& fine, int iStart,
                                const int nComp, const amrex::IntVect ratio,
                                const amrex::Geometry& cgeom,
                                const amrex::Geometry& fgeom,
                                const amrex::iMultiFab& fstatus) {
  BL_PROFILE("sum_two_lev_interface_node");
  sum_fine_to_coarse_lev_bny_node(coarse, fine, iStart, nComp, ratio);
  sum_coarse_to_fine_lev_bny_node(coarse, fine, iStart, nComp, ratio, cgeom,
                                  fgeom, fstatus);
  // amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  // amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);
  // amrex::FabArray<FAB> ftmp(f.boxArray(), f.DistributionMap(), nComp, 0);
  // ftmp.setVal(0.0);
  // amrex::UNodeBilinear<typename FAB::value_type> mapper;
  // interp_from_coarse_to_fine(c, ftmp, 0, nComp, ratio, cgeom, fgeom,
  // &mapper); for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
  //   FAB& fab = f[mfi];
  //   const auto& box = mfi.validbox();
  //   const auto& data = fab.array();
  //   const auto& statusArr = fstatus[mfi].array();
  //   const auto& tmp = ftmp[mfi].array();
  //   const auto lo = amrex::lbound(box);
  //   const auto hi = amrex::ubound(box);
  //   for (int iVar = 0; iVar < f.nComp(); iVar++)
  //     for (int k = lo.z; k <= hi.z; ++k)
  //       for (int j = lo.y; j <= hi.y; ++j)
  //         for (int i = lo.x; i <= hi.x; ++i) {
  //           if (bit::is_lev_edge(statusArr(i, j, k))) {
  //             data(i, j, k, iVar) = data(i, j, k, iVar) + tmp(i, j, k, iVar);
  //           }
  //         }
  // }
  // average_down_nodal(f, c, ratio);
}

// This function is called recursively to combine active patahces into larger
// boxes. The number of boxes should be minimized.
inline void get_boxlist_from_region(amrex::BoxList& bl, GridInfo& gridInfo,
                                    amrex::IntVect imin, amrex::IntVect imax) {

  amrex::IntVect patchSize;
  for (int i = 0; i < nDim; ++i) {
    patchSize[i] = gridInfo.get_patch_size(i);
  }

  amrex::IntVect patchLen = imax - imin + 1;

  for (int i = 0; i < nDim; ++i) {
    if ((patchLen[i] % patchSize[i]) != 0)
      amrex::Abort("Error: the sub-region range is wrong!");
    patchLen[i] /= patchSize[i];
  }

  if (patchLen == 1) {
    if (gridInfo.get_status(AMREX_D_DECL(imin[ix_], imin[iy_], imin[iz_])) ==
        GridInfo::iPicOn_) {
      bl.push_back(amrex::Box(imin, imax, { AMREX_D_DECL(0, 0, 0) }));
    }
    return;
  }

  amrex::IntVect dhalf;
  for (int i = 0; i < nDim; ++i) {
    dhalf[i] = ceil(patchLen[i] / 2.0) * patchSize[i];
  }

  amrex::BoxList blLoc;
  AMREX_D_TERM(for (int i = imin[ix_]; i <= imax[ix_]; i += dhalf[ix_]),
               for (int j = imin[iy_]; j <= imax[iy_]; j += dhalf[iy_]),
               for (int k = imin[iz_]; k <= imax[iz_]; k += dhalf[iz_])) {
    amrex::IntVect iminsub = { AMREX_D_DECL(i, j, k) };
    amrex::IntVect imaxsub = iminsub + dhalf - 1;

    for (int i = 0; i < nDim; ++i) {
      if (imaxsub[i] > imax[i])
        imaxsub[i] = imax[i];
    }
    get_boxlist_from_region(blLoc, gridInfo, iminsub, imaxsub);
  }

  // The number '3' is chosen based on numerical experiments.
  for (int i = 0; i < 3; ++i)
    blLoc.simplify();

  bl.join(blLoc);

  for (int i = 0; i < 3; ++i)
    bl.simplify();
}

#endif

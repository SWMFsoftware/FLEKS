#include "../../srcInterface/multi_ipic3d_domain.h"
#include "Domain.h"

using namespace amrex;

void Domain::find_output_list(const PlotWriter& writerIn,
                              long int& nPointAllProc,
                              PlotWriter::VectorPointList& pointList_II,
                              std::array<double, nDimMax>& xMin_D,
                              std::array<double, nDimMax>& xMax_D) {
  const auto plo = geom.ProbLo();
  const auto plh = geom.ProbHi();

  Real xMinL_D[nDimMax] = { plh[ix_], plh[iy_], plh[iz_] };
  Real xMaxL_D[nDimMax] = { plo[ix_], plo[iy_], plo[iz_] };

  const auto dx = geom.CellSize();

  int iBlock = 0;
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
    FArrayBox& fab = nodeE[mfi];
    const Box& box = mfi.validbox();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    int iMax = hi.x, jMax = hi.y, kMax = hi.z;

    const bool isNode = true;
    if (isNode) {
      // Avoid double counting the shared edges.
      iMax--;
      jMax--;
      kMax--;

      if ((!geom.isPeriodic(ix_)) && box.bigEnd(ix_) == hi.x)
        iMax++;
      if ((!geom.isPeriodic(iy_)) && box.bigEnd(iy_) == hi.y)
        jMax++;
      if ((!geom.isPeriodic(iz_)) && box.bigEnd(iz_) == hi.z)
        kMax++;
    }

    for (int k = lo.z; k <= kMax; ++k) {
      const double zp = k * dx[iz_] + plo[iz_];
      for (int j = lo.y; j <= jMax; ++j) {
        const double yp = j * dx[iy_] + plo[iy_];
        for (int i = lo.x; i <= iMax; ++i) {
          const double xp = i * dx[ix_] + plo[ix_];

          if (writerIn.is_inside_plot_region(i, j, k, xp, yp, zp)) {
            pointList_II.push_back({ (double)i, (double)j, (double)k, xp, yp,
                                     zp, (double)iBlock });
            if (xp < xMinL_D[ix_])
              xMinL_D[ix_] = xp;
            if (yp < xMinL_D[iy_])
              xMinL_D[iy_] = yp;
            if (zp < xMinL_D[iz_])
              xMinL_D[iz_] = zp;

            if (xp > xMaxL_D[ix_])
              xMaxL_D[ix_] = xp;
            if (yp > xMaxL_D[iy_])
              xMaxL_D[iy_] = yp;
            if (zp > xMaxL_D[iz_])
              xMaxL_D[iz_] = zp;
          }
        }
      }
    }
    iBlock++;
  }

  long nPoint = pointList_II.size();
  ParallelDescriptor::ReduceLongSum(nPoint);

  ParallelDescriptor::ReduceRealMin(xMinL_D, nDimMax);
  ParallelDescriptor::ReduceRealMax(xMaxL_D, nDimMax);

  for (int iDim = 0; iDim < nDimMax; ++iDim) {
    xMin_D[iDim] = xMinL_D[iDim];
    xMax_D[iDim] = xMaxL_D[iDim];
  }
}

void Domain::get_field_var(const VectorPointList& pointList_II,
                           const std::vector<std::string>& sVar_I,
                           MDArray<double>& var_II) {
  const int iBlk_ = 6;

  long nPoint = pointList_II.size();
  int nVar = sVar_I.size();

  int iBlockCount = 0;
  long iPoint = 0;
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
    FArrayBox& fab = nodeE[mfi];
    const Box& box = mfi.validbox();

    while (iPoint < nPoint) {
      const int ix = pointList_II[iPoint][ix_];
      const int iy = pointList_II[iPoint][iy_];
      const int iz = pointList_II[iPoint][iz_];
      const int iBlock = pointList_II[iPoint][iBlk_];

      if (iBlock == iBlockCount) {

        for (int iVar = 0; iVar < nVar; ++iVar) {
          var_II(iPoint, iVar) = get_var(sVar_I[iVar], ix, iy, iz, mfi);
        }
        iPoint++;
      } else {
        break;
      }
    }

    iBlockCount++;
  }
}

double Domain::get_var(std::string var, const int ix, const int iy,
                       const int iz, const MFIter& mfi) {

  auto get_is = [var]() {
    std::string::size_type pos;
    std::stringstream ss;
    int is;
    pos = var.find_first_of("0123456789");
    ss << var.substr(pos);
    ss >> is;
    return is;
  };

  double value;
  if (var.substr(0, 1) == "X") {
    const auto plo = geom.ProbLo();
    const auto dx = geom.CellSize();
    value = ix * dx[ix_] + plo[ix_];
  } else if (var.substr(0, 1) == "Y") {
    const auto plo = geom.ProbLo();
    const auto dx = geom.CellSize();
    value = iy * dx[iy_] + plo[iy_];
  } else if (var.substr(0, 1) == "Z") {
    const auto plo = geom.ProbLo();
    const auto dx = geom.CellSize();
    value = iz * dx[iz_] + plo[iz_];
  } else if (var.substr(0, 2) == "Ex") {
    const amrex::Array4<amrex::Real const>& arr = nodeE[mfi].array();
    value = arr(ix, iy, iz, ix_);
  } else if (var.substr(0, 2) == "Ey") {
    const amrex::Array4<amrex::Real const>& arr = nodeE[mfi].array();
    value = arr(ix, iy, iz, iy_);
  } else if (var.substr(0, 2) == "Ez") {
    const amrex::Array4<amrex::Real const>& arr = nodeE[mfi].array();
    value = arr(ix, iy, iz, iz_);
  } else if (var.substr(0, 2) == "Bx") {
    const amrex::Array4<amrex::Real const>& arr = nodeB[mfi].array();
    value = arr(ix, iy, iz, ix_);
  } else if (var.substr(0, 2) == "By") {
    const amrex::Array4<amrex::Real const>& arr = nodeB[mfi].array();
    value = arr(ix, iy, iz, iy_);
  } else if (var.substr(0, 2) == "Bz") {
    const amrex::Array4<amrex::Real const>& arr = nodeB[mfi].array();
    value = arr(ix, iy, iz, iz_);
  } else if (var.substr(0, 4) == "rhoS" || var.substr(0, 3) == "uxS" ||
             var.substr(0, 3) == "uyS" || var.substr(0, 3) == "uzS" ||
             var.substr(0, 4) == "pXXS" || var.substr(0, 4) == "pYYS" ||
             var.substr(0, 4) == "pZZS" || var.substr(0, 4) == "pXYS" ||
             var.substr(0, 4) == "pXZS" || var.substr(0, 4) == "pYZS") {
    int iVar;
    if (var.substr(0, 4) == "rhoS")
      iVar = iRho_;
    if (var.substr(0, 3) == "uxS")
      iVar = iUx_;
    if (var.substr(0, 3) == "uyS")
      iVar = iUy_;
    if (var.substr(0, 3) == "uzS")
      iVar = iUz_;
    if (var.substr(0, 4) == "pXXS")
      iVar = iPxx_;
    if (var.substr(0, 4) == "pYYS")
      iVar = iPyy_;
    if (var.substr(0, 4) == "pZZS")
      iVar = iPzz_;
    if (var.substr(0, 4) == "pXYS")
      iVar = iPxy_;
    if (var.substr(0, 4) == "pXZS")
      iVar = iPxz_;
    if (var.substr(0, 4) == "pYZS")
      iVar = iPyz_;

    const amrex::Array4<amrex::Real const>& arr =
        nodePlasma[get_is()][mfi].array();
    value = arr(ix, iy, iz, iVar);
  } else if (var.substr(0, 2) == "pS") {
    const amrex::Array4<amrex::Real const>& arr =
        nodePlasma[get_is()][mfi].array();
    value = (arr(ix, iy, iz, iPxx_) + arr(ix, iy, iz, iPyy_) +
             arr(ix, iy, iz, iPzz_)) /
            3.0;
  } else {
    value = 0;
  }

  return value;
}

void find_output_list_caller(const PlotWriter& writerIn,
                             long int& nPointAllProc,
                             PlotWriter::VectorPointList& pointList_II,
                             std::array<double, nDimMax>& xMin_D,
                             std::array<double, nDimMax>& xMax_D) {
  MPICs->find_output_list(writerIn, nPointAllProc, pointList_II, xMin_D,
                          xMax_D);
}

void get_field_var_caller(const VectorPointList& pointList_II,
                          const std::vector<std::string>& sVar_I,
                          MDArray<double>& var_II) {
  MPICs->get_field_var(pointList_II, sVar_I, var_II);
}

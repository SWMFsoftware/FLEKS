#ifndef AMREX_LINEAR_SOLVER_H_
#define AMREX_LINEAR_SOLVER_H_

#include <AMReX_GMRES.H>
#include <AMReX_MultiFab.H>

class Pic;

namespace amrex {

/**
 * @brief A wrapper class that satisfies the requirements for amrex::GMRES
 * to use MultiFab as the vector type and Pic::update_E_matvec as the operator.
 */
class AmrexLinearSolver {
public:
  using RT = Real;
  using V = MultiFab;

  AmrexLinearSolver() : m_pic(nullptr), m_iLev(-1) {}
  AmrexLinearSolver(Pic* pic, int iLev) : m_pic(pic), m_iLev(iLev) {}

  void define(Pic* pic, int iLev) {
    m_pic = pic;
    m_iLev = iLev;
  }

  // Required by amrex::GMRES
  void apply(V& Ax, V const& x);

  void assign(V& lhs, V const& rhs) {
    MultiFab::Copy(lhs, rhs, 0, 0, lhs.nComp(), 0);
  }

  RT dotProduct(V const& v1, V const& v2) {
    return MultiFab::Dot(v1, 0, v2, 0, v1.nComp(), 0);
  }

  void increment(V& lhs, V const& rhs, RT a) {
    MultiFab::Saxpy(lhs, a, rhs, 0, 0, lhs.nComp(), 0);
  }

  void linComb(V& lhs, RT a, V const& rhs_a, RT b, V const& rhs_b) {
    MultiFab::LinComb(lhs, a, rhs_a, 0, b, rhs_b, 0, 0, lhs.nComp(), 0);
  }

  V makeVecRHS();
  V makeVecLHS();

  RT norm2(V const& v) { return v.norm2(0, v.nComp()); }

  void precond(V& lhs, V const& rhs) {
    // No preconditioning for now, just copy
    MultiFab::Copy(lhs, rhs, 0, 0, lhs.nComp(), 0);
  }

  void scale(V& v, RT fac) { v.mult(fac); }

  void setToZero(V& v) { v.setVal(0.0); }

private:
  Pic* m_pic;
  int m_iLev;
};

} // namespace amrex

#endif

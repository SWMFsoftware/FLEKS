#include "AmrexLinearSolver.h"
#include "Pic.h"

namespace amrex {

void AmrexLinearSolver::apply(V& Ax, V const& x) {
    // We assume Ax and x have the correct number of components and ghost cells.
    // Pic::update_E_matvec should be updated to take MultiFab.
    m_pic->update_E_matvec(x, Ax, m_iLev, true);
}

MultiFab AmrexLinearSolver::makeVecRHS() {
    // Vector for RHS usually doesn't need ghost cells.
    return MultiFab(m_pic->get_n_grids(m_iLev), m_pic->DistributionMap(m_iLev), 3, 0);
}

MultiFab AmrexLinearSolver::makeVecLHS() {
    // Vector for LHS (solution) needs ghost cells for matrix-vector product.
    return MultiFab(m_pic->get_n_grids(m_iLev), m_pic->DistributionMap(m_iLev), 3, m_pic->get_n_ghost());
}


} // namespace amrex

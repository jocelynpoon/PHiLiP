#ifndef __LINEAR_SOLVER_H__
#define __LINEAR_SOLVER_H__

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include "parameters/all_parameters.h"

namespace PHiLiP {

    // Still need to make a LinearSolver class for our problems
    // Note that right hand side should be const
    //  however, the Trilinos wrapper gives and error when trying to
    //  map it. This is probably because the Trilinos function 
    //  does not take right_hand_side as a const
    template <typename real>
    std::pair<unsigned int, double>
    //LinearSolver<real>::solve_linear (
    solve_linear (
        dealii::TrilinosWrappers::SparseMatrix &system_matrix,
        dealii::Vector<real> &right_hand_side,
        dealii::Vector<real> &solution,
        const Parameters::LinearSolverParam &param);

} // PHiLiP namespace

#endif

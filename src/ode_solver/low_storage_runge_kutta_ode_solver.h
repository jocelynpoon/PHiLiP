#ifndef __LOW_STORAGE_RUNGE_KUTTA_ODESOLVER__
#define __LOW_STORAGE_RUNGE_KUTTA_ODESOLVER__

#include "JFNK_solver/JFNK_solver.h"
#include "dg/dg_base.hpp"
#include "ode_solver_base.h"
#include "runge_kutta_methods/low_storage_rk_tableau_base.h"
#include "relaxation_runge_kutta/empty_RRK_base.h"

namespace PHiLiP {
namespace ODE {

/// Runge-Kutta ODE solver (explicit or implicit) derived from ODESolver.
#if PHILIP_DIM==1
template <int dim, typename real, int n_rk_stages, typename MeshType = dealii::Triangulation<dim>>
#else
template <int dim, typename real, int n_rk_stages, typename MeshType = dealii::parallel::distributed::Triangulation<dim>>
#endif
class LowStorageRungeKuttaODESolver: public ODESolverBase <dim, real, MeshType>
{
public:
    LowStorageRungeKuttaODESolver(std::shared_ptr< DGBase<dim, real, MeshType> > dg_input,
            std::shared_ptr<LowStorageRKTableauBase<dim,real,MeshType>> rk_tableau_input,
            std::shared_ptr<EmptyRRKBase<dim,real,MeshType>> RRK_object_input); ///< Constructor.

    /// Function to evaluate solution update
    double err_time_step(real dt, const bool pseudotime);

    void step_in_time(real dt, const bool pseudotime);

    /// Function to allocate the ODE system
    void allocate_ode_system ();

protected:
    /// Stores Butcher tableau a and b, which specify the RK method
    std::shared_ptr<LowStorageRKTableauBase<dim,real,MeshType>> butcher_tableau;

    /// Stores functions related to relaxation Runge-Kutta (RRK).
    /// Functions are empty by default.
    std::shared_ptr<EmptyRRKBase<dim,real,MeshType>> relaxation_runge_kutta;

    /// Implicit solver for diagonally-implicit RK methods, using Jacobian-free Newton-Krylov 
    /** This is initialized for any RK method, but solution-sized vectors are 
     *  only initialized if there is an implicit solve
     */
    JFNKSolver<dim,real,MeshType> solver;
    
    /// Storage for the derivative at each Runge-Kutta stage
    std::vector<dealii::LinearAlgebra::distributed::Vector<double>> rk_stage;
    
    /// Indicator for zero diagonal elements; used to toggle implicit solve.
    std::vector<bool> butcher_tableau_aii_is_zero;
/*
    std::vector<dealii::LinearAlgebra::distributed::Vector<double>> storage_register_2;
    std::vector<dealii::LinearAlgebra::distributed::Vector<double>> storage_register_1;
    std::vector<dealii::LinearAlgebra::distributed::Vector<double>> storage_register_3;
    std::vector<dealii::LinearAlgebra::distributed::Vector<double>> storage_register_4;
    std::vector<dealii::LinearAlgebra::distributed::Vector<double>> rhs;
*/ 
    
    real w;

    double epsilon[3];
    int num_delta;
    int rk_order;
    double atol;
    double rtol;
    //double beta_controller[3];
};

} // ODE namespace
} // PHiLiP namespace

#endif


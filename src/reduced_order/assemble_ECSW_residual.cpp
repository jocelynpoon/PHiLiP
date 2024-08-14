#include "assemble_ECSW_residual.h"
#include <iostream>

#include "flow_solver/flow_solver.h"
#include "flow_solver/flow_solver_factory.h"

namespace PHiLiP {
namespace HyperReduction {
using Eigen::MatrixXd;

template <int dim, int nstate>
AssembleECSWRes<dim,nstate>::AssembleECSWRes(
    const PHiLiP::Parameters::AllParameters *const parameters_input,
    const dealii::ParameterHandler &parameter_handler_input,
    std::shared_ptr<DGBase<dim,double>> &dg_input, 
    std::shared_ptr<ProperOrthogonalDecomposition::PODBase<dim>> pod, 
    MatrixXd snapshot_parameters_input,
    Parameters::ODESolverParam::ODESolverEnum ode_solver_type)
        : AssembleECSWBase<dim, nstate>(parameters_input, parameter_handler_input, dg_input, pod, snapshot_parameters_input, ode_solver_type)
{
}

template <int dim, int nstate>
void AssembleECSWRes<dim,nstate>::build_problem(){
    std::cout << "Solve for A and b for the NNLS Problem from POD Snapshots"<< std::endl;
    MatrixXd snapshotMatrix = this->pod->getSnapshotMatrix();
    const Epetra_CrsMatrix epetra_pod_basis = this->pod->getPODBasis()->trilinos_matrix();
    Epetra_CrsMatrix epetra_system_matrix = this->dg->system_matrix.trilinos_matrix();

    // Get dimensions of the problem
    int num_snaps_POD = snapshotMatrix.cols(); // Number of snapshots used to build the POD basis
    int n_reduced_dim_POD = epetra_pod_basis.NumGlobalCols(); // Reduced subspace dimension
    int N_FOM_dim = epetra_pod_basis.NumGlobalRows(); // Length of solution vector
    int num_elements_N_e = this->dg->triangulation->n_active_cells(); // Number of elements (equal to N if there is one DOF per cell)

    // Create empty and temporary C and d structs
    Epetra_MpiComm epetra_comm(MPI_COMM_WORLD);
    int training_snaps;
    // Check if all or a subset of the snapshots will be used for training
    if (this->all_parameters->hyper_reduction_param.num_training_snaps != 0) {
        std::cout << "LIMITED NUMBER OF SNAPSHOTS"<< std::endl;
        training_snaps = this->all_parameters->hyper_reduction_param.num_training_snaps;
    }
    else{
        training_snaps = num_snaps_POD;
    }
    Epetra_Map RowMap((n_reduced_dim_POD*training_snaps), 0, epetra_comm); // Number of rows in residual based training matrix = n * (number of training snapshots)
    Epetra_Map ColMap(num_elements_N_e, 0, epetra_comm);

    Epetra_CrsMatrix C(Epetra_DataAccess::Copy, RowMap, num_elements_N_e);
    Epetra_Vector d(RowMap);

    // Loop through the given number of training snapshots to find residuals
    const unsigned int max_dofs_per_cell = this->dg->dof_handler.get_fe_collection().max_dofs_per_cell();
    std::vector<dealii::types::global_dof_index> current_dofs_indices(max_dofs_per_cell); 
    int row_num = 0;
    int snap_num = 0;
    for(auto snap_param : this->snapshot_parameters.rowwise()){
        std::cout << "Snapshot Parameter Values" << std::endl;
        std::cout << snap_param << std::endl;
        dealii::LinearAlgebra::ReadWriteVector<double> snapshot_s;
        snapshot_s.reinit(num_elements_N_e);
        // Extract snapshot from the snapshotMatrix
        for (int snap_row = 0; snap_row < num_elements_N_e; snap_row++){
            snapshot_s(snap_row) = snapshotMatrix(snap_row, snap_num);
        }
        dealii::LinearAlgebra::distributed::Vector<double> reference_solution(this->dg->solution);
        reference_solution.import(snapshot_s, dealii::VectorOperation::values::insert);
        
        // Modifiy parameters for snapshot location and create new flow solver
        Parameters::AllParameters params = this->reinitParams(snap_param);
        std::unique_ptr<FlowSolver::FlowSolver<dim,nstate>> flow_solver = FlowSolver::FlowSolverFactory<dim,nstate>::select_flow_case(&params, this->parameter_handler);
        this->dg = flow_solver->dg;

        // Set solution to snapshot and re-compute the residual/Jacobian
        this->dg->solution = reference_solution;
        const bool compute_dRdW = true;
        this->dg->assemble_residual(compute_dRdW);
        Epetra_Vector epetra_right_hand_side(Epetra_DataAccess::Copy, epetra_system_matrix.RowMap(), this->dg->right_hand_side.begin());

        // Compute test basis
        epetra_system_matrix = this->dg->system_matrix.trilinos_matrix();
        std::shared_ptr<Epetra_CrsMatrix> epetra_test_basis = this->local_generate_test_basis(epetra_system_matrix, epetra_pod_basis);

        // Loop through the elements
        for (const auto &cell : this->dg->dof_handler.active_cell_iterators())
        {
            int cell_num = cell->active_cell_index();
            const int fe_index_curr_cell = cell->active_fe_index();
            const dealii::FESystem<dim,dim> &current_fe_ref = this->dg->fe_collection[fe_index_curr_cell];
            const int n_dofs_curr_cell = current_fe_ref.n_dofs_per_cell();

            current_dofs_indices.resize(n_dofs_curr_cell);
            cell->get_dof_indices(current_dofs_indices);

            // Create L_e matrix for current cell
            Epetra_Map LeRowMap(n_dofs_curr_cell, 0, epetra_comm);
            Epetra_Map LeColMap(N_FOM_dim, 0, epetra_comm);
            Epetra_CrsMatrix L_e(Epetra_DataAccess::Copy, LeRowMap, N_FOM_dim);
            double posOne = 1.0;

            for(int i = 0; i < n_dofs_curr_cell; i++){
                const int col = current_dofs_indices[i];
                L_e.InsertGlobalValues(i, 1, &posOne , &col);
            }
            L_e.FillComplete(LeColMap, LeRowMap);

            // Extract residual contributions of the current cell into global dimension
            Epetra_Vector local_r(LeRowMap);
            Epetra_Vector global_r_e(LeColMap);
            L_e.Multiply(false, epetra_right_hand_side, local_r);
            L_e.Multiply(true, local_r, global_r_e);

            // Find reduced-order representation of contribution
            Epetra_Map cseRowMap(n_reduced_dim_POD, 0, epetra_comm);
            Epetra_Vector c_se(cseRowMap);

            epetra_test_basis->Multiply(true, global_r_e, c_se);
            double *c_se_array = new double[n_reduced_dim_POD];

            c_se.ExtractCopy(c_se_array);
            
            // Sub into entries of C and d
            for (int k = 0; k < n_reduced_dim_POD; ++k){
                int place = row_num+k;
                C.InsertGlobalValues(place, 1, &c_se_array[k], &cell_num);
                d.SumIntoGlobalValues(1, &c_se_array[k], &place);
            }
            
        }
        row_num+=n_reduced_dim_POD;
        snap_num+=1;

        // Check if number of training snapshots has been reached
        if (this->all_parameters->hyper_reduction_param.num_training_snaps != 0) {
            std::cout << "LIMITED NUMBER OF SNAPSHOTS"<< std::endl;
            if (snap_num > (this->all_parameters->hyper_reduction_param.num_training_snaps-1)){
                break;
            }
        }
    }

    C.FillComplete(ColMap, RowMap);

    // Sub temp C and d into class A and b
    this->A->reinit(C);
    this->b.reinit(d.GlobalLength());
    for(int z = 0 ; z < d.GlobalLength() ; z++){
        this->b(z) = d[z];
    }
}

#if PHILIP_DIM==1
    template class AssembleECSWRes<PHILIP_DIM, PHILIP_DIM>;
#endif

#if PHILIP_DIM!=1
    template class AssembleECSWRes<PHILIP_DIM, PHILIP_DIM+2>;
#endif

} // HyperReduction namespace
} // PHiLiP namespace
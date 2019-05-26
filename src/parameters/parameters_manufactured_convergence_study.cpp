#include "parameters_manufactured_convergence_study.h"

namespace PHiLiP {
namespace Parameters {

// Manufactured Solution inputs
ManufacturedConvergenceStudyParam::ManufacturedConvergenceStudyParam () {}

void ManufacturedConvergenceStudyParam::declare_parameters (dealii::ParameterHandler &prm)
{
    prm.enter_subsection("manufactured solution convergence study");
    {
        //prm.declare_entry("output", "quiet",
        //                  dealii::Patterns::Selection("quiet|verbose"),
        //                  "State whether output from solver runs should be printed. "
        //                  "Choices are <quiet|verbose>.");
        prm.declare_entry("grid_type", "hypercube",
                          dealii::Patterns::Selection("hypercube|sinehypercube|read_grid"),
                          "Enum of generated grid. "
                          "If read_grid, must have grids xxxx#.msh, where # is the grid numbering from 0 to number_of_grids-1."
                          "Choices are <hypercube|sinehypercube|read_grid>.");

        prm.declare_entry("input_grids", "xxxx",
                          dealii::Patterns::FileName(dealii::Patterns::FileName::FileType::input),
                          "Prefix of Gmsh grids xxxx#.msh used in the grid convergence if read_grid is chosen as the grid_type. ");

        prm.declare_entry("output_meshes", "true",
                          dealii::Patterns::Bool(),
                          "Writes out meshes used for the simulation."
                          "Output will be Gmsh grids named grid-#.msh");

        prm.declare_entry("random_distortion", "0.0",
                          dealii::Patterns::Double(0.0, 0.5),
                          "Randomly disturb grid."
                          "Displaces node by percentage of longest associated edge.");

        prm.declare_entry("initial_grid_size", "2",
                          dealii::Patterns::Integer(),
                          "Initial grid of size (initial_grid_size)^dim");
        prm.declare_entry("number_of_grids", "4",
                          dealii::Patterns::Integer(),
                          "Number of grids in grid study");
        prm.declare_entry("grid_progression", "1.5",
                          dealii::Patterns::Double(),
                          "Multiplier on grid size. "
                          "nth-grid will be of size (initial_grid^grid_progression)^dim");

        prm.declare_entry("degree_start", "0",
                          dealii::Patterns::Integer(),
                          "Starting degree for convergence study");
        prm.declare_entry("degree_end", "3",
                          dealii::Patterns::Integer(),
                          "Last degree used for convergence study");
    }
    prm.leave_subsection();
}

void ManufacturedConvergenceStudyParam ::parse_parameters (dealii::ParameterHandler &prm)
{
    prm.enter_subsection("manufactured solution convergence study");
    {
        //const std::string output_string = prm.get("output");
        //if (output_string == "verbose") output = verbose;
        //if (output_string == "quiet") output = quiet;
        const std::string grid_string = prm.get("grid_type");
        if (grid_string == "hypercube") grid_type = GridEnum::hypercube;
        if (grid_string == "sinehypercube") grid_type = GridEnum::sinehypercube;
        if (grid_string == "read_grid") {
            grid_type = GridEnum::read_grid;
            input_grids = prm.get("input_grids");
        }

        random_distortion           = prm.get_double("random_distortion");
        output_meshes               = prm.get_bool("output_meshes");

        degree_start                = prm.get_integer("degree_start");
        degree_end                  = prm.get_integer("degree_end");

        initial_grid_size           = prm.get_integer("initial_grid_size");
        number_of_grids             = prm.get_integer("number_of_grids");
        grid_progression            = prm.get_double("grid_progression");
    }
    prm.leave_subsection();
}

} // Parameters namespace
} // PHiLiP namespace

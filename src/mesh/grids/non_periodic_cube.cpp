#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <Sacado.hpp>
#include "non_periodic_cube.h"

namespace PHiLiP::Grids {

template<int dim, typename TriangulationType>
void non_periodic_cube(
    TriangulationType&  grid,
    double              domain_left,
    double              domain_right,
    bool                colorize,
    const int           left_boundary_id) 
{
    dealii::GridGenerator::hyper_cube(grid, domain_left, domain_right, colorize);
    std::cout << left_boundary_id << std::endl;
    // Other tests cases set boundary outside of the flow_solver case so 
    // if left_boundary_id is not set, it is skipped so it can be set elsewhere.
    if (left_boundary_id != 9999) {
        for (auto cell = grid.begin_active(); cell != grid.end(); ++cell) {
            // Set a dummy material ID - Fails without this for some reason
            cell->set_material_id(9002);
            if (cell->face(0)->at_boundary()) cell->face(0)->set_boundary_id(left_boundary_id);
            if (cell->face(1)->at_boundary()) cell->face(1)->set_boundary_id(1001);
        }
    } else if (dim == 2) {
        double xmax = 2.0;
        double xmin = 0.0;
        double ymax = 1.0;
        double ymin = 0.0;

        unsigned int n_subdivisions_x = 100;
        unsigned int n_subdivisions_y = 50;
        
        dealii::Point<dim> p1;
        dealii::Point<dim> p2;
        p1[0] = xmin; p1[1] = ymin;
        p2[0] = xmax; p2[1] = ymax;
        
        std::vector<unsigned int> n_subdivisions(2);

        n_subdivisions[0] = n_subdivisions_x;
        n_subdivisions[1] = n_subdivisions_y;

        dealii::GridGenerator::subdivided_hyper_rectangle(grid, n_subdivisions, p1, p2, true);

        // Set boundary type and design type
        for (typename dealii::parallel::distributed::Triangulation<dim>::active_cell_iterator cell = grid.begin_active(); cell != grid.end(); ++cell) {
            for (unsigned int face = 0; face < dealii::GeometryInfo<2>::faces_per_cell; ++face) {
                if (cell->face(face)->at_boundary()) {
                    unsigned int current_id = cell->face(face)->boundary_id();
                    if (current_id == 0) {
                        cell->face(face)->set_boundary_id(1007); // left bc, post shock condition
                    }
                    else if (current_id == 1) {
                        cell->face(face)->set_boundary_id(1002); // right bc, subsonic outlet
                    }
                    else if (current_id == 2) {
                        cell->face(face)->set_boundary_id(1001); // bottom bc, wall
                    }
                    else if (current_id == 3) {
                        cell->face(face)->set_boundary_id(1001); // top bc, wall
                    }
                }
            }
        }
    }
}

#if PHILIP_DIM==1
template void non_periodic_cube<1, dealii::Triangulation<1>>(
    dealii::Triangulation<1>&   grid,
    double                      domain_left,
    double                      domain_right,
    bool                        colorize,
    const int                   left_boundary_id);
#else
template void non_periodic_cube<2, dealii::parallel::distributed::Triangulation<2>>(
    dealii::parallel::distributed::Triangulation<2>&    grid,
    double                                              domain_left,
    double                                              domain_right,
    bool                                                colorize,
    const int                                           left_boundary_id);
#endif
} // namespace PHiLiP::Grids

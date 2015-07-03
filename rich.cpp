#include <iostream>
#include <fstream>
#include "polyfit.hpp"
#include "source/newtonian/two_dimensional/physical_geometry.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/misc/mesh_generator.hpp"
#include "fermi_table.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "source/newtonian/two_dimensional/simple_cfl.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include "source/misc/simple_io.hpp"
#include "temperature_appendix.hpp"
#include "units.hpp"
#include "bracketed.hpp"
#include "safe_map_read.hpp"
#include "vector_io.hpp"
#include "vector_utils.hpp"
#include "get_composition_data.hpp"
#include "interpolator.hpp"
#include "initial_data.hpp"
#include "create_pressure_reference.hpp"
#include "rectangle_stretch.hpp"
#include "source/tessellation/right_rectangle.hpp"
#include "source/newtonian/test_2d/clip_grid.hpp"
#include "sim_data.hpp"
#include "my_main_loop.hpp"

using namespace std;

int main(void)
{
  Units units;
  SimData sim_data(InitialData("radius_list.txt",
			       "density_list.txt",
			       "temperature_list.txt",
			       "velocity_list.txt"),
		   units);
  hdsim& sim = sim_data.getSim();
  my_main_loop(sim,sim_data.getEOS());
  return 0;
}

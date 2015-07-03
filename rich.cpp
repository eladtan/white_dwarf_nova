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
#include "source/newtonian/test_2d/multiple_diagnostics.hpp"
#include "nuclear_burn.hpp"
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
#include "write_cycle.hpp"
#include "sim_data.hpp"

using namespace std;
using namespace simulation2d;

int main(void)
{
  Units units;
  SimData sim_data(InitialData("radius_list.txt",
			       "density_list.txt",
			       "temperature_list.txt",
			       "velocity_list.txt"),
		   units);
  hdsim& sim = sim_data.getSim();
  write_snapshot_to_hdf5(sim,"initial.h5");
  const double tf = 10;
  SafeTimeTermination term_cond(tf,1e6);
  vector<DiagnosticFunction*> diag_list = VectorInitialiser<DiagnosticFunction*>()
     [new ConsecutiveSnapshots(new ConstantTimeInterval(tf/1000),
			       new Rubric("snapshot_",".h5"),
			       vector<DiagnosticAppendix*>
			       (1,new TemperatureAppendix(sim_data.getEOS())))]
     [new WriteTime("time.txt")]
     [new WriteCycle("cycle.txt")]
    ();
  MultipleDiagnostics diag(diag_list);
  NuclearBurn manip(string("alpha_table"),
		    string("ghost"),
		    sim_data.getEOS());
  main_loop(sim,
	    term_cond,
	    &hdsim::TimeAdvance,
	    &diag,
	    &manip);
  write_snapshot_to_hdf5(sim,"final.h5",
			 vector<DiagnosticAppendix*>
			 (1,new TemperatureAppendix(sim_data.getEOS())));
  return 0;
}

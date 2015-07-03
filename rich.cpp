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
#include "source/newtonian/two_dimensional/source_terms/CenterGravity.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylindrical_complementary.hpp"
#include "source/newtonian/two_dimensional/source_terms/SeveralSources.hpp"
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
#include "generate_atomic_properties.hpp"
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
#include "create_grid.hpp"
#include "inner_bc.hpp"
#include "lazy_cell_updater.hpp"
#include "calc_init_cond.hpp"
#include "circular_section.hpp"
#include "lazy_extensive_updater.hpp"
#include "monopole_self_gravity.hpp"
#include "write_cycle.hpp"

using namespace std;
using namespace simulation2d;

class SimData
{
public:

  SimData(const InitialData& id, const Units& u):
    pg_(Vector2D(0,0), Vector2D(1,0)),
    outer_(Vector2D(-0.5*id.radius_mid.front(),0.9*id.radius_mid.front()),
	   Vector2D(0.5*id.radius_mid.front(),1.2*id.radius_mid.back())),
    //  tess_(cartesian_mesh(100,100,outer_.getBoundary().first,outer_.getBoundary().second),
    tess_(create_grid(outer_.getBoundary(),1e-2,0.9*id.radius_list.front()),
	  outer_),
    eos_("eos_tab.coded",1,1,0,generate_atomic_properties()),
    rs_(),
    point_motion_(),
    gravity_acc_(u.gravitation_constant*1.816490e33*u.gram,
		 0, Vector2D(0,0)),
    gravity_force_(gravity_acc_),
    msg_(linspace(id.radius_list.front(),id.radius_list.back(),100),
	 u.gravitation_constant,
	 pair<double,double>(M_PI*0.46,M_PI*0.54)),
    geom_force_(pg_.getAxis()),
    force_(VectorInitialiser<SourceTerm*>(&gravity_force_)
	   (&geom_force_)(&msg_)()),
    tsf_(0.3),
    fc_(rs_,string("ghost"),id.radius_mid.back()),
    eu_(),
    cu_(),
    sim_(tess_,
	 outer_,
	 pg_,
	 calc_init_cond(tess_,eos_,id,CircularSection(id.radius_mid.front(),
						      id.radius_mid.back(),
						      0.46*M_PI,
						      0.54*M_PI)),
	 eos_,
	 point_motion_,
	 force_,
	 tsf_,
	 fc_,
	 eu_,
	 cu_) {}

  hdsim& getSim(void)
  {
    return sim_;
  }

  const FermiTable& getEOS(void) const
  {
    return eos_;
  }

private:
  const CylindricalSymmetry pg_;
  const SquareBox outer_;
  VoronoiMesh tess_;
  const FermiTable eos_;
  const Hllc rs_;
  Eulerian point_motion_;
  CenterGravity gravity_acc_;
  ConservativeForce gravity_force_;
  MonopoleSelfGravity msg_;
  CylindricalComplementary geom_force_;
  SeveralSources force_;
  const SimpleCFL tsf_;
  //  const SimpleFluxCalculator fc_;
  const InnerBC fc_;
  //  const SimpleExtensiveUpdater eu_;
  const LazyExtensiveUpdater eu_;
  //  const SimpleCellUpdater cu_;
  const LazyCellUpdater cu_;
  hdsim sim_;
};

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

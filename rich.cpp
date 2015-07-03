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

using namespace std;
using namespace simulation2d;

namespace {

  template<typename T> vector<T> polyfit_sc(const pair<vector<T>,vector<T> >& x_y, int deg)
  {
    return polyfit(x_y.first, x_y.second, deg);
  }

  vector<pair<double,double> > calc_mass_radius_list(const Tessellation& tess,
						     const vector<ComputationalCell>& cells,
						     const CacheData& cd)
  {
    vector<pair<double, double> > res;
    for(size_t i=0;i<cells.size();++i){
      if(cells[i].stickers.find("ghost")->second)
	continue;
      const double radius = abs(tess.GetCellCM(static_cast<int>(i)));
      const double mass = cd.volumes[i]*cells[i].density;
      res.push_back(pair<double,double>(radius,mass));
    }
    return res;
  }

  vector<double> calc_mass_in_shells(const vector<pair<double,double> >& mass_radius_list,
				     const vector<double>& sample_points)
  {
    vector<double> res(sample_points.size(),0);
    for(vector<pair<double,double> >::const_iterator it=
	  mass_radius_list.begin();
	it!=mass_radius_list.end();
	++it){
      for(size_t i=0;i<sample_points.size();++i){
	if(sample_points[i]>it->first){
	  res[i] += it->second;
	  //	  break;
	}
      }
    }
    return res;
  }

  /*
  vector<double> cumsum(const vector<double>& v)
  {
    vector<double> res(v.size(),0);
    res[0] = v[0];
    for(size_t i=1;i<v.size();++i)
      res[i] = res[i-1] + v[i];
    return res;
  }
  */
  
  class MonopoleSelfGravity: public SourceTerm
  {
  public:

    MonopoleSelfGravity(const vector<double>& sample_radii,
			double gravitation_constant,
			const pair<double,double>& section_angles):
      sample_radii_(sample_radii),
      gravitation_constant_(gravitation_constant),
      section2shell_(2./(cos(section_angles.first)-cos(section_angles.second))) {}

    vector<Extensive> operator()
    (const Tessellation& tess,
     const PhysicalGeometry& /*pg*/,
     const CacheData& cd,
     const vector<ComputationalCell>& cells,
     const vector<Extensive>& /*fluxes*/,
     const vector<Vector2D>& /*point_velocities*/,
     const double /*t*/) const
    {
      const vector<double> mass_sample =
	calc_mass_in_shells(calc_mass_radius_list(tess,cells,cd),
			    sample_radii_);
      const Interpolator radius_mass_interp(sample_radii_,
					    mass_sample);

      /*
      {
	ofstream f("mass_radius.txt");
	for(size_t i=0;i<mass_sample.size();++i){
	  f << sample_radii_[i] << " " << mass_sample[i]*section2shell_ << endl;
	}
	f.close();
      }
      */

      vector<Extensive> res(static_cast<size_t>(tess.GetPointNo()));
      for(size_t i=0;i<res.size();++i){
	if(cells[i].stickers.find("ghost")->second)
	  continue;
	const Vector2D r = tess.GetCellCM(static_cast<int>(i));
	const double radius = abs(r);
	const double mass = radius_mass_interp(radius)*section2shell_;
	const Vector2D acc = (-1)*gravitation_constant_*r*mass/pow(radius,3);
	const double volume = cd.volumes[i];
	res[i].mass = 0;
	res[i].momentum = volume*cells[i].density*acc;
	res[i].energy = volume*cells[i].density*ScalarProd(acc,cells[i].velocity);
      }
      return res;
    }

  private:
    const vector<double> sample_radii_;
    const double gravitation_constant_;
    const double section2shell_;
  };

  class WriteCycle: public DiagnosticFunction
  {
  public:

    WriteCycle(const string& fname):
      fname_(fname) {}

    void operator()(const hdsim& sim)
    {
      write_number(sim.getCycle(),fname_);
    }

  private:
    const string fname_;
  };
}

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

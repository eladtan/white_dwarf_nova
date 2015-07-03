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
#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
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

using namespace std;
using namespace simulation2d;

namespace {

  map<string,vector<double> > get_composition_data(void)
  {
    map<string,vector<double> > res;
    const map<string,pair<double,double> > atomic_properties = generate_atomic_properties();
    for(map<string,pair<double,double> >::const_iterator it =
	  atomic_properties.begin();
	it != atomic_properties.end(); ++it)
      res[it->first] = decapitate(load_txt(string("tracer_")+it->first+".txt"));
    return res;
  }

  class Interpolator
  {
  public:

    Interpolator(const vector<double>& x_list,
		 const vector<double>& y_list):
      x_list_(x_list), y_list_(y_list)
    {
      assert(is_strictly_increasing(x_list_));
      assert(x_list_.size()==y_list_.size());
    }

    double operator()(double x) const
    {
      const double x_list_front = x_list_.front();
      const double x_list_back = x_list_.back();
      assert(x>x_list_front);
      assert(x<x_list_back);       
      //      assert(x>x_list_.front());
      //      assert(x<x_list_.back());
      for(size_t i=1;i<x_list_.size();++i){
	if(x_list_.at(i)>x)
	  return y_list_.at(i-1) + (y_list_.at(i)-y_list_.at(i-1))*
	    (x-x_list_.at(i-1))/(x_list_.at(i)-x_list_.at(i-1));
      }
      throw "point outside bound";
    }

  private:
    const vector<double> x_list_;
    const vector<double> y_list_;
  };

  vector<double> mid_array(const vector<double>& v)
  {
    assert(v.size()>0);
    vector<double> res(v.size()-1);
    for(size_t i=0;i<res.size();++i)
      res.at(i) = 0.5*(v.at(i)+v.at(i+1));
    return res;
  }

  class InitialData
  {
  public:

    const vector<double> radius_list;
    const vector<double> radius_mid;
    const vector<double> density_list;
    const vector<double> temperature_list;
    const vector<double> velocity_list;
    const map<string,vector<double> > tracers_list;

    InitialData(const string& radius_file,
		const string& density_file,
		const string& temperature_file,
		const string& velocity_file):
      radius_list(load_txt(radius_file)),
      radius_mid(mid_array(radius_list)),
      density_list(decapitate(load_txt(density_file))),
      temperature_list(decapitate(load_txt(temperature_file))),
      velocity_list(load_txt(velocity_file)),
      tracers_list(get_composition_data()) {}
  };

  vector<double> create_pressure_reference(const FermiTable& eos,
					   const InitialData& id)
  {
    vector<double> res(id.radius_list.size());
    for(size_t i=0;i<id.density_list.size();++i){
      map<string,double> tracer;
      for(map<string,vector<double> >::const_iterator it=
	    id.tracers_list.begin();
	  it!=id.tracers_list.end();
	  ++it)
	tracer[it->first] = (it->second).at(i);
      res[i] = eos.dt2p(id.density_list.at(i),
			id.temperature_list.at(i),
			tracer);
    }
    return res;
  }

  pair<Vector2D, Vector2D> stretch(const pair<Vector2D,Vector2D>& boundaries,
				   double ratio)
  {
    const Vector2D centre = 0.5*(boundaries.first + boundaries.second);
    return pair<Vector2D,Vector2D>(centre+ratio*(boundaries.first-centre),
				   centre+ratio*(boundaries.second-centre));
  }

  bool is_in(const Vector2D& point,
	     const pair<Vector2D,Vector2D>& boundaries)
  {
    return ((boundaries.first.x<point.x) &&
	    (boundaries.second.x>point.x) &&
	    (boundaries.first.y<point.y) &&
	    (boundaries.second.y>point.y));
      
  }

  vector<Vector2D> rectangular_clip(const vector<Vector2D>& point_list,
				    const pair<Vector2D,Vector2D>& boundaries)
  {
    vector<Vector2D> res;
    for(size_t i=0;i<point_list.size();++i){
      if(is_in(point_list.at(i),boundaries))
	res.push_back(point_list.at(i));
    }
    return res;
  }

  vector<Vector2D> create_grid(const pair<Vector2D,Vector2D>& boundaries,
			       const double dq,
			       const double r_min)
  {
    vector<Vector2D> res;
    //   const double r_min = 0.1*boundaries.first.y;
    const double r_max = boundaries.second.y;
    for(double r=r_min;
	r<r_max;
	r*=(1+dq)){
      for(double q=0;q<2*M_PI;q+=dq)
	res.push_back(Vector2D(r*cos(q),r*sin(q)));
    }
    res = rectangular_clip(res,stretch(boundaries,0.99));
    ofstream f("mesh_points.txt");
    for(size_t i=0;i<res.size();++i)
      f << res.at(i).x << " " << res.at(i).y << endl;
    f.close();
    return res;
  }

  double calc_tracer_flux(const Edge& edge,
			  const Tessellation& tess,
			  const vector<ComputationalCell>& cells,
			  const string& name,
			  const Conserved& hf)
  {
    if(hf.Mass>0 && 
       edge.neighbors.first>0 && 
       edge.neighbors.first<tess.GetPointNo()){
      assert(cells.at(static_cast<size_t>(edge.neighbors.first)).tracers.count(name)==1);
      return hf.Mass*
	cells.at(static_cast<size_t>(edge.neighbors.first)).tracers.find(name)->second;
    }
    if(hf.Mass<0 && 
       edge.neighbors.second>0 && 
       edge.neighbors.second<tess.GetPointNo()){
      assert(cells.at(static_cast<size_t>(edge.neighbors.second)).tracers.count(name)==1);
      return hf.Mass*
	cells.at(static_cast<size_t>(edge.neighbors.second)).tracers.find(name)->second;
    }
    return 0;    
  }

  class InnerBC: public FluxCalculator
  {
  public:

    InnerBC(const RiemannSolver& rs,
	    const string& ghost,
	    const double radius):
      rs_(rs), ghost_(ghost), radius_(radius) {}

    vector<Extensive> operator()
    (const Tessellation& tess,
     const vector<Vector2D>& point_velocities,
     const vector<ComputationalCell>& cells,
     const EquationOfState& eos,
     const double /*time*/,
     const double /*dt*/) const
    {
      vector<Extensive> res(tess.getAllEdges().size());
      for(size_t i=0;i<tess.getAllEdges().size();++i){
	const Conserved hydro_flux =
	  calcHydroFlux(tess,point_velocities,
			cells, eos, i);
	res.at(i).mass = hydro_flux.Mass;
	res.at(i).momentum = hydro_flux.Momentum;
	res.at(i).energy = hydro_flux.Energy;
	for(map<string,double>::const_iterator it =
	      cells.front().tracers.begin();
	    it!=cells.front().tracers.end();
	    ++it)
	  res.at(i).tracers[it->first] =
	    calc_tracer_flux(tess.getAllEdges().at(i),
			     tess,cells,it->first,hydro_flux);
      }
      return res;
    }

  private:
    const RiemannSolver& rs_;
    const string ghost_;
    const double radius_;

    const Conserved calcHydroFlux
    (const Tessellation& tess,
     const vector<Vector2D>& point_velocities,
     const vector<ComputationalCell>& cells,
     const EquationOfState& eos,
     const size_t i) const
    {
      const Edge& edge = tess.GetEdge(static_cast<int>(i));
      const pair<bool,bool> flags
	(edge.neighbors.first>=0 && edge.neighbors.first<tess.GetPointNo(),
	 edge.neighbors.second>=0 && edge.neighbors.second<tess.GetPointNo());
      assert(flags.first || flags.second);
      if(!flags.first){
	const size_t right_index = static_cast<size_t>(edge.neighbors.second);
	const ComputationalCell& right_cell = cells.at(right_index);
	if(right_cell.stickers.find(ghost_)->second)
	  return Conserved();
	const Vector2D p = Parallel(edge);
	const Primitive right = convert_to_primitive(right_cell, eos);
	const Primitive left = reflect(right,p);
	const Vector2D n = remove_parallel_component
	  (tess.GetMeshPoint(edge.neighbors.second)-
	   edge.vertices.first,p);
	return rotate_solve_rotate_back
	  (rs_, left, right, 0, n, p);
      }
      if(!flags.second){
	const size_t left_index =
	  static_cast<size_t>(edge.neighbors.first);
	const ComputationalCell& left_cell = cells.at(left_index);
	if(left_cell.stickers.find(ghost_)->second)
	  return Conserved();
	const Primitive left = convert_to_primitive(left_cell, eos);
	const Vector2D p = Parallel(edge);
	const Primitive right = reflect(left, p);
	const Vector2D n = remove_parallel_component
	  (edge.vertices.second -
	   tess.GetMeshPoint(edge.neighbors.first),p);
	return rotate_solve_rotate_back
	  (rs_, left, right, 0, n, p);
      }
      const size_t left_index =
	static_cast<size_t>(edge.neighbors.first);
      const size_t right_index =
	static_cast<size_t>(edge.neighbors.second);
      const ComputationalCell& left_cell = cells.at(left_index);
      const ComputationalCell& right_cell = cells.at(right_index);
      if(left_cell.stickers.find(ghost_)->second &&
	 right_cell.stickers.find(ghost_)->second)
	return Conserved();
      const Vector2D p = Parallel(edge);
      const Vector2D n =
	tess.GetMeshPoint(edge.neighbors.second) -
	tess.GetMeshPoint(edge.neighbors.first);
      const double velocity = Projection
	(tess.CalcFaceVelocity
	 (point_velocities.at(left_index),
	  point_velocities.at(right_index),
	  tess.GetCellCM(edge.neighbors.first),
	  tess.GetCellCM(edge.neighbors.second),
	  calc_centroid(edge)),n);
      if(left_cell.stickers.find(ghost_)->second){
	const Primitive right =
	  convert_to_primitive(right_cell, eos);
	const Primitive left =
	  (abs(tess.GetMeshPoint(static_cast<int>(left_index)))>radius_) ?
	  right : reflect(right,p);
	return rotate_solve_rotate_back
	  (rs_,left,right,velocity,n,p);
      }
      if(right_cell.stickers.find(ghost_)->second){
	const Primitive left =
	  convert_to_primitive(left_cell, eos);
	const Primitive right =
	  (abs(tess.GetMeshPoint(static_cast<int>(right_index)))>radius_) ?
	  left : reflect(left,p);
	return rotate_solve_rotate_back
	  (rs_,left,right,velocity,n,p);
      }
      const Primitive left =
	convert_to_primitive(left_cell, eos);
      const Primitive right =
	convert_to_primitive(right_cell, eos);
      return rotate_solve_rotate_back
	(rs_,left,right,velocity,n,p);
    }
  };

  class LazyCellUpdater: public CellUpdater
  {
  public:

    LazyCellUpdater(void) {}

    vector<ComputationalCell> operator()
    (const Tessellation& /*tess*/,
     const PhysicalGeometry& /*pg*/,
     const EquationOfState& eos,
     const vector<Extensive>& extensives,
     const vector<ComputationalCell>& old,
     const CacheData& cd) const
    {
      vector<ComputationalCell> res = old;
      for(size_t i=0;i<extensives.size();++i){
	if(old.at(i).stickers.find("ghost")->second)
	  continue;
	const double volume = cd.volumes[i];
	res.at(i).density = extensives.at(i).mass/volume;
	res.at(i).velocity = extensives.at(i).momentum/extensives.at(i).mass;
	const double total_energy = extensives.at(i).energy/extensives.at(i).mass;
	const double kinetic_energy = 0.5*ScalarProd(res.at(i).velocity, res.at(i).velocity);
	const double thermal_energy = total_energy - kinetic_energy;
	for(map<string,double>::const_iterator it =
	      extensives.at(i).tracers.begin();
	    it!=extensives.at(i).tracers.end();
	    ++it)
	  res.at(i).tracers[it->first] = it->second/extensives.at(i).mass;
	res.at(i).pressure = eos.de2p(res.at(i).density, thermal_energy, res.at(i).tracers);
      }
      return res;
    }
  };

  vector<ComputationalCell> calc_init_cond(const Tessellation& tess,
					   const FermiTable& eos,
					   const InitialData& id,
					   const Shape2D& cd)
  {
    save_txt("pressure_reference.txt",create_pressure_reference(eos,id));
    vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
    const Interpolator density_interpolator(id.radius_mid,
					    id.density_list);
    const Interpolator temperature_interpolator(id.radius_mid,
						id.temperature_list);
    const Interpolator velocity_interpolator(id.radius_list,
					     id.velocity_list);
    map<string,Interpolator*> tracer_intepolators;
    for(map<string,vector<double> >::const_iterator it=
	  id.tracers_list.begin();
	it!=id.tracers_list.end(); ++it)
      tracer_intepolators[it->first] = new Interpolator(id.radius_mid,
							it->second);
    for(size_t i=0;i<res.size();++i){
      res.at(i).density = id.density_list.back();
      res.at(i).velocity = Vector2D(0,0);
      res.at(i).stickers["ghost"] = true;
      for(map<string,Interpolator*>::const_iterator it=
	    tracer_intepolators.begin();
	  it!=tracer_intepolators.end();
	  ++it)
	res.at(i).tracers[it->first] = 0;
      res.at(i).tracers["He4"] = 1;
      res.at(i).pressure = eos.dt2p(res.at(i).density,
				    id.temperature_list.back(),
				    res.at(i).tracers);
      const Vector2D r = tess.GetCellCM(static_cast<int>(i));
      const double radius = abs(r);
      /*
	if(radius<id.radius_list.front() || radius>id.radius_list.back())
	continue;
      */
      if(!cd(r))
	continue;
      res.at(i).stickers["ghost"] = false;
      const double density = density_interpolator(radius);
      const double temperature = temperature_interpolator(radius);
      const double velocity = velocity_interpolator(radius);
      for(map<string,Interpolator*>::const_iterator it=
	    tracer_intepolators.begin();
	  it!=tracer_intepolators.end();
	  ++it)
	res.at(i).tracers[it->first] = (*(it->second))(radius);
      const double pressure = eos.dt2p(density, temperature, res.at(i).tracers);
      res.at(i).density = density;
      res.at(i).pressure = pressure;
      res.at(i).velocity = r*velocity/radius;
    }
    for(map<string,Interpolator*>::iterator it=
	  tracer_intepolators.begin();
	it!=tracer_intepolators.end();
	++it)
      delete it->second;
    return res;
  }

  class CircularSection: public Shape2D
  {
  public:

    CircularSection(const double radius_in,
		    const double radius_out,
		    const double angle_left,
		    const double angle_right):
      radius_in_(radius_in),
      radius_out_(radius_out),
      angle_left_(angle_left),
      angle_right_(angle_right) {}

    bool operator()(const Vector2D& r) const
    {
      const double radius = abs(r);
      const double angle = atan2(r.y,r.x);
      return (radius_in_<radius &&
	      radius_out_>radius &&
	      angle_left_<angle &&
	      angle_right_>angle);
    }

  private:
    const double radius_in_;
    const double radius_out_;
    const double angle_left_;
    const double angle_right_;
  };

  class LazyExtensiveUpdater: public ExtensiveUpdater
  {
  public:

    LazyExtensiveUpdater(void) {}

    void operator()
    (const vector<Extensive>& fluxes,
     const PhysicalGeometry& /*pg*/,
     const Tessellation& tess,
     const double dt,
     const CacheData& cd,
     const vector<ComputationalCell>& cells,
     vector<Extensive>& extensive) const
    {
      const vector<Edge>& edge_list = tess.getAllEdges();
      for(size_t i=0;i<edge_list.size();++i){
	const Edge& edge = edge_list.at(i);
	const Extensive delta = dt*cd.areas[i]*fluxes.at(i);
	if(bracketed(0,edge.neighbors.first,tess.GetPointNo()) &&
	   !safe_map_read
	   (cells.at(static_cast<size_t>(edge.neighbors.first)).stickers,string("ghost")))
	  extensive.at(static_cast<size_t>(edge.neighbors.first)) -= delta;
	if(bracketed(0,edge.neighbors.second,tess.GetPointNo()) &&
	   !safe_map_read
	   (cells.at(static_cast<size_t>(edge.neighbors.second)).stickers,string("ghost")))
	  extensive.at(static_cast<size_t>(edge.neighbors.second)) += delta;
      }
    }
  };

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

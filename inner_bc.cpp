#include "inner_bc.hpp"
#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"
#include "safe_retrieve.hpp"

namespace {
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
}

InnerBC::InnerBC(const RiemannSolver& rs,
		 const string& ghost,
		 const double radius,
		 const double acceleration):
  rs_(rs),
  ghost_(ghost),
  radius_(radius),
  a_(acceleration)
{
  assert(radius_>0);
  assert(a_>0);
}

vector<Extensive> InnerBC::operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& /*extensives*/,
   const EquationOfState& eos,
   const double /*time*/,
   const double dt) const
{
  vector<Extensive> res(tess.getAllEdges().size());
  for(size_t i=0;i<tess.getAllEdges().size();++i){
    const Conserved hydro_flux =
      calcHydroFlux(tess,point_velocities,
		    cells, eos, i, dt*a_);
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

namespace {

  Primitive boost(const Primitive& origin,
		  const Vector2D& v)
  {
    Primitive res = origin;
    res.Velocity += v;
    return res;
  }
  
  Conserved support_riemann(const RiemannSolver& rs,
			    const Vector2D& rmp,
			    const Edge& edge,
			    const Primitive& cell,
			    const Vector2D& support,
			    bool left_real)
  {
    const Vector2D& p = Parallel(edge);
    Conserved res = left_real ?
      rotate_solve_rotate_back
      (rs,
       cell,
       boost(reflect(cell,p),support),
       0,
       remove_parallel_component
       (edge.vertices.second-rmp,p),
       p) :
      rotate_solve_rotate_back
      (rs,
       boost(reflect(cell,p),support),
       cell,
       0,
       remove_parallel_component
       (rmp-edge.vertices.second,p),
       p);
    res.Mass = 0;
    res.Energy = 0;
    return res;
  }

  Conserved regular_riemann
  (const RiemannSolver& rs,
   const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const Edge& edge)
  {
    const size_t left_index =
      static_cast<size_t>(edge.neighbors.first);
    const size_t right_index =
      static_cast<size_t>(edge.neighbors.second);
    const Primitive left =
      convert_to_primitive(cells[left_index], eos);
    const Primitive right =
      convert_to_primitive(cells[right_index], eos);
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
    return rotate_solve_rotate_back
      (rs,left,right,velocity,n,p);
  }
}

const Conserved InnerBC::calcHydroFlux
(const Tessellation& tess,
 const vector<Vector2D>& point_velocities,
 const vector<ComputationalCell>& cells,
 const EquationOfState& eos,
 const size_t i,
 double support) const
{
  const Edge& edge = tess.GetEdge(static_cast<int>(i));
  const pair<bool,bool> flags
    (edge.neighbors.first>=0 && edge.neighbors.first<tess.GetPointNo(),
     edge.neighbors.second>=0 && edge.neighbors.second<tess.GetPointNo());
  assert(flags.first || flags.second);
  if(!flags.first)
    return support_riemann
      (rs_,
       tess.GetMeshPoint(edge.neighbors.second),
       edge,
       convert_to_primitive
       (cells.at
	(static_cast<size_t>(edge.neighbors.second)),
	eos),
       Vector2D(0,0),
       false);
  if(!flags.second)
    return support_riemann
      (rs_,
       tess.GetMeshPoint(edge.neighbors.first),
       edge,
       convert_to_primitive
       (cells.at
	(static_cast<size_t>(edge.neighbors.first)),
	eos),
       Vector2D(0,0),
       true);
  const size_t left_index =
    static_cast<size_t>(edge.neighbors.first);
  const size_t right_index =
    static_cast<size_t>(edge.neighbors.second);
  const ComputationalCell& left_cell = cells.at(left_index);
  const ComputationalCell& right_cell = cells.at(right_index);
  if(safe_retrieve(left_cell.stickers,ghost_)){
    if(safe_retrieve(right_cell.stickers,ghost_))
      return Conserved();
    return support_riemann
      (rs_,
       tess.GetMeshPoint(edge.neighbors.second),
       edge,
       convert_to_primitive(right_cell,eos),
       Vector2D(0,support),
       false);
  }
  if(safe_retrieve(right_cell.stickers,ghost_))
    return support_riemann
      (rs_,
       tess.GetMeshPoint(edge.neighbors.first),
       edge,
       convert_to_primitive(left_cell,eos),
       Vector2D(0,support),
       true);
  return regular_riemann
    (rs_,
     tess,
     point_velocities,
     cells,
     eos,
     edge);
}

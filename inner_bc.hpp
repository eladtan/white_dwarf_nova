#ifndef INNER_BC_HPP
#define INNER_BC_HPP 1

#include "source/newtonian/two_dimensional/flux_calculator_2d.hpp"
#include "source/newtonian/common/riemann_solver.hpp"
#include "core_atmosphere_gravity.hpp"
#include "source/newtonian/two_dimensional/HydroBoundaryConditions.hpp"

class InnerBC: public HydroBoundaryConditions
{
public:

  InnerBC(const RiemannSolver& rs,
	  const string& ghost,
	  const EquationOfState& eos_);

  vector< pair< size_t, Extensive > > operator() 
  (const Tessellation& tess, 
   const vector<ComputationalCell>& cells) const;

private:
  const RiemannSolver& rs_;
  const string ghost_;
  const EquationOfState& eos_;

  const Conserved calcHydroFlux
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const size_t i,
   const CoreAtmosphereGravity::AccelerationCalculator& ac) const;
};

#endif // INNER_BC_HPP

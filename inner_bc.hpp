#ifndef INNER_BC_HPP
#define INNER_BC_HPP 1

#include "source/newtonian/two_dimensional/flux_calculator_2d.hpp"
#include "source/newtonian/common/riemann_solver.hpp"
#include "core_atmosphere_gravity.hpp"

class InnerBC: public FluxCalculator
{
public:

  InnerBC(const RiemannSolver& rs,
	  const string& ghost,
	  const CoreAtmosphereGravity& cag);

  vector<Extensive> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& extensives,
   const CacheData& cd,
   const EquationOfState& eos,
   const double /*time*/,
   const double /*dt*/) const;

private:
  const RiemannSolver& rs_;
  const string ghost_;
  const CoreAtmosphereGravity& cag_;

  const Conserved calcHydroFlux
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const size_t i,
   const CoreAtmosphereGravity::AccelerationCalculator& ac) const;
};

#endif // INNER_BC_HPP

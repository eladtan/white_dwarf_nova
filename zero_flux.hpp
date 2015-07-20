#ifndef ZERO_FLUX_HPP
#define ZERO_FLUX_HPP 1

#include "source/newtonian/two_dimensional/flux_calculator_2d.hpp"

class ZeroFlux: public FluxCalculator
{
public:
  ZeroFlux(void);

  vector<Extensive> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& extensives,
   const EquationOfState& eos,
   const double /*time*/,
   const double /*dt*/) const;
};

#endif // ZERO_FLUX_HPP

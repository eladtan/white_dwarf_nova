#include "zero_flux.hpp"

ZeroFlux::ZeroFlux(void) {}

vector<Extensive> ZeroFlux::operator()
  (const Tessellation& tess,
   const vector<Vector2D>& /*point_velocities*/,
   const vector<ComputationalCell>& /*cells*/,
   const vector<Extensive>& /*extensives*/,
   const EquationOfState& /*eos*/,
   const double /*time*/,
   const double /*dt*/) const
{
  return vector<Extensive>(tess.getAllEdges().size());
}

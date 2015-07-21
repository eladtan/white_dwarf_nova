#include "no_cell_update.hpp"

NoCellUpdate::NoCellUpdate(void) {}

vector<ComputationalCell> NoCellUpdate::operator()
(const Tessellation& /*tess*/,
 const PhysicalGeometry& /*pg*/,
 const EquationOfState& /*eos*/,
 const vector<Extensive>& /*extensives*/,
 const vector<ComputationalCell>& old,
 const CacheData& /*cd*/) const
{
  return old;
}

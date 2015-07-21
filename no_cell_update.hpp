#ifndef NO_CELL_UPDATE_HPP
#define NO_CELL_UPDATE_HPP 1

#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"

class NoCellUpdate: public CellUpdater
{
public:
  NoCellUpdate(void);

  vector<ComputationalCell> operator()
  (const Tessellation& /*tess*/,
   const PhysicalGeometry& /*pg*/,
   const EquationOfState& eos,
   const vector<Extensive>& extensives,
   const vector<ComputationalCell>& old,
   const CacheData& cd) const;
};

#endif // NO_CELL_UPDATE_HPP

#ifndef REFLECTIVE_GHOST_THROUGHOUT_HPP
#define REFLECTIVE_GHOST_THROUGHOUT_HPP 1

#include "source/newtonian/two_dimensional/GhostPointGenerator.hpp"

using std::pair;

class ReflectiveGhostThroughout: public GhostPointGenerator
{
public:

  boost::container::flat_map<size_t, ComputationalCell>
  operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells) const;

  pair<ComputationalCell,ComputationalCell>
  GetGhostGradient
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   const vector<pair<ComputationalCell,ComputationalCell> >& gradients,
   size_t ghost_index) const;
};

#endif // REFLECTIVE_GHOST_THROUGHOUT_HPP

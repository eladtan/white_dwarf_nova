#include "reflective_ghost_throughout.hpp"

namespace {
  ComputationalCell reverse_normal_velocity
  (const ComputationalCell& c,
   const Edge& e)
  {
    const Vector2D p =
      normalize(e.vertices.second - e.vertices.first);
    ComputationalCell res = c;
    res.velocity = 2*ScalarProd(res.velocity,p)*p - res.velocity;
    return res;
  }
}

boost::container::flat_map<size_t, ComputationalCell>
ReflectiveGhostThroughout::operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells) const
{
  boost::container::flat_map<size_t, ComputationalCell> res;
  const vector<pair<size_t, size_t> > ghosts =
    GetOuterEdgesIndeces(tess);
  for(size_t i=0;i<ghosts.size();++i){
    const Edge& edge = tess.GetEdge
      (static_cast<int>(ghosts[i].first));
    const ComputationalCell c =
      reverse_normal_velocity
      (cells.at
       (static_cast<size_t>
	(ghosts[i].second==2 ?
	 edge.neighbors.first :
	 edge.neighbors.second)),
       edge);
    res.insert
      (pair<size_t,ComputationalCell>
       (static_cast<size_t>
	(ghosts[i].second == 1 ?
	 edge.neighbors.first :
	 edge.neighbors.second),
	c));
  }
  return res;
}

pair<ComputationalCell,ComputationalCell>
ReflectiveGhostThroughout::GetGhostGradient
(const Tessellation& /*tess*/,
 const vector<ComputationalCell>& cells,
 const vector<pair<ComputationalCell,ComputationalCell> >& /*gradients*/,
 size_t /*ghost_index*/) const
{
  return pair<ComputationalCell,ComputationalCell>
    (0*cells.front(),
     0*cells.front());
}

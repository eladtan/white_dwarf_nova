#include "create_pressure_reference.hpp"

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

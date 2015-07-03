#include "get_composition_data.hpp"
#include "generate_atomic_properties.hpp"
#include "vector_io.hpp"
#include "vector_utils.hpp"

map<string,vector<double> > get_composition_data(void)
{
  map<string,vector<double> > res;
  const map<string,pair<double,double> > atomic_properties = generate_atomic_properties();
  for(map<string,pair<double,double> >::const_iterator it =
	atomic_properties.begin();
      it != atomic_properties.end(); ++it)
    res[it->first] = decapitate(load_txt(string("tracer_")+it->first+".txt"));
  return res;
}

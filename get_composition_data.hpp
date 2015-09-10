#ifndef GET_COMPOSITION_DATA_HPP
#define GET_COMPOSITION_DATA_HPP 1

#include <map>
#include <string>
#include <vector>
#include "boost/container/flat_map.hpp"

using std::string;
using std::vector;
using std::pair;

boost::container::flat_map<string,vector<double> > 
get_composition_data(void);

#endif // GET_COMPOSITION_DATA_HPP

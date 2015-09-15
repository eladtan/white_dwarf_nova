#ifndef GENERATE_ATOMIC_PROPERTIES_HPP
#define GENERATE_ATOMIC_PROPERTIES_HPP 1

#include <map>
#include <string>
#include "boost/container/flat_map.hpp"

using std::string;
using std::pair;

boost::container::flat_map<string,pair<double,double> > generate_atomic_properties(void);

#endif // GENERATE_ATOMIC_PROPERTIES_HPP

#ifndef INITIAL_DATA_HPP
#define INITIAL_DATA_HPP 1

#include <vector>
#include <map>
#include <string>

using std::vector;
using std::map;
using std::string;

class InitialData
{
public:

  const vector<double> radius_list;
  const vector<double> radius_mid;
  const vector<double> density_list;
  const vector<double> temperature_list;
  const vector<double> velocity_list;
  const map<string,vector<double> > tracers_list;

  InitialData(const string& radius_file,
	      const string& density_file,
	      const string& temperature_file,
	      const string& velocity_file);
};

#endif // INITIAL_DATA_HPP

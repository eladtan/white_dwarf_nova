#ifndef NUCLEAR_BURN_HPP
#define NUCLEAR_BURN_HPP 1

#include <map>
#include <string>
#include "source/newtonian/test_2d/main_loop_2d.hpp"

using std::map;
using std::string;
using std::pair;

class NuclearBurn: public Manipulate
{
public:
  NuclearBurn(const map<string,pair<double,double> >& atomic_propeties,
	      const string& ignore_label);

  void operator()(hdsim& sim);
  
private:

  const map<string,pair<double,double> > atomic_properties_;
  mutable double t_prev_;
  const string ignore_label_;
  const FermiTable& eos_;
};

#endif // NUCLEAR_BURN_HPP

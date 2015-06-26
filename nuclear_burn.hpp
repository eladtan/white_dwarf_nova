#ifndef NUCLEAR_BURN_HPP
#define NUCLEAR_BURN_HPP 1

#include <map>
#include <string>
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "fermi_table.hpp"

using std::map;
using std::string;
using std::pair;

class NuclearBurn: public Manipulate
{
public:
  NuclearBurn(const string& ignore_label,
	      const FermiTable& eos);

  void operator()(hdsim& sim);
  
private:

  mutable double t_prev_;
  const string ignore_label_;
  const FermiTable& eos_;
};

#endif // NUCLEAR_BURN_HPP

#ifndef FILTERED_CONSERVED_HPP
#define FILTERED_CONSERVED_HPP 1

#include <string>
#include "source/newtonian/test_2d/main_loop_2d.hpp"

using std::string;

class FilteredConserved: public DiagnosticFunction
{
public:

  FilteredConserved(const string& fname);

  void operator()(const hdsim& sim);

  ~FilteredConserved(void);

private:
  mutable vector<pair<double,Extensive> > data_;
  const string fname_;
};

#endif // FILTERED_CONSERVED_HPP

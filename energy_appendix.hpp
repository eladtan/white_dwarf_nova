#ifndef ENERGY_APPENDIX_HPP
#define ENERGY_APPENDIX_HPP 1

#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "fermi_table.hpp"

class EnergyAppendix: public DiagnosticAppendix
{
public:

  EnergyAppendix(const FermiTable& eos);

  string getName(void) const;

  vector<double> operator()(const hdsim& sim) const;

private:
  const FermiTable& eos_;
};

#endif // ENERGY_APPENDIX_HPP

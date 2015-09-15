#ifndef CREATE_PRESSURE_REFERENCE_HPP
#define CREATE_PRESSURE_REFERENCE_HPP 1

#include <vector>
#include "fermi_table.hpp"
#include "initial_data.hpp"

using std::vector;

vector<double> create_pressure_reference(const FermiTable& eos,
					 const InitialData& id);

#endif // CREATE_PRESSURE_REFERENCE_HPP

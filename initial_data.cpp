#include "initial_data.hpp"
#include "vector_io.hpp"
#include "vector_utils.hpp"
#include "get_composition_data.hpp"

InitialData::InitialData(const string& radius_file,
			 const string& density_file,
			 const string& temperature_file,
			 const string& velocity_file):
  radius_list(load_txt(radius_file)),
  radius_mid(mid_array(radius_list)),
  density_list(decapitate(load_txt(density_file))),
  temperature_list(decapitate(load_txt(temperature_file))),
  velocity_list(load_txt(velocity_file)),
  tracers_list(get_composition_data()) {}

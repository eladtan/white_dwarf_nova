#include "units.hpp"
#include "sim_data.hpp"
#include "my_main_loop.hpp"

using namespace std;

int main(void)
{
  const Units units;
  const InitialData id
    ("radius_list.txt",
     "density_list.txt",
     "temperature_list.txt",
     "velocity_list.txt");
  SimData sim_data(id,
		   units,
		   CircularSection(id.radius_mid.front(),
				   id.radius_mid.back(),
				   0.46*M_PI,
				   0.54*M_PI));
  hdsim& sim = sim_data.getSim();
  my_main_loop(sim,sim_data.getEOS());
  return 0;
}

#include "units.hpp"
#include "sim_data.hpp"
#include "my_main_loop.hpp"

using namespace std;

int main(void)
{
  Units units;
  SimData sim_data(InitialData("radius_list.txt",
			       "density_list.txt",
			       "temperature_list.txt",
			       "velocity_list.txt"),
		   units);
  hdsim& sim = sim_data.getSim();
  my_main_loop(sim,sim_data.getEOS());
  return 0;
}

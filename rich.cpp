#include <ctime>
#include "units.hpp"
#include "sim_data.hpp"
#include "my_main_loop.hpp"

using namespace std;

int main(void)
{
  const clock_t begin = clock();
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

  const clock_t end = clock();
  ofstream f("wall_time.txt");
  f << static_cast<double>(end-begin)/CLOCKS_PER_SEC << endl;
  f.close();

  return 0;
}

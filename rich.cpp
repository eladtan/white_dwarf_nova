#include "units.hpp"
#include "sim_data.hpp"
#include "my_main_loop.hpp"
#include "nuclear_burn.hpp"

using namespace std;

namespace {
  map<string,double> create_tracers(void)
  {
    map<string,double> res;
    res["He4"] = 1;
    res["C12"] = 0;
    res["O16"] = 0;
    res["Ne20"] = 0;
    res["Mg24"] = 0;
    res["Si28"] = 0;
    res["S32"] = 0;
    res["Ar36"] = 0;
    res["Ca40"] = 0;
    res["Ti44"] = 0;
    res["Cr48"] = 0;
    res["Fe52"] = 0;
    res["Ni56"] = 0;
    return res;
  }
}

int main(void)
{
  /*
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
  */

  FermiTable eos("eos_tab.coded",1,1,0,generate_atomic_properties());
  const map<string,double> tracers = create_tracers();
  NuclearBurn nb(string("alpha_table"),
		 string("ghost"),
		 eos);
  const double density = 1.9745e6;
  const double temperature = 2.0116e8;
  const double e = eos.dt2e(density,
			    temperature,
			    tracers);
  const double q = nb.calcEnergyDepositionRate
    (density,
     temperature,
     tracers);
  cout << e << " " << q << endl;
  return 0;
}

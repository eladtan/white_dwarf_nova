#include "my_main_loop.hpp"
#include "temperature_appendix.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include "write_cycle.hpp"
#include "source/newtonian/test_2d/multiple_diagnostics.hpp"
#include "nuclear_burn.hpp"
#include "multiple_manipulation.hpp"
#include "equilibrium_retouch.hpp"

using namespace simulation2d;

void my_main_loop(hdsim& sim, const FermiTable& eos)
{
  write_snapshot_to_hdf5(sim,"initial.h5",
			 vector<DiagnosticAppendix*>
			 (1,new TemperatureAppendix(eos)));
  const double tf = 10;
  SafeTimeTermination term_cond(tf, 1e6);
  vector<DiagnosticFunction*> diag_list = VectorInitialiser<DiagnosticFunction*>()
    [new ConsecutiveSnapshots(new ConstantTimeInterval(tf/1000),
			      new Rubric("snapshot_",".h5"),
			      vector<DiagnosticAppendix*>
			      (1,new TemperatureAppendix(eos)))]
    [new WriteTime("time.txt")]
    [new WriteCycle("cycle.txt")]
    ();
  MultipleDiagnostics diag(diag_list);
  /*
  NuclearBurn manip(string("alpha_table"),
		    string("ghost"),
		    eos);
  */
  MultipleManipulation manip
    (VectorInitialiser<Manipulate*>(new NuclearBurn(string("alpha_table"),
				       string("ghost"),
				       eos))
     (new EquilibriumRetouch)());
    main_loop(sim,
	    term_cond,
	    &hdsim::TimeAdvance,
	    &diag,
	    &manip);
  write_snapshot_to_hdf5(sim,"final.h5",
			 vector<DiagnosticAppendix*>
			 (1,new TemperatureAppendix(eos)));
}

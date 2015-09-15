#include "write_cycle.hpp"
#include "source/misc/simple_io.hpp"

WriteCycle::WriteCycle(const string& fname):
  fname_(fname) {}

void WriteCycle::operator()(const hdsim& sim)
{
  write_number(sim.getCycle(),fname_);
}

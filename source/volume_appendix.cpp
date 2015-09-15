#include "volume_appendix.hpp"

VolumeAppendix::VolumeAppendix(void) {}

string VolumeAppendix::getName(void) const
{
  return "volumes";
}

vector<double> VolumeAppendix::operator()(const hdsim& sim) const
{
  vector<double> res;
  const CacheData& cd = sim.getCacheData();
  for(size_t i=0;i<cd.volumes.size();++i)
    res.push_back(cd.volumes[i]);
  return res;
}

#include "generate_atomic_properties.hpp"

boost::container::flat_map<string,pair<double,double> > generate_atomic_properties(void)
{
  boost::container::flat_map<string,pair<double,double> > res;
  res["He4"] = pair<double,double>(4,2);
  res["C12"] = pair<double,double>(12,6);
  res["O16"] = pair<double,double>(16,8);
  res["Ne20"] = pair<double,double>(20,10);
  res["Mg24"] = pair<double,double>(24,12);
  res["Si28"] = pair<double,double>(28,14);
  res["S32"] = pair<double,double>(32,16);
  res["Ar36"] = pair<double,double>(36,18);
  res["Ca40"] = pair<double,double>(40,20);
  res["Ti44"] = pair<double,double>(44,22);
  res["Cr48"] = pair<double,double>(48,44);
  res["Fe52"] = pair<double,double>(52,26);
  res["Ni56"] = pair<double,double>(56,28);
  return res;
}

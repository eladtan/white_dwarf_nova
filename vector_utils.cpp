#include "vector_utils.hpp"
#include <cassert>

using std::size_t;

bool is_strictly_increasing(const vector<double>& v)
{
  for(size_t i=1;i<v.size();++i){
    if(v.at(i)<=v.at(i-1))
      return false;
  }
  return true;
}

vector<double> decapitate(const vector<double>& v)
{
  assert(v.size()>0);
  vector<double> res(v.size()-1);
  for(size_t i=0;i<res.size();++i)
    res.at(i) = v.at(i);
  return res;
}

vector<double> mid_array(const vector<double>& v)
{
  assert(v.size()>0);
  vector<double> res(v.size()-1);
  for(size_t i=0;i<res.size();++i)
    res.at(i) = 0.5*(v.at(i)+v.at(i+1));
  return res;
}

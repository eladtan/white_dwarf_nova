#include "vector_io.hpp"
#include <fstream>

using std::ifstream;
using std::ofstream;

vector<double> load_txt(const string& fname)
{
  vector<double> res;

  ifstream f(fname.c_str());
  assert(f);
  double buf = 0;
  while(f>>buf)
    res.push_back(buf);
  f.close();
  return res;
}

void save_txt(const string& fname, const vector<double>& v)
{
  ofstream f(fname.c_str());
  assert(f);
  const size_t v_size = v.size();
  for(size_t i=0;i<v_size;++i)
    f << v.at(i) << "\n";
  f.close();
}

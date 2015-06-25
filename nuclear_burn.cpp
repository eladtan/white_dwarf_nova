#include "nuclear_burn.hpp"

extern "C" {
  void burn_step_(int* indexeos,
		  double* density,
		  double* energy,
		  double* tburn,
		  double* xn,
		  double* atomw,
		  double* atomn,
		  double* dedtmp,
		  int* matters,
		  double* dt,
		  double* qrec,
		  int* nse,
		  double* tmp_nse,
		  int* key_done,
		  char* screen_type);
}

namespace {

  pair<double,vector<double> > burn_step_wrapper(double density,
						 double energy,
						 double tburn,
						 vector<double> xn,
						 pair<double,double> az,
						 double dt)
  {
    int indexeos = 0;
    double dedtmp = 0;
    int matters = static_cast<int>(xn.size());
    double qrec = 0;
    int nse = 0;
    double tmp_nse = 1e10;
    char screen_type[80] = "default";
    int key_done = 0;
    burn_step_(&indexeos,
	       &density,
	       &energy,
	       &tburn,
	       &xn[0],
	       &az.first,
	       &az.second,
	       &dedtmp,
	       &matters,
	       &dt,
	       &qrec,
	       &nse,
	       &tmp_nse,
	       &key_done,
	       screen_type);
    assert(key_done==1);
    return pair<double,vector<double> >(qrec,xn);
  }

  template<class S, class T> const T& safe_retrieve(const map<S,T>& m,
						    const S& s)
  {
    /*
    map<S,T>::const_iterator it = m.find(s);
    assert(it!=m.end());
    return it->second;
    */
    assert(m.find(s)!=m.end());
    return m.find(s)->second;
  }

  pair<double,double> calc_average_atomic_properties
  (const map<string,pair<double,double> >& atomic_properties,
   const map<string,double>& compositions)
  {
    double total = 0;
    double aa = 0;
    double zz = 0;
    for(map<string,pair<double,double> >::const_iterator it=
	  atomic_properties.begin();
	it!=atomic_properties.end();
	++it){
      const double mass_frac = safe_retrieve(compositions,it->first);
      const double A = it->second.first;
      const double Z = it->second.second;
      total += mass_frac;
      aa += mass_frac/A;
      zz += mass_frac*Z/A;
    }
    return pair<double,double>(total/aa,zz/aa);
  }
}

NuclearBurn::NuclearBurn
(const map<string,pair<double,double> >& atomic_properties,
 const string& ignore_label):
  atomic_properties_(atomic_properties),
  t_prev_(0),
  ignore_label_(ignore_label) {}

void NuclearBurn::operator()(hdsim& sim)
{
  const double dt = sim.getTime() - t_prev_;
  vector<ComputationalCell>& cells = sim.getAllCells();
  for(size_t i=0;i<cells.size();++i){
    ComputationalCell& cell = cells[i];
    if(safe_retrieve(cell.stickers,ignore_label_))
      continue;
  }
  sim.recalculateExtensives();
}

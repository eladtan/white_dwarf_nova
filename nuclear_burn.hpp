#ifndef NUCLEAR_BURN_HPP
#define NUCLEAR_BURN_HPP 1

class NuclearBurn
{
public:
  NuclearBurn(const map<string,pair<double,double> >& atomic_propeties,
	      const string& ignore_label);

  void operator()(hdsim& sim);
  
private:

  const map<string,pair<double,double> > atomic_properties_;
  mutable double t_prev_;
  const string ignore_label_;
  const FermiTable& eos_;
};

#endif // NUCLEAR_BURN_HPP

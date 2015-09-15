#ifndef INTERPOLATOR_HPP
#define INTERPOLATOR_HPP 1

#include <vector>

using std::vector;

class Interpolator
{
public:

  Interpolator(const vector<double>& x_list,
	       const vector<double>& y_list);

  double operator()(double x) const;

private:
  const vector<double> x_list_;
  const vector<double> y_list_;
};

#endif // INTERPOLATOR_HPP

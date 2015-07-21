#include "interpolator.hpp"
#include "vector_utils.hpp"
#include <cassert>

using std::size_t;

Interpolator::Interpolator(const vector<double>& x_list,
			   const vector<double>& y_list):
  x_list_(x_list), y_list_(y_list)
{
  assert(is_strictly_increasing(x_list_));
  assert(x_list_.size()==y_list_.size());
}

double Interpolator::linear(double x) const
{
  const double x_list_front = x_list_.front();
  const double x_list_back = x_list_.back();
  assert(x>x_list_front);
  assert(x<x_list_back);       
  //      assert(x>x_list_.front());
  //      assert(x<x_list_.back());
  for(size_t i=1;i<x_list_.size();++i){
    if(x_list_.at(i)>x)
      return y_list_.at(i-1) + (y_list_.at(i)-y_list_.at(i-1))*
	(x-x_list_.at(i-1))/(x_list_.at(i)-x_list_.at(i-1));
  }
  throw "point outside bound";
}

namespace {
  double hsf(double x)
  {
    if(x>0)
      return 1;
    if(x<0)
      return 0;
    return 0.5;
  }
}

double Interpolator::flat(double x) const
{
  const double x_list_front = x_list_.front();
  const double x_list_back = x_list_.back();
  assert(x>x_list_front);
  assert(x<x_list_back);       
  //      assert(x>x_list_.front());
  //      assert(x<x_list_.back());
  for(size_t i=1;i<x_list_.size();++i){
    if(x_list_.at(i)>x)
      return y_list_.at(i-1) + (y_list_.at(i)-y_list_.at(i-1))*
	hsf(x-0.5*(x_list_.at(i)+x_list_.at(i-1)));
  }
  throw "point outside bound";
}

double Interpolator::operator()(double x, bool t) const
{
  if(t)
    return linear(x);
  else
    return flat(x);
}

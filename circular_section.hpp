#ifndef CIRCULAR_SECTION_HPP
#define CIRCULAR_SECTION_HPP 1

#include "source/tessellation/shape_2d.hpp"

using std::pair;

class CircularSection: public Shape2D
{
public:

  CircularSection(const double radius_in,
		  const double radius_out,
		  const double angle_left,
		  const double angle_right);

  bool operator()(const Vector2D& r) const;

  pair<double,double> getAngles(void) const;

private:
  const double radius_in_;
  const double radius_out_;
  const double angle_left_;
  const double angle_right_;
};

#endif // CIRCULAR_SECTION_HPP

#include <cmath>
#include "circular_section.hpp"

CircularSection::CircularSection
(const double radius_in,
 const double radius_out,
 const double angle_left,
 const double angle_right):
  radius_in_(radius_in),
  radius_out_(radius_out),
  angle_left_(angle_left),
  angle_right_(angle_right) {}

bool CircularSection::operator()(const Vector2D& r) const
{
  const double radius = abs(r);
  const double angle = atan2(r.y,r.x);
  return (radius_in_<radius &&
	  radius_out_>radius &&
	  angle_left_<angle &&
	  angle_right_>angle);
}

pair<double,double> CircularSection::getAngles(void) const
{
  return pair<double,double>(angle_left_, angle_right_);
}

pair<double,double> CircularSection::getRadii(void) const
{
  return pair<double,double>(radius_in_, radius_out_);
}

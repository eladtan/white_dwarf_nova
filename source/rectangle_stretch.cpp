#include "rectangle_stretch.hpp"

pair<Vector2D, Vector2D> rectangle_stretch
(const pair<Vector2D,Vector2D>& boundaries,
 double ratio)
{
  const Vector2D centre = 0.5*(boundaries.first + boundaries.second);
  return pair<Vector2D,Vector2D>(centre+ratio*(boundaries.first-centre),
				 centre+ratio*(boundaries.second-centre));
}

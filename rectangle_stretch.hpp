#ifndef RECTANGLE_STRETCH_HPP
#define RECTANGLE_STRETCH_HPP 1

#include <utility>
#include "source/tessellation/geometry.hpp"

using std::pair;

pair<Vector2D, Vector2D> rectangle_stretch
(const pair<Vector2D,Vector2D>& boundaries,
 double ratio);

#endif // RECTANGLE_STRETCH_HPP

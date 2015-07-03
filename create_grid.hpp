#ifndef CREATE_GRID_HPP
#define CREATE_GRID_HPP 1

#include <vector>
#include "source/tessellation/geometry.hpp"

using std::vector;
using std::pair;

vector<Vector2D> create_grid(const pair<Vector2D,Vector2D>& boundaries,
			     const double dq,
			     const double r_min);

#endif // CREATE_GRID_HPP

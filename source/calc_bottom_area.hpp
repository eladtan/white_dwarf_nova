#ifndef CALC_BOTTOM_AREA_HPP
#define CALC_BOTTOM_AREA_HPP 1

#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/tessellation/shape_2d.hpp"

double calc_bottom_area(const Tessellation& tess,
			const Shape2D& shape,
			const PhysicalGeometry& pg);

#endif // CALC_BOTTOM_AREA_HPP

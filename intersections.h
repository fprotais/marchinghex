#pragma once
#include <array>
#include <ultimaille/geometry.h>

using namespace UM; 
namespace intersections {
	bool point_is_in_tet(const vec3& A, const vec3& B, const vec3& C, const vec3& D, const vec3& P, std::array<double, 4>& l, const double dillatation = 1E-8);
	bool point_is_in_triangle(const vec2& A, const vec2& B, const vec2& C, const vec2& P, std::array<double, 3>& l);
}






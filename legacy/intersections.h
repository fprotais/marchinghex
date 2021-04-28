#pragma once
#include <array>
#include <ultimaille/vec.h>

namespace intersections {
	bool point_is_in_tet(const UM::vec3& A, const UM::vec3& B, const UM::vec3& C, const UM::vec3& D, const UM::vec3& P, std::array<double, 4>& l, const double dillatation = 1E-8);
	bool point_is_in_triangle(const UM::vec2& A, const UM::vec2& B, const UM::vec2& C, const UM::vec2& P, std::array<double, 3>& l);
	bool inter_sec_triangle(const UM::vec3& A, const UM::vec3& B, const UM::vec3& C, const UM::vec3& s1, const UM::vec3& s2, double& lA, double& lB, double& lC);
	
	/*
		C---D
		|	|
		|	|
		A---B
	*/
	bool inter_sec_quad(const UM::vec3& A, const UM::vec3& B, const UM::vec3& C, const UM::vec3& D, const UM::vec3& s1, const UM::vec3& s2, UM::vec3& inter);
	bool point_is_cube(const UM::vec3& O, const UM::vec3& dO, const UM::vec3& p);
}






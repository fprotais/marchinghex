#include "intersections.h"

using namespace UM; 


namespace intersections {
	bool point_is_in_tet(const vec3& A, const vec3& B, const vec3& C, const vec3& D, const vec3& P, std::array<double, 4>& l, const double dillatation) {
		constexpr int find[4][3] = { { 1, 2, 3 },{ 3, 2, 0 },{ 0, 1, 3 },{ 0, 2, 1 } };
		vec3 center = 0.25 * (A + B + C + D);
		double dillat = dillatation * (A - B).norm();
		center *= dillat;
		std::array<vec3, 4> p = { { (1 + dillat) * A - center, (1 + dillat) * B - center, (1 + dillat) * C - center, (1 + dillat) * D - center } };

		double total_volume = dot(p[3] - p[0], cross(p[1] - p[0], p[2] - p[0]));

		FOR(f, 4) {
			const vec3 vA = p[find[f][0]];
			const vec3 vB = p[find[f][1]];
			const vec3 vC = p[find[f][2]];

			l[f] = dot(vC - P, cross(vA - P, vB - P)) / total_volume;
		}
		return (l[0] >= 0 && l[1] >= 0 && l[2] >= 0 && l[3] >= 0);

	}

	inline double det(const vec2& v1, const  vec2& v2) {
		return v1.x * v2.y - v1.y * v2.x;
	}
	bool point_is_in_triangle(const vec2& A, const vec2& B, const vec2& C, const vec2& P, std::array<double, 3>& l) {
		double area = det(B - A, C - A);
		l[0] = det(C - B, P - B) / area;
		l[1] = det(A - C, P - C) / area;
		l[2] = det(B - A, P - A) / area;

		return l[0]>= 0 && l[1]>=0 && l[2]>=0;
	}

}



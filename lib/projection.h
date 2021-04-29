#include <ultimaille/all.h>
#include <array>

double point_segment_squared_distance(const UM::vec3& P, const std::array<UM::vec3, 2>& AB, UM::vec3& closest_point, std::array<double, 2>& l);

double point_triangle_squared_distance(const UM::vec3& P, const std::array<UM::vec3, 3>& ABC, UM::vec3& closest_point, std::array<double, 3>& l);

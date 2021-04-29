#pragma once
#include <ultimaille/all.h>

double hausdorff_dist(const std::vector<UM::vec3>& P, const std::vector<std::array<UM::vec3, 3>>& mesh, std::vector<double>& dist2mesh);

double hausdorff_dist(const UM::Hexahedra& hex, const UM::Tetrahedra& tet, int tri_sampling = 5);

double compute_scaled_jacobian(const UM::Hexahedra& hex, UM::CellAttribute<double>& cell_min_sj);

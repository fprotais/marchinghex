#pragma once
#include <array>
#include <vector>
#include <ultimaille/all.h>

namespace marchingcube {

	void hexify(const std::vector<UM::vec3>& verts, const std::vector<int>& tets, std::vector<UM::vec3>& mc_verts, std::vector<int>& mc_hexes);

	// moved output in argument to decrease allocation
	void extract_polyedra_dual(const std::array<bool, 8>& corner_is_in, std::array<int, 24>& he_array);
	struct Polyedra {
		std::vector<std::vector<int>> facets_;
		std::vector<UM::vec3> verts_;
	};
	struct Triangulated_Polyedra {
		std::vector<std::array<int, 3>> triangles_;
		std::vector<UM::vec3> verts_;
		std::vector<std::array<int, 3>> triangles_adjacent_;
		std::vector<int> triangles_original_facets_;
	};
	struct Hexadreized_Polyedra {
		std::vector<std::array<int, 8>> hexaedra_;
		std::vector<UM::vec3> verts_;
	};
	void extract_polyedra(const std::array<int, 24>& he_array, const std::array<bool, 8>& corner_is_in, const std::array <UM::vec3, 8>& corners, const std::array<UM::vec3, 24>& corner_bnd_alternative, Polyedra& P);

	void compute_polyedra_convex_hull(const Polyedra& polyedra, Triangulated_Polyedra& TP);

	UM::vec3 center_of_polyedra_facet_with_flattening(const Triangulated_Polyedra& TP, const int facet);

	void extract_hexes_from_polyedra(const Triangulated_Polyedra& TPol, const Polyedra& P, Hexadreized_Polyedra& HP);

}

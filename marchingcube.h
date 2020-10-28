#pragma once
#include <array>
#include <vector>
#include <ultimaille/geometry.h>
#include <ultimaille/volume.h>

using namespace UM;
namespace marchingcube {
	//void hexify(const Tetrahedra& m, Hexahedra& hex);

	void hexify(const std::vector<vec3>& verts, const std::vector<int>& tets, std::vector<vec3>& mc_verts, std::vector<int>& mc_hexes);

	std::array<int, 24> extract_polyedra_dual(const std::array<bool, 8>& corner_is_in);
	struct Polyedra {
		std::vector<std::vector<int>> facets_;
		std::vector<vec3> verts_;
	};
	struct Triangulated_Polyedra {
		std::vector<std::array<int, 3>> triangles_;
		std::vector<vec3> verts_;
		std::vector<std::array<int, 3>> triangles_adjacent_;
		std::vector<int> triangles_original_facets_;
	};
	struct Hexadreized_Polyedra {
		std::vector<std::array<int, 8>> hexaedra_;
		std::vector<vec3> verts_;
	};
	Polyedra extract_polyedra(const std::array<size_t, 24>& he_array, const std::array<bool, 8>& corner_is_in, const std::array <vec3, 8>& corners, const std::array<vec3, 24>& corner_bnd_alternative);
		
	Triangulated_Polyedra compute_polyedra_convex_hull(const Polyedra& polyedra);

	vec3 center_of_polyedra_facet_with_flattening(const Triangulated_Polyedra& TP, const int facet);

	Hexadreized_Polyedra extract_hexes_from_polyedra(const Triangulated_Polyedra& TPol, const Polyedra& P);





}



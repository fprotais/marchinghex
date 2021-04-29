#pragma once
#include <array>
#include <vector>
#include <ultimaille/all.h>


namespace Marchinghex {
	struct Polyedra {
		std::vector<std::vector<int>> facets_;
		std::vector<UM::vec3> verts_;
		std::vector<int> vertex_id_;
	};
	struct Triangulated_Polyedra {
		std::vector<std::array<int, 3>> triangles_;
		std::vector<UM::vec3> verts_;
		std::vector<std::array<int, 3>> triangles_adjacent_;
		std::vector<int> triangles_original_facets_;
	};
	enum VERT_TYPE { VERT_IS_INSIDE = 0, VERT_IS_ON_BOUNDARY, VERT_IS_ON_FEATURE, VERT_IS_FEATURE_POINT };

	struct Hexadreized_Polyedra {
		std::vector<std::array<int, 8>> hexaedra_;
		std::vector<UM::vec3> verts_;
		std::vector<VERT_TYPE> vert_type_;
	}; 
	void extract_polyedra_dual(const std::array<bool, 8>& corner_is_in, std::array<int, 24>& he_array);

	void extract_polyedra(const std::array<int, 24>& he_array, const std::array<bool, 8>& corner_is_in, const std::array <UM::vec3, 8>& corners, const std::array<UM::vec3, 24>& corner_bnd_alternative, Polyedra& P);

	void compute_polyedra_convex_hull(const Polyedra& polyedra, Triangulated_Polyedra& TP);

	UM::vec3 center_of_polyedra_facet_with_flattening(const Triangulated_Polyedra& TP, const int facet);

	static const UM::vec3 NO_WISH(1E10, 1E10, 1E10);
	static const std::array<UM::vec3, 7> DEFAULT_NO_WISH = { NO_WISH,NO_WISH,NO_WISH,NO_WISH,NO_WISH,NO_WISH,NO_WISH };
	void extract_hexes_from_polyedra(const Triangulated_Polyedra& TPol, const Polyedra& P, Hexadreized_Polyedra& hexes, const std::array<UM::vec3, 7>& wishes = DEFAULT_NO_WISH);
}


class grid_pattern_extractor {
public:
	grid_pattern_extractor(const UM::Hexahedra& grid)
	{
		grid_.points.data->assign(grid.points.data->begin(), grid.points.data->end());
		grid_.cells.assign(grid.cells.begin(), grid.cells.end());
		init();

	}
	grid_pattern_extractor(const UM::Tetrahedra& m, const UM::Hexahedra& grid)
		: grid_pattern_extractor(grid) {
		set_domain(m);
	}
	void set_domain(const UM::Tetrahedra& m);
	void set_boundary(const UM::Triangles& bnd);
	void set_hard_edges(const UM::PolyLine& hd, bool detect_features_points = true);
	void set_feature_points(const std::vector<UM::vec3>& points);
	// only extract inside hexes, not parallel
	void extract_patterns(UM::Hexahedra& hex);
	void extract_patterns(UM::Hexahedra& hex, std::vector<bool>& is_in, bool do_projection = true);
	void extract_patterns(UM::Hexahedra& hex, std::vector<bool>& is_in, std::vector<UM::vec3>& vert_wish, std::vector<Marchinghex::VERT_TYPE>& wish_type);

	int debug_ = 1;
private:
	void init();
	void grid_pre_processing();
	void extract(const int h, const std::array<bool, 8>& c_is_in, Marchinghex::Hexadreized_Polyedra& hexes, bool proj);

	UM::Hexahedra grid_;
	std::vector<std::array<bool, 8>> vert_is_in_;
	std::vector<std::array<UM::vec3, 7>> wishes_;
	std::vector<std::array<UM::vec3, 24>>  bnd_intersections_;
	std::vector<std::array<UM::vec3, 24>>  default_bnd_intersections_;

};
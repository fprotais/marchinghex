#pragma once
#include "utils/basic.h"
#include <array>
#include <vector>
namespace marchingcube {
	void hexify(const GEO::Mesh& m, GEO::Mesh& hex);


	std::array<size_t, 24> extract_polyedra_dual(const std::array<bool, 8>& corner_is_in);
	struct Polyedra {
		std::vector<std::vector<size_t>> facets_;
		std::vector<GEO::vec3> verts_;
	};
	struct Triangulated_Polyedra {
		std::vector<std::array<size_t, 3>> triangles_;
		std::vector<GEO::vec3> verts_;
		std::vector<std::array<size_t, 3>> triangles_adjacent_;
		std::vector<size_t> triangles_original_facets_;
	};
	struct Hexadreized_Polyedra {
		std::vector<std::array<size_t, 8>> hexaedra_;
		std::vector<GEO::vec3> verts_;
	};
	Polyedra extract_polyedra(const std::array<size_t, 24>& he_array, const std::array<bool, 8>& corner_is_in, const std::array <GEO::vec3, 8>& corners, const std::array<GEO::vec3, 24>& corner_bnd_alternative);
		
	Triangulated_Polyedra compute_polyedra_convex_hull(const Polyedra& polyedra);

	GEO::vec3 center_of_polyedra_facet_with_flattening(const Triangulated_Polyedra& TP, const size_t facet);

	Hexadreized_Polyedra extract_hexes_from_polyedra(const Triangulated_Polyedra& TPol, const Polyedra& P);





}


#include "utils/MeshGeo.h"
#include "polycube_util.h"
namespace Hexdom {
	inline void init_config_marching_cube_on_a_cell(std::array <GEO::vec3, 8>& corners,  std::array<GEO::vec3, 24>& corner_bnd_alternative) {
		corners = { vec3(0,0,0), vec3(1,0,0),vec3(0,1,0), vec3(1,1,0), vec3(0,0,1), vec3(1,0,1), vec3(0,1,1), vec3(1,1,1) };
		double shift = 0.5;
		corner_bnd_alternative = {
			vec3(shift , 0, 0),vec3(0, shift, 0),vec3(0, 0, shift),
			vec3(1 - shift, 0, 0),vec3(1, shift, 0),vec3(1, 0, shift),
			vec3(shift, 1, 0),vec3(0, 1 - shift, 0),vec3(0, 1, shift),
			vec3(1 - shift, 1, 0),vec3(1, 1 - shift, 0),vec3(1, 1, shift),
			vec3(shift, 0, 1),vec3(0, shift, 1),vec3(0, 0, 1 - shift),
			vec3(1 - shift, 0, 1),vec3(1, shift, 1),vec3(1, 0, 1 - shift),
			vec3(shift, 1, 1),vec3(0, 1 - shift, 1),vec3(0, 1, 1 - shift),
			vec3(1 - shift, 1, 1),vec3(1, 1 - shift, 1),vec3(1, 1, 1 - shift),
		};

	}
	inline void test_marching_cube_on_a_cell(std::array<bool, 8> config) {
		std::cerr << "-> ";
		FOR(i, 8) std::cerr << config[i] << " ";
		std::cerr << std::endl;
		std::array <GEO::vec3, 8> corners; 
		std::array<GEO::vec3, 24> corner_bnd_alternative;
		init_config_marching_cube_on_a_cell(corners, corner_bnd_alternative);
		
		std::array<size_t, 24> he_array = marchingcube::extract_polyedra_dual(config);
		marchingcube::Polyedra P = marchingcube::extract_polyedra(he_array, config, corners, corner_bnd_alternative);
		GEO::Mesh poly;
		poly.vertices.create_vertices((index_t)P.verts_.size());
		FOR(v, poly.vertices.nb()) X(poly)[v] = P.verts_[v];
		FOR(f, P.facets_.size()) {
			poly.facets.create_facets(1, (index_t)P.facets_[f].size());
			FOR(fv, P.facets_[f].size()) poly.facets.set_vertex(f, fv, (index_t)P.facets_[f][fv]);
		}
		DROP(poly, "polyedra");
		marchingcube::Triangulated_Polyedra TP = marchingcube::compute_polyedra_convex_hull(P);
		GEO::Mesh Tpoly;
		Tpoly.vertices.create_vertices((index_t)TP.verts_.size());
		FOR(v, Tpoly.vertices.nb()) X(Tpoly)[v] = TP.verts_[v];
		Tpoly.facets.create_triangles((index_t)TP.triangles_.size());
		FOR(f, Tpoly.facets.nb()) {
			FOR(fv, 3) Tpoly.facets.set_vertex(f, fv, (index_t)TP.triangles_[f][fv]);
		}
		GEO::Attribute<index_t> or_facet(Tpoly.facets.attributes(), "facet");
		FOR(f, Tpoly.facets.nb()) {
			or_facet[f] = (index_t)TP.triangles_original_facets_[f];
		}
		DROP(Tpoly, "triangulated");
		marchingcube::Hexadreized_Polyedra HP = marchingcube::extract_hexes_from_polyedra(TP, P);
		GEO::Mesh hexPoly;
		hexPoly.vertices.create_vertices((index_t)HP.verts_.size());
		FOR(v, hexPoly.vertices.nb()) X(hexPoly)[v] = HP.verts_[v];
		hexPoly.cells.create_hexes((index_t)HP.hexaedra_.size());
		FOR(h, hexPoly.cells.nb()) FOR(hv, 8)
			hexPoly.cells.set_vertex(h, hv, (index_t)HP.hexaedra_[h][hv]);
		DROP(hexPoly, "Hexadreized");

	}
}


#include <iostream>
#include <cstdlib>
#include <ultimaille/all.h>


#include <chrono>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include "marchingcube.h"



using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)

inline void add_2_mesh(UM::Polygons & m, const Marchinghex::Polyedra & P) {
	int off_v = m.points.create_points(P.verts_.size());
	FOR(v, P.verts_.size()) m.points[off_v + v] = P.verts_[v];
	FOR(f, P.facets_.size()) {
		int off_c = m.create_facets(1, P.facets_[f].size());
		FOR(i, P.facets_[f].size()) m.vert(off_c, i) = off_v + P.facets_[f][i];
	}
}
inline void add_2_mesh(UM::Triangles& m, const Marchinghex::Triangulated_Polyedra& CH, FacetAttribute<int>& ch_facet_id) {
	int off_v = m.points.create_points(CH.verts_.size());
	FOR(v, CH.verts_.size()) m.points[off_v + v] = CH.verts_[v];
	int off_c = m.create_facets(CH.triangles_.size());
	FOR(f, CH.triangles_.size()) {
		FOR(i, 3) m.vert(off_c + f, i) = off_v + CH.triangles_[f][i];
		ch_facet_id[off_c + f] = CH.triangles_original_facets_[f];
	}
}
inline void add_2_mesh(UM::Hexahedra& m, const Marchinghex::Hexadreized_Polyedra& H) {
	int off_v = m.points.create_points(H.verts_.size());
	FOR(v, H.verts_.size()) m.points[off_v + v] = H.verts_[v];
	int off_c = m.create_cells(H.hexaedra_.size());
	FOR(c, H.hexaedra_.size()) {
		FOR(i, 8) m.vert(off_c + c, i) = off_v + H.hexaedra_[c][i];
	}
}


static void init_config_marching_cube_on_a_cell(std::array <vec3, 8>&corners, std::array<vec3, 24>&corner_bnd_alternative) {
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

void add_marching_hex_config_to_meshes(const std::array<bool, 8>& config, Hexahedra& m_hexes, Polygons& m_poly, Triangles& m_convex_hull, FacetAttribute<int>& ch_facet_id, const vec3& shift) {

	std::array<vec3, 8> corners;
	std::array<vec3, 24> corner_bnd_alternative;
	init_config_marching_cube_on_a_cell(corners, corner_bnd_alternative);
	FOR(i, 8) corners[i] += shift;
	FOR(i, 24) corner_bnd_alternative[i] += shift;


	std::array<int, 24> he_array;
	Marchinghex::extract_polyedra_dual(config, he_array);
	Marchinghex::Polyedra polyedra;
	Marchinghex::extract_polyedra(he_array, config, corners, corner_bnd_alternative, polyedra);
	Marchinghex::Triangulated_Polyedra polyedra_convex_hull;
	Marchinghex::compute_polyedra_convex_hull(polyedra, polyedra_convex_hull);
	Marchinghex::Hexadreized_Polyedra hexes;
	Marchinghex::extract_hexes_from_polyedra(polyedra_convex_hull, polyedra, hexes);

	add_2_mesh(m_poly, polyedra);

	add_2_mesh(m_convex_hull, polyedra_convex_hull, ch_facet_id);

	add_2_mesh(m_hexes, hexes);

}



static void print_all_config(const std::string& ext = ".geogram") {
	Hexahedra allhexes;
	Polygons allpolyedra;
	Triangles allconvexhull;
	FacetAttribute<int> ch_facet_id(allconvexhull);
	std::filesystem::create_directory("allhex");
	std::filesystem::create_directory("allpoly");
	std::filesystem::create_directory("allhull");
	FOR(i7, 2) FOR(i6, 2) FOR(i5, 2) FOR(i4, 2) FOR(i3, 2) FOR(i2, 2) FOR(i1, 2) FOR(i0, 2) {
		std::array<bool, 8> config = { (bool)i0,(bool)i1,(bool)i2,(bool)i3,(bool)i4,(bool)i5,(bool)i6,(bool)i7 };
		int num_config= 0;
		FOR(i, 8) if (config[i]) num_config += 1 << i;
		FOR(i, 8) std::cerr << config[i] << " ";
		std::cerr << "\t| number :" << num_config << std::endl;
		vec3 shift(8 * i3 + 4 * i2 + 2 * i1 + i0, 8 * i7 + 4 * i6 + 2 * i5 + i4, 0);
		shift = 1.5 * shift;
		Hexahedra lochexes;
		Polygons locpolyedra;
		Triangles locconvexhull;
		FacetAttribute<int> loc_facet_id(locconvexhull);

		add_marching_hex_config_to_meshes(config, lochexes, locpolyedra, locconvexhull, loc_facet_id, { 0,0,0 });
		add_marching_hex_config_to_meshes(config, allhexes, allpolyedra, allconvexhull, ch_facet_id, shift);
		write_by_extension("allhex/hexes_c" + std::to_string(num_config) + ext, lochexes);
		write_by_extension("allpoly/polyedra_c" + std::to_string(num_config) + ext, locpolyedra);
		write_by_extension("allhull/convexhulls_c" + std::to_string(num_config) + ext, locconvexhull, { {},{ {"id", loc_facet_id.ptr} }, {}});

	}
	write_by_extension("hexes" + ext, allhexes);
	write_by_extension("polyedra" + ext, allpolyedra);
	write_by_extension("convexhulls" + ext, allconvexhull, { {},{ {"id", ch_facet_id.ptr} }, {} });
}






void  extract(const std::array<vec3, 8>& hex, const std::array<bool, 8>& c_is_in, std::array<vec3, 24>& corner_bnd_alternative, Marchinghex::Hexadreized_Polyedra& hexes, const std::string& name = "extract") {
	std::array<int, 24> he_array;
	Marchinghex::extract_polyedra_dual(c_is_in, he_array);
	Marchinghex::Polyedra polyedra;
	Marchinghex::Triangulated_Polyedra polyedra_convex_hull;

	Marchinghex::extract_polyedra(he_array, c_is_in, hex, corner_bnd_alternative, polyedra);
	Marchinghex::compute_polyedra_convex_hull(polyedra, polyedra_convex_hull);
	Marchinghex::extract_hexes_from_polyedra(polyedra_convex_hull, polyedra, hexes, Marchinghex::DEFAULT_NO_WISH);

}

inline int bin_to_int(int x, int y, int z) {
	int res = 0;
	res += (x % 2) << 0;
	res += (y % 2) << 1;
	res += (z % 2) << 2;
	return res;
}
inline std::array<bool, 8> inverse_cell(const std::array<bool, 8>& is_in) {
	std::array<bool, 8> res;
	FOR(i, 8) res[i] = !is_in[i];
	return res;
}


void visualise_splitting(const std::array<bool, 8>& config, Hexahedra& cube, Hexahedra& mc_m, Hexahedra& comp_mc_m, const vec3& shift) {
	
	std::vector<std::array<vec3, 24>> cba;
	std::vector<std::array<bool, 8>> configs;
	std::vector<std::array<bool, 3>> splits;

	{
		std::array<vec3, 8> hex;
		std::array<vec3, 24> corner_bnd_alternative;

		init_config_marching_cube_on_a_cell(hex, corner_bnd_alternative);
		FOR(i, 8) hex[i] += shift;
		FOR(i, 24) corner_bnd_alternative[i] += shift;
		cube.points.create_points(8);
		FOR(i, 8) {
			cube.points[i] = hex[i];
		}
		cube.create_cells(1);
		FOR(i, 8) {
			cube.vert(0, i) = i;
		}
		cba.push_back(corner_bnd_alternative);
		configs.push_back(config);
		splits.push_back({ 1,1,1 });
	}

	constexpr int facet_vertex[6][4] = { {0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, {0,1,2,3}, {4,5,6,7} };
	constexpr int d_jump[3][3] = {
	{1,2,4},
	{2,1,4},
	{4,1,2},
	};

	FOR(dim_cut, 3) {

		FOR(h, cube.ncells()) if (splits[h][dim_cut]) {
			splits[h][dim_cut] = false;
			std::array<vec3, 8> hex;
			FOR(hv, 8) hex[hv] = cube.points[cube.vert(h, hv)];

			std::array<vec3, 8> hex1, hex2;
			std::array<bool, 8> config1, config2;
			std::array<vec3, 24> cba1, cba2;

			FOR(i, 4) {
				int o_i_1 = facet_vertex[2 * dim_cut][i];
				int o_i_2 = facet_vertex[2 * dim_cut + 1][i];
				hex1[o_i_1] = hex[o_i_1];
				config1[o_i_1] = configs[h][o_i_1];
				hex2[o_i_2] = hex[o_i_2];
				config2[o_i_2] = configs[h][o_i_2];
				FOR(d, 3) cba1[3 * o_i_1 + d] = cba[h][3 * o_i_1 + d];
				FOR(d, 3) cba2[3 * o_i_2 + d] = cba[h][3 * o_i_2 + d];
			}

			FOR(d2, 2) FOR(d1, 2) {

				int i = d2 * d_jump[dim_cut][1] + d1 * d_jump[dim_cut][2];
				int fi = d1 * 2 + d2;
				int o_i_1 = facet_vertex[2 * dim_cut][fi];
				int o_i_2 = facet_vertex[2 * dim_cut + 1][fi];
				vec3 new_v = 0.5 * (hex[o_i_1] + hex[o_i_2]);
				if (configs[h][o_i_2] && !configs[h][o_i_1]) new_v = 0.5 * (hex[o_i_2] + cba[h][3 * o_i_2 + dim_cut]);
				if (!configs[h][o_i_2] && configs[h][o_i_1]) new_v = 0.5 * (hex[o_i_1] + cba[h][3 * o_i_1 + dim_cut]);
				hex1[d_jump[dim_cut][0] + i] = new_v;
				hex2[i] = new_v;
				bool has_in = configs[h][o_i_2] || configs[h][o_i_1];
				config1[d_jump[dim_cut][0] + i] = has_in;
				config2[i] = has_in;

				FOR(d, 3) cba1[3 * (d_jump[dim_cut][0] + i) + d] = vec3(1e100, 1e100, 1e100);
				FOR(d, 3) cba2[3 * (i)+d] = vec3(1e100, 1e100, 1e100);

				cba1[3 * (d_jump[dim_cut][0] + i) + dim_cut] = cba[h][3 * o_i_2 + dim_cut];
				cba2[3 * (i)+dim_cut] = cba[h][3 * o_i_1 + dim_cut];
			}
			FOR(z, 2) FOR(y, 2) FOR(x, 2) FOR(d, 3) {

				int s[3] = { 0,0,0 }; s[d] = 1;
				if (cba1[3 * bin_to_int(x, y, z) + d][0] == 1e100)
					cba1[3 * bin_to_int(x, y, z) + d] = 0.5 * (hex1[bin_to_int(x, y, z)]
						+ hex1[bin_to_int(x + s[0], y + s[1], z + s[2])]);
				if (cba2[3 * bin_to_int(x, y, z) + d][0] == 1e100)
					cba2[3 * bin_to_int(x, y, z) + d] = 0.5 * (hex2[bin_to_int(x, y, z)]
						+ hex2[bin_to_int(x + s[0], y + s[1], z + s[2])]);
			}

			int off_set_vert = cube.points.create_points(16);
			FOR(hv, 8) {
				cube.points[off_set_vert + hv] = hex1[hv];
				cube.points[off_set_vert + 8 + hv] = hex2[hv];
			}
			int off_set_cells = cube.create_cells(1);
			FOR(hv, 8) {
				cube.vert(h, hv) = off_set_vert + hv;
				cube.vert(off_set_cells, hv) = off_set_vert + 8 + hv;
			}
			cba[h] = cba1;
			cba.push_back(cba2);
			configs[h] = config1;
			configs.push_back(config2);
			splits.push_back(splits[h]);
		}
	}

	FOR(h, cube.ncells()) {
		Marchinghex::Hexadreized_Polyedra mc_hexes;
		std::array<vec3, 8> hex;
		FOR(hv, 8) hex[hv] = cube.points[cube.vert(h, hv)];
		extract(hex, configs[h], cba[h], mc_hexes);
		add_2_mesh(mc_m, mc_hexes);
		extract(hex, inverse_cell(configs[h]), cba[h], mc_hexes);
		add_2_mesh(comp_mc_m, mc_hexes);

	}
}


void visualise_splitting(const std::string& ext = ".mesh") {
	Hexahedra allhexes_pri;
	Hexahedra allhexes_compl;

	std::filesystem::create_directory("split");
	FOR(i7, 2) FOR(i6, 2) FOR(i5, 2) FOR(i4, 2) FOR(i3, 2) FOR(i2, 2) FOR(i1, 2) FOR(i0, 2) {
		std::array<bool, 8> config = { (bool)i0,(bool)i1,(bool)i2,(bool)i3,(bool)i4,(bool)i5,(bool)i6,(bool)i7 };
		int num_config = 0;
		FOR(i, 8) if (config[i]) num_config += 1 << i;
		FOR(i, 8) std::cerr << config[i] << " ";
		std::cerr << "\t| number :" << num_config << std::endl;
		vec3 shift(8 * i3 + 4 * i2 + 2 * i1 + i0, 8 * i7 + 4 * i6 + 2 * i5 + i4, 0);
		shift = 1.5 * shift;
		Hexahedra lochexes_pri;
		Hexahedra lochexes_dual;
		Hexahedra lochexes_cube;
		Hexahedra tmp_cube;


		visualise_splitting(config, lochexes_cube, lochexes_pri, lochexes_dual, { 0,0,0 });
		visualise_splitting(config, tmp_cube, allhexes_pri, allhexes_compl, shift);
		write_by_extension("split/hexes_pri_c" + std::to_string(num_config) + ext, lochexes_pri);
		write_by_extension("split/hexes_compl_c" + std::to_string(num_config) + ext, lochexes_dual);
		write_by_extension("split/grid_c" + std::to_string(num_config) + ext, lochexes_cube);


	}
	write_by_extension("split_hexes_pri" + ext, allhexes_pri);
	write_by_extension("split_hexes_compl" + ext, allhexes_compl);


}

int main(int argc, char** argv) {
	print_all_config(".mesh");
	visualise_splitting(".mesh");

    return 0;
}

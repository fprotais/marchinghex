#include "marchingcube.h"
#include "intersections.h"
#include <map>
#include <chrono>


#define FOR(i, n) for(int i = 0; i < n; i++)
static constexpr int NOT_AN_ID = -1;
using namespace UM;
namespace Marchinghex {
	static void renumber_facet(std::vector<int>& facet, const std::vector<int>& vert_id, const int facet_number) {
		constexpr int facet_order[6][12] = {
			{0,2 + 8,14 + 8,4,13 + 8,19 + 8,6,20 + 8,8 + 8,2,7 + 8,1 + 8},
			{1,4 + 8,10 + 8,3,11 + 8,23 + 8,7,22 + 8,16 + 8,5,17 + 8,5 + 8},

			{0,0 + 8,3 + 8,1,5 + 8,17 + 8,5,15 + 8,12 + 8,4,14 + 8,2 + 8},
			{2,8 + 8,20 + 8,6,18 + 8,21 + 8,7,23 + 8,11 + 8,3,9 + 8,6 + 8},

			{0,1 + 8,7 + 8,2,6 + 8,9 + 8,3,10 + 8,4 + 8,1,3 + 8,0 + 8},
			{4,12 + 8,15 + 8,5,16 + 8,22 + 8,7,21 + 8,18 + 8,6,19 + 8,13 + 8}
		};
		std::vector<int> new_facet(facet.size());
		int pos = 0;
		FOR(i, 12) {
			FOR(fv, facet.size()) if (vert_id[facet[fv]] == facet_order[facet_number][i]) {
				new_facet[pos++] = facet[fv];
				break;
			}
		}
		facet.clear();
		facet.assign(new_facet.begin(), new_facet.end());
	}

	inline std::pair<int, int> ordered_pair(const int v1, const int v2) {
		if (v1 > v2) return { v1, v2 };
		return { v2, v1 };
	}

	inline int prec_on_dual(const std::array<int, 24>& he_array, const int he) {
		if (he > 24) return NOT_AN_ID;
		if (he_array[he] == NOT_AN_ID) return NOT_AN_ID;
		int curr_he = he;
		int i = 0;
		while (he_array[curr_he] != he) {
			curr_he = he_array[curr_he];
			if (i++ > 24) return NOT_AN_ID;
		}
		return curr_he;
	}
	void extract_polyedra_dual(const std::array<bool, 8>& corner_is_in, std::array<int, 24>& he_array) {
		he_array = {
			2,0,1,
			4,5,3,
			7,8,6,
			11,9,10,
			13,14,12,
			17,15,16,
			20,18,19,
			22,23,21
		};
		// corner a, corner b, he a, he b
		constexpr int edges[12][4] = {
			{0,1,0,3},
			{2,3,6,9},
			{4,5,12,15},
			{6,7,18,21},
			{0,2,1,7},
			{1,3,4,10},
			{4,6,13,19},
			{5,7,16,22},
			{0,4,2,14},
			{1,5,5,17},
			{2,6,8,20},
			{3,7,11,23},
		};

		FOR(i, 12) if (!corner_is_in[edges[i][0]] && !corner_is_in[edges[i][1]]) {
			int a = prec_on_dual(he_array, edges[i][2]);
			int b = prec_on_dual(he_array, edges[i][3]);
			he_array[a] = he_array[edges[i][3]];
			he_array[b] = he_array[edges[i][2]];
			he_array[edges[i][2]] = NOT_AN_ID;
			he_array[edges[i][3]] = NOT_AN_ID;
		}
	}

	void extract_polyedra(const std::array<int, 24>& he_array, const std::array<bool, 8>& corner_is_in, const std::array<vec3, 8>& corners, const std::array<vec3, 24>& corner_bnd_alternative, Polyedra& polyedra) {

		std::vector<bool> mark(24, false);
		std::vector<int> vert_id;
		std::vector<std::array<int, 3>> vert_facets;
		constexpr std::array<int, 3> array_vert_facets[8] = { {0,2,4}, {1,2,4}, {0,3,4}, {1,3,4}, {0,2,5}, {1,2,5}, {0,3,5}, {1,3,5} };
		int last_facet = 5;
		FOR(he, 24) if ((!mark[he]) && (he_array[he] != NOT_AN_ID)) {
			int c = he / 3;
			if (corner_is_in[c]) {
				polyedra.verts_.push_back(corners[c]);
				vert_facets.push_back(array_vert_facets[c]);
				mark[he + 1] = true;
				mark[he + 2] = true;
				vert_id.push_back(c);
			}
			else {
				int curr_he = he;
				last_facet++;
				do {
					mark[curr_he] = true;
					polyedra.verts_.push_back(corner_bnd_alternative[curr_he]);
					vert_id.push_back(8 + curr_he);
					int c_local = curr_he / 3;
					int he_local = curr_he % 3;
					vert_facets.push_back({ last_facet, 0, 0 });
					vert_facets.back()[1] = array_vert_facets[c_local][(he_local + 1) % 3];
					vert_facets.back()[2] = array_vert_facets[c_local][(he_local + 2) % 3];
					curr_he = he_array[curr_he];
				} while (curr_he != he);
			}

		}


		polyedra.facets_.resize(last_facet + 1);
		FOR(v, polyedra.verts_.size()) FOR(i, 3) polyedra.facets_[vert_facets[v][i]].push_back(v);
		FOR(f, 6) renumber_facet(polyedra.facets_[f], vert_id, f);
		polyedra.vertex_id_.assign(vert_id.begin(), vert_id.end());
	}

	void compute_polyedra_convex_hull(const Polyedra& polyedra, Triangulated_Polyedra& triangles) {

		triangles.verts_.resize(polyedra.verts_.size());
		FOR(v, polyedra.verts_.size()) triangles.verts_[v] = polyedra.verts_[v];
		// build an initial triangulation
		FOR(f, polyedra.facets_.size()) {
			if (polyedra.facets_[f].size() < 3) continue;
			int off_t = (int)triangles.triangles_.size();
			triangles.triangles_.resize(off_t + polyedra.facets_[f].size() - 2);
			triangles.triangles_original_facets_.resize(off_t + polyedra.facets_[f].size() - 2);
			triangles.triangles_adjacent_.resize(off_t + polyedra.facets_[f].size() - 2, { NOT_AN_ID,NOT_AN_ID ,NOT_AN_ID });
			FOR(t, polyedra.facets_[f].size() - 2) {
				triangles.triangles_[off_t + t][0] = polyedra.facets_[f][0];
				triangles.triangles_[off_t + t][1] = polyedra.facets_[f][t + 1];
				triangles.triangles_[off_t + t][2] = polyedra.facets_[f][t + 2];
				if (t > 0) triangles.triangles_adjacent_[off_t + t][0] = off_t + t - 1;
				if (t < polyedra.facets_[f].size() - 3) triangles.triangles_adjacent_[off_t + t][2] = off_t + t + 1;
				triangles.triangles_original_facets_[off_t + t] = f;
			}
		}

		// we do edge swaps until it is convex
		bool good = true;
		do {
			good = true;
			FOR(i, 3 * triangles.triangles_.size()) {
				int f = i / 3;
				int fe = i % 3;
				int f2 = triangles.triangles_adjacent_[f][fe];
				if (f2 == NOT_AN_ID) continue;
				vec3 e = triangles.verts_[triangles.triangles_[f][(fe + 1) % 3]]
					- triangles.verts_[triangles.triangles_[f][fe]];
				vec3 n1 = cross(
					triangles.verts_[triangles.triangles_[f][1]] - triangles.verts_[triangles.triangles_[f][0]],
					triangles.verts_[triangles.triangles_[f][2]] - triangles.verts_[triangles.triangles_[f][0]]
				);
				vec3 n2 = cross(
					triangles.verts_[triangles.triangles_[f2][1]] - triangles.verts_[triangles.triangles_[f2][0]],
					triangles.verts_[triangles.triangles_[f2][2]] - triangles.verts_[triangles.triangles_[f2][0]]
				);
				if (e * cross(n1, n2) < -1E-8) {
					int tri1[3];
					FOR(li, 3) tri1[li] = triangles.triangles_[f][(fe + li) % 3];
					int fe2 = triangles.triangles_[f2][0] == tri1[0] ? 2 : triangles.triangles_[f2][1] == tri1[0] ? 0 : 1;
					int tri2[3];
					FOR(li, 3) tri2[li] = triangles.triangles_[f2][(fe2 + li) % 3];
					triangles.triangles_[f][0] = tri1[0];
					triangles.triangles_[f][1] = tri2[2];
					triangles.triangles_[f][2] = tri1[2];
					triangles.triangles_[f2][0] = tri2[0];
					triangles.triangles_[f2][1] = tri1[2];
					triangles.triangles_[f2][2] = tri2[2];

					// yepeeeeee
					int f3 = triangles.triangles_adjacent_[f][(fe + 1) % 3];
					int f4 = triangles.triangles_adjacent_[f][(fe + 2) % 3];
					int f5 = triangles.triangles_adjacent_[f2][(fe2 + 1) % 3];
					int f6 = triangles.triangles_adjacent_[f2][(fe2 + 2) % 3];

					triangles.triangles_adjacent_[f][0] = f5;
					triangles.triangles_adjacent_[f][1] = f2;
					triangles.triangles_adjacent_[f][2] = f4;
					triangles.triangles_adjacent_[f2][0] = f3;
					triangles.triangles_adjacent_[f2][1] = f;
					triangles.triangles_adjacent_[f2][2] = f6;

					if (f3 != NOT_AN_ID) FOR(fe3, 3) if (triangles.triangles_adjacent_[f3][fe3] == f) triangles.triangles_adjacent_[f3][fe3] = f2;
					if (f4 != NOT_AN_ID) FOR(fe4, 3) if (triangles.triangles_adjacent_[f4][fe4] == f) triangles.triangles_adjacent_[f4][fe4] = f;
					if (f5 != NOT_AN_ID) FOR(fe5, 3) if (triangles.triangles_adjacent_[f5][fe5] == f2) triangles.triangles_adjacent_[f5][fe5] = f;
					if (f6 != NOT_AN_ID) FOR(fe6, 3) if (triangles.triangles_adjacent_[f6][fe6] == f2) triangles.triangles_adjacent_[f6][fe6] = f2;
					good = false;
					break;

				}
			}
		} while (!good);
	}

	vec3 center_of_polyedra_facet_with_flattening(const Triangulated_Polyedra& TP, const int Pfacet) {
		std::vector<vec2> new_pos(TP.verts_.size());
		//if (Pfacet < 6) 
		{
			std::vector<bool> mark_v(TP.verts_.size(), false);
			FOR(f, TP.triangles_.size()) if (TP.triangles_original_facets_[f] == Pfacet)
				FOR(fv, 3) mark_v[TP.triangles_[f][fv]] = true;
			int nb_vert = 0;
			vec3 res;
			FOR(v, TP.verts_.size()) if (mark_v[v]) {
				res += TP.verts_[v];
				nb_vert++;
			}
			return res / nb_vert;
		}

		std::vector<int> pile;
		std::vector<bool> facet_is_plan(TP.triangles_.size(), false);
		std::vector<bool> vertex_is_plan(TP.verts_.size(), false);
		//plan 1st triangle
		{
			int f0 = NOT_AN_ID;
			FOR(t, TP.triangles_.size()) if (TP.triangles_original_facets_[t] == Pfacet) {
				f0 = t;
				break;
			}
			assert(f0 != NOT_AN_ID);
			facet_is_plan[f0] = true;
			vec3 e1 = TP.verts_[TP.triangles_[f0][1]] - TP.verts_[TP.triangles_[f0][0]];
			vec3 e2 = TP.verts_[TP.triangles_[f0][2]] - TP.verts_[TP.triangles_[f0][0]];

			double theta = std::acos(e1 * e2 / (e2.norm() * e1.norm()));

			new_pos[TP.triangles_[f0][0]] = vec2(0, 0);
			new_pos[TP.triangles_[f0][1]] = e1.norm() * vec2(1, 0);
			new_pos[TP.triangles_[f0][2]] = e2.norm() * vec2(std::cos(theta), std::sin(theta));

			FOR(fv, 3) vertex_is_plan[TP.triangles_[f0][fv]] = true;
			FOR(fc, 3) {
				int f2 = TP.triangles_adjacent_[f0][fc];
				if (f2 != NOT_AN_ID && TP.triangles_original_facets_[f2] == Pfacet)
					pile.push_back(f2);
			}
		}

		while (!pile.empty()) {
			int f1 = pile.back(); pile.pop_back();
			if (facet_is_plan[f1]) continue;
			facet_is_plan[f1] = true;
			int fv_no_plannar = 2;
			FOR(fv, 2) if (!vertex_is_plan[TP.triangles_[f1][fv]]) fv_no_plannar = fv;
			vec3 e1 = TP.verts_[TP.triangles_[f1][(fv_no_plannar + 2) % 3]] - TP.verts_[TP.triangles_[f1][(fv_no_plannar + 1) % 3]];
			vec3 e2 = TP.verts_[TP.triangles_[f1][(fv_no_plannar + 0) % 3]] - TP.verts_[TP.triangles_[f1][(fv_no_plannar + 1) % 3]];
			vec2 n1 = new_pos[TP.triangles_[f1][(fv_no_plannar + 2) % 3]] - new_pos[TP.triangles_[f1][(fv_no_plannar + 1) % 3]];
			double theta = std::acos(e1 * e2 / (e2.norm() * e1.norm()));

			vec2 n2 = vec2(n1.x * std::cos(theta) - n1.y * std::sin(theta), n1.x * std::sin(theta) + n1.y * std::cos(theta));
			n2 = n2 * e2.norm() / e1.norm();
			new_pos[TP.triangles_[f1][fv_no_plannar]] = new_pos[TP.triangles_[f1][(fv_no_plannar + 1) % 3]] + n2;

			vertex_is_plan[TP.triangles_[f1][fv_no_plannar]] = true;
			FOR(fc, 3) {
				int f2 = TP.triangles_adjacent_[f1][fc];
				if (f2 == NOT_AN_ID) continue;
				if (facet_is_plan[f2]) continue;
				if (TP.triangles_original_facets_[f2] == Pfacet)
					pile.push_back(f2);
			}

		}
		vec2 plannar_center;
		int nb_verts = 0;
		FOR(v, TP.verts_.size()) if (vertex_is_plan[v]) {
			plannar_center += new_pos[v];
			nb_verts++;
		}
		plannar_center = plannar_center / nb_verts;

		FOR(t, TP.triangles_.size()) if (TP.triangles_original_facets_[t] == Pfacet) {
			std::array<double, 3> l = { 0,0,0 };
			bool is_in = intersections::point_is_in_triangle(
				new_pos[TP.triangles_[t][0]],
				new_pos[TP.triangles_[t][1]],
				new_pos[TP.triangles_[t][2]],
				plannar_center,
				l
			);
			if (is_in) {
				vec3 center;
				FOR(fv, 3) center += l[fv] * TP.verts_[TP.triangles_[t][fv]];
				return center;
			}

		}

		// in case center is outside 
		std::vector<bool> mark_v(TP.verts_.size(), false);
		FOR(f, TP.triangles_.size()) if (TP.triangles_original_facets_[f] == Pfacet)
			FOR(fv, 3) mark_v[TP.triangles_[f][fv]] = true;
		int nb_vert = 0;
		vec3 res;
		FOR(v, TP.verts_.size()) if (mark_v[v]) {
			res += TP.verts_[v];
			nb_vert++;
		}
		return res / nb_vert;

	}

	inline bool is_wish(const vec3& v) {
		return !(v.x == NO_WISH.x && v.y == NO_WISH.y && v.z == NO_WISH.z);
	}
	void extract_hexes_from_polyedra(const Triangulated_Polyedra& TPol, const Polyedra& P, Hexadreized_Polyedra& hexes, const std::array<vec3, 7>& wishes) {
		hexes.hexaedra_.resize(TPol.verts_.size());

		const int nb_pol_verts = (int)TPol.verts_.size();

		int nb_pol_facets = (int)P.facets_.size();
		if (TPol.triangles_.size() == 0) return;

		std::vector<int> connected_composant;
		int nb_of_composant = 0;
		{
			DisjointSet ds((int)P.verts_.size());
			FOR(f, P.facets_.size()) FOR(fv, P.facets_[f].size())
				ds.merge((int)P.facets_[f][fv], (int)(P.facets_[f][(fv + 1) % P.facets_[f].size()]));
			nb_of_composant = ds.get_sets_id(connected_composant);
		}
		hexes.verts_.resize(nb_pol_verts + nb_pol_facets + nb_of_composant);
		FOR(v, nb_pol_verts) hexes.verts_[v] = TPol.verts_[v];
		hexes.vert_type_.resize(hexes.verts_.size(), VERT_IS_INSIDE);
		for (int i = 0; i < nb_pol_verts; i++) {
			if (P.vertex_id_[i] > 7) hexes.vert_type_[i] = VERT_IS_ON_BOUNDARY;
		}
		std::vector<std::vector<int>> P_to_TP_facets(P.facets_.size());
		FOR(t, TPol.triangles_.size()) P_to_TP_facets[TPol.triangles_original_facets_[t]].push_back(t);
		std::vector<vec3> absolute_center(nb_of_composant);
		bool has_feature_corner = false;
		std::vector<int> nb_of_center(nb_of_composant, 0);
		FOR(f, P.facets_.size()) {
			if (P.facets_[f].size() < 3) continue;
			if (f>5) hexes.vert_type_[nb_pol_verts + nb_of_composant + f] = VERT_IS_ON_BOUNDARY;

			vec3 center = center_of_polyedra_facet_with_flattening(TPol, f);
			if (f > 5 && is_wish(wishes[0])) {
				if (nb_of_composant < 2) {
					center = wishes[0];
					hexes.vert_type_[nb_pol_verts + nb_of_composant + f] = VERT_IS_FEATURE_POINT;
					has_feature_corner = true; 
				}
				else std::cerr << "A hard corner is unmatchable due to resolution" << std::endl;
			}
			hexes.verts_[nb_pol_verts + nb_of_composant + f] = center;
			absolute_center[connected_composant[P.facets_[f][0]]] += center; nb_of_center[connected_composant[P.facets_[f][0]]]++;
		}
		FOR(i, nb_of_composant) assert(nb_of_center[i] != 0);
		FOR(i, nb_of_composant) absolute_center[i] = absolute_center[i] / nb_of_center[i];
		FOR(i, nb_of_composant) hexes.verts_[nb_pol_verts + i] = absolute_center[i];

		std::map<std::pair<int, int>, int> edge_mid_vertices;
		std::vector<bool> vertex_is_hardedge(hexes.verts_.size(), false);


		{
			std::vector<vec3> new_verts;
			std::vector<VERT_TYPE> new_verts_type;
			std::map<std::pair<int, int>, vec3> edge_mid_vertices_wish;
			std::map<std::pair<int, int>, VERT_TYPE> edge_mid_vertices_type;

			FOR(f, 6/*P.facets_.size()*/) FOR(fv, P.facets_[f].size()) {
				const int v1 = P.facets_[f][fv];
				const int v2 = P.facets_[f][(fv + 1) % P.facets_[f].size()];
				bool has_wish = false;
				if (P.vertex_id_[v1] > 7 && P.vertex_id_[v2] > 7) {
					if (is_wish(wishes[f + 1])) {
						has_wish = true;
						edge_mid_vertices_wish[ordered_pair(v1, v2)] = wishes[f + 1];
						edge_mid_vertices_type[ordered_pair(v1, v2)] = VERT_IS_ON_FEATURE;
					}
					else {
						edge_mid_vertices_type[ordered_pair(v1, v2)] = VERT_IS_ON_BOUNDARY;
					}
				}
				else {
					edge_mid_vertices_type[ordered_pair(v1, v2)] = VERT_IS_INSIDE;
				}
				edge_mid_vertices[ordered_pair(v1, v2)] = NOT_AN_ID;
			}
			for (auto& edge_mid_vertex : edge_mid_vertices) {
				edge_mid_vertex.second = (int)hexes.verts_.size() + (int)new_verts.size();
				new_verts_type.push_back(edge_mid_vertices_type[edge_mid_vertex.first]);
				if (edge_mid_vertices_wish.find(edge_mid_vertex.first) == edge_mid_vertices_wish.end()) {
					new_verts.push_back(0.5 * (P.verts_[edge_mid_vertex.first.first] + P.verts_[edge_mid_vertex.first.second]));
					vertex_is_hardedge.push_back(false);
				}
				else {
					new_verts.push_back(edge_mid_vertices_wish[edge_mid_vertex.first]);
					vertex_is_hardedge.push_back(true);
				}
			}
			hexes.verts_.resize(nb_pol_verts + nb_pol_facets + nb_of_composant + new_verts.size());
			FOR(v, new_verts.size()) hexes.verts_[nb_pol_verts + nb_pol_facets + nb_of_composant + v] = new_verts[v];
			hexes.vert_type_.resize(nb_pol_verts + nb_pol_facets + nb_of_composant + new_verts.size());
			FOR(v, new_verts.size()) hexes.vert_type_[nb_pol_verts + nb_pol_facets + nb_of_composant + v] = new_verts_type[v];
			


		}
		bool setting_middle_point_on_feature = false;
		if (!has_feature_corner) FOR(f, P.facets_.size()) if (f > 5) {
			std::vector<vec3> set_vertices;
			FOR(fv, P.facets_[f].size()) {
				std::pair<int, int> p = ordered_pair(P.facets_[f][fv], P.facets_[f][(fv + 1) % P.facets_[f].size()]);
				if (vertex_is_hardedge[edge_mid_vertices[p]])
					set_vertices.push_back(hexes.verts_[edge_mid_vertices[p]]);
			}
			if (set_vertices.size() > 1) {
				setting_middle_point_on_feature = true;
				vec3 new_center;
				FOR(i, set_vertices.size()) new_center += set_vertices[i];
				new_center = new_center / (double)set_vertices.size();
				hexes.verts_[nb_pol_verts + nb_of_composant + f] = new_center;
				hexes.vert_type_[nb_pol_verts + nb_of_composant + f] = VERT_IS_ON_FEATURE;

			}
		}

		std::vector<std::array<int, 3>> vert_neigh_verts(P.verts_.size(), { NOT_AN_ID, NOT_AN_ID, NOT_AN_ID });
		std::vector<std::array<int, 3 >> vert_neigh_facets(P.verts_.size(), { NOT_AN_ID, NOT_AN_ID, NOT_AN_ID });
		FOR(f, P.facets_.size()) FOR(fv, P.facets_[f].size()) {
			FOR(i, 3) if (vert_neigh_verts[P.facets_[f][fv]][i] == -1) {
				vert_neigh_verts[P.facets_[f][fv]][i] = P.facets_[f][(fv + 1) % P.facets_[f].size()];
				vert_neigh_facets[P.facets_[f][fv]][i] = f;
				break;
			}
		}
		std::vector<std::vector<int>> facet_of_vertex(nb_pol_verts);
		FOR(f, nb_pol_facets) FOR(fv, P.facets_[f].size()) facet_of_vertex[P.facets_[f][fv]].push_back(f);
		FOR(v, nb_pol_verts) {

			assert(facet_of_vertex[v].size() == 3);

			// renumbering so that everything is oriented well
			std::array<int, 3> vecfacets = vert_neigh_facets[v];
			std::array<int, 3> other_vertex = vert_neigh_verts[v];
			bool v3_is_in_f1 = false;
			FOR(fv, P.facets_[vecfacets[0]].size()) if (P.facets_[vecfacets[0]][fv] == other_vertex[2]) v3_is_in_f1 = true;
			if (!v3_is_in_f1) {
				std::swap(vecfacets[1], vecfacets[2]);
				std::swap(other_vertex[1], other_vertex[2]);
			}

			hexes.hexaedra_[v][0] = v;
			hexes.hexaedra_[v][1] = edge_mid_vertices[ordered_pair(v, other_vertex[0])];
			hexes.hexaedra_[v][2] = edge_mid_vertices[ordered_pair(v, other_vertex[1])];
			hexes.hexaedra_[v][3] = nb_pol_verts + nb_of_composant + vecfacets[1];
			hexes.hexaedra_[v][4] = edge_mid_vertices[ordered_pair(v, other_vertex[2])];
			hexes.hexaedra_[v][5] = nb_pol_verts + nb_of_composant + vecfacets[0];
			hexes.hexaedra_[v][6] = nb_pol_verts + nb_of_composant + vecfacets[2];
			hexes.hexaedra_[v][7] = nb_pol_verts + connected_composant[v];

		}
	}

}

static void simplifying_vertices(Hexahedra& hex, double tol = 1E-6) {
	//removing isolated verts
	std::vector<bool> tokill(hex.nverts(), true);
	FOR(h, hex.ncells()) FOR(hv,8)  tokill[hex.vert(h,hv)] = false;
	hex.delete_vertices(tokill);

	// colocating vertices 
	std::vector<int> old2new;
	std::vector<vec3> vecarray(hex.points.begin(), hex.points.end());
	colocate(vecarray, old2new, tol);
	DisjointSet ds(hex.nverts());
	FOR(v, hex.nverts()) ds.merge(v, old2new[v]);
	int nb_group = ds.get_sets_id(old2new);

	std::vector<vec3> new_verts(nb_group);
	FOR(v, hex.nverts()) new_verts[old2new[v]] = hex.points[v];
	hex.points.data->clear();
	hex.points.data->assign(new_verts.begin(), new_verts.end());

	FOR(h, hex.ncells())  FOR(hc, 8) hex.vert(h, hc) = old2new[hex.vert(h, hc)];
}



inline int bin_to_int(int x, int y, int z) {
	int res = 0;
	res += (x % 2) << 0;
	res += (y % 2) << 1;
	res += (z % 2) << 2;
	return res;
}
void  grid_pattern_extractor::init() {
	vert_is_in_.resize(grid_.ncells(), { { 0,0,0,0,0,0,0,0 } });
	wishes_.resize(grid_.ncells(), Marchinghex::DEFAULT_NO_WISH);
	bnd_intersections_.resize(grid_.ncells());
	default_bnd_intersections_.resize(grid_.ncells());
	FOR(c, grid_.ncells()) {
		FOR(z, 2) FOR(y, 2) FOR(x, 2) FOR(d, 3) {
			int s[3] = { 0,0,0 }; s[d] = 1;
			default_bnd_intersections_[c][3 * bin_to_int(x, y, z) + d] = 0.5 * (grid_.points[grid_.vert(c, bin_to_int(x, y, z))]
				+ grid_.points[grid_.vert(c, bin_to_int(x + s[0], y + s[1], z + s[2]))]);
			bnd_intersections_[c][3 * bin_to_int(x, y, z) + d] = default_bnd_intersections_[c][3 * bin_to_int(x, y, z) + d];
		}
	}

}

void grid_pattern_extractor::set_domain(const UM::Tetrahedra& m) {
	std::cerr << "Matching domain with hexes..." << std::endl;

	std::vector<BBox3> boxes(grid_.nverts());
	FOR(v, grid_.nverts()) {
		boxes[v].add(grid_.points[v]); 
	} 
	HBoxes looker(boxes);
	int nb_cells_done = 0;
	#pragma omp parallel
	{
		std::vector<bool> vert_is_in(grid_.nverts(), false);
		#pragma omp for nowait
		FOR(c, m.ncells()) {
			//#pragma omp atomic
			//++nb_cells_done;
			//if (debug_) if (nb_cells_done % 1000 == 0) {
			//	#pragma omp critical
			//	std::cerr << " Matching cells		 :		c = " << (nb_cells_done / 1000) * 1000 << "    /    " << m.ncells() << std::endl;
			//}

			vec3 A = m.points[m.vert(c, 0)];
			vec3 B = m.points[m.vert(c, 1)];
			vec3 C = m.points[m.vert(c, 2)];
			vec3 D = m.points[m.vert(c, 3)];
			BBox3 box;
			FOR(cv, 4) box.add(m.points[m.vert(c, cv)]);
			box.dilate(1E-6);
			std::vector<int> intersect;
			looker.intersect(box, intersect);
			for (int i : intersect) {
				std::array<double, 4> l;
				if (intersections::point_is_in_tet(A, B, C, D, grid_.points[i], l)) {
					vert_is_in[i] = true;
				}
			}
		}
		#pragma omp critical
		{
			FOR(h, grid_.ncells()) FOR(cv, 8) {
				vert_is_in_[h][cv] = vert_is_in_[h][cv] || vert_is_in[grid_.vert(h, cv)];
			}
		}

	}
		
	


}
inline std::pair<int, int> orderedpair(int a, int b) {
	if (a > b) return { a, b };
	else return { b, a };
}
void grid_pattern_extractor::set_boundary(const UM::Triangles& bnd) {
	std::cerr << "Matching boundary with hexes..." << std::endl;
	constexpr int hexedges[12][2] = {
		{0,1}, {2,3}, {4,5}, {6,7},
		{0,2}, {1,3}, {4,6}, {5,7},
		{0,4}, {1,5}, {2,6}, {3,7},
	};
	constexpr int edge_dim[12] = { 0,0,0,0,1,1,1,1,2,2,2,2 };
	std::vector<std::pair<int, int>> edge;
	std::map<std::pair<int, int>, vec3> edge_inter;
	FOR(h, grid_.ncells()) FOR(he, 12) edge_inter[orderedpair(grid_.vert(h, hexedges[he][0]), grid_.vert(h, hexedges[he][1]))] = 0.5*(grid_.points[grid_.vert(h, hexedges[he][0])] + grid_.points[grid_.vert(h, hexedges[he][1])]);
	std::vector<BBox3> boxes(edge_inter.size());

	int id = 0;
	for (auto p : edge_inter) {
		edge.push_back(p.first);
		boxes[id].add(grid_.points[p.first.first]);
		boxes[id++].add(grid_.points[p.first.second]);
	}

	HBoxes looker(boxes);
	FOR(t, bnd.nfacets()) {
		if (debug_)
			if (t % 1000 == 0)std::cerr << " Matching triangles		 :		t = " << t << "    /    " << bnd.nfacets() << std::endl;
		vec3 A = bnd.points[bnd.vert(t, 0)];
		vec3 B = bnd.points[bnd.vert(t, 1)];
		vec3 C = bnd.points[bnd.vert(t, 2)];
		BBox3 box;
		FOR(tv, 3) box.add(bnd.points[bnd.vert(t, tv)]);
		box.dilate(1E-6);

		std::vector<int> intersect;
		looker.intersect(box, intersect);
		for (int i : intersect) {
			vec3 X1 = grid_.points[edge[i].first];
			vec3 X2 = grid_.points[edge[i].second];

			double l1,l2,l3;
			if (intersections::inter_sec_triangle(A, B, C, X1, X2, l1, l2, l3)) {
				vec3 bndcut = l1 * A + l2 * B + l3 * C;

				if ((bndcut - X1).norm() < 0.1) bndcut = 0.9 * X1 + 0.1 * X2;
				if ((bndcut - X2).norm() < 0.1) bndcut = 0.1 * X1 + 0.9 * X2;

				edge_inter[edge[i]] = bndcut;
			}
		}
	}
	FOR(h, grid_.ncells()) {
		FOR(he, 12) {
			int l1 = hexedges[he][0];
			int l2 = hexedges[he][1];
			int d = edge_dim[he];
			bnd_intersections_[h][3 * l1 + d] = edge_inter[orderedpair(grid_.vert(h, l1), grid_.vert(h, l2))];
			bnd_intersections_[h][3 * l2 + d] = edge_inter[orderedpair(grid_.vert(h, l1), grid_.vert(h, l2))];
		}
	}
}

static void extract_feature_points(const UM::PolyLine& hd, std::vector<vec3>& points) {
	std::vector<int> vertex_hd(hd.nverts(), 0);
	FOR(e, hd.nsegments()) FOR(ev, 2) vertex_hd[hd.vert(e, ev)]++;
	FOR(v, hd.nverts()) if (vertex_hd[v] > 2) {
		points.push_back(hd.points[v]);
	}
}
void grid_pattern_extractor::set_feature_points(const std::vector<UM::vec3>& points) {
	std::vector<BBox3> boxes(grid_.ncells());
	FOR(h, grid_.ncells()) {
		FOR(hv, 8) boxes[h].add(grid_.points[grid_.vert(h, hv)]);
	}
	HBoxes looker(boxes);

	for (vec3 P : points) {
		BBox3 box; box.add(P);
		box.dilate(1E-6);
		std::vector<int> intersect;
		looker.intersect(box, intersect);
		for (int i : intersect) {
			std::array<vec3, 8> hex;
			FOR(j, 8) hex[j] = grid_.points[grid_.vert(i, j)];
			if (intersections::point_is_in_hex(hex, P)) {
				wishes_[i][0] = P;
				break;
			}
		}
	}

}

void grid_pattern_extractor::set_hard_edges(const UM::PolyLine& hd, bool detect_features_points) {
	std::cerr << "Matching hardedges with hexes..." << std::endl;
	if (detect_features_points) {
		std::vector<vec3> features_points;
		extract_feature_points(hd, features_points);
		set_feature_points(features_points);
	}
	
	std::vector<std::array<vec3, 4>> Qs;
	VolumeConnectivity vec(grid_);
	std::vector<int> map2Qs(grid_.ncells() * 6, -1);
	FOR(h, grid_.ncells()) FOR(cf, 6) if (map2Qs[6 * h + cf] == -1) {
		std::array<vec3, 4> Q;
		FOR(cfv, 4) Q[cfv] = grid_.points[grid_.facet_vert(h, cf, cfv)];
		map2Qs[6 * h + cf] = Qs.size();
		Qs.push_back(Q);
		int opp_cf2 = vec.adjacent[grid_.facet(h, cf)];
		if (opp_cf2 == -1) continue;
		map2Qs[opp_cf2] = map2Qs[6 * h + cf];
	}
	std::vector<BBox3> boxes(Qs.size());
	FOR(q, Qs.size()) {
		FOR(qv, 4) boxes[q].add(Qs[q][qv]);
	}
	HBoxes looker(boxes);
	std::vector<vec3> Q_intersect(Qs.size(), Marchinghex::NO_WISH);

	FOR(e, hd.nsegments()) {
		vec3 A = hd.points[hd.vert(e, 0)];
		vec3 B = hd.points[hd.vert(e, 1)];
		BBox3 box; box.add(A); box.add(B);
		box.dilate(1E-6);
		std::vector<int> intersect;
		looker.intersect(box, intersect);
		for (int i : intersect) {
			vec3 inter;
			if (intersections::inter_sec_quad(Qs[i][0], Qs[i][1], Qs[i][2], Qs[i][3], A, B, inter)) {
				Q_intersect[i] = inter;
			}
		}
	}
	FOR(h, grid_.ncells()) FOR(hf, 6) {
		wishes_[h][hf + 1] = Q_intersect[map2Qs[6 * h + hf]];
	}
}



inline double grid_tolerance(const Hexahedra& m) {
	double min_edge = 1E-6;
	FOR(h, m.ncells()) FOR(hv1, 8) FOR(hv2, 8) if (hv1 != hv2) {
		min_edge = std::min(min_edge, 0.1 * (m.points[m.vert(h, hv1)] - m.points[m.vert(h, hv2)]).norm());
	}
	return min_edge;
}


void grid_pattern_extractor::extract_patterns(UM::Hexahedra& hex) {
	std::cerr << "Running marchinghex..." << std::endl;
	hex.clear();
	FOR(h, grid_.ncells()) {
		if (debug_)
		if (h % 1000 == 0)std::cerr << " Extracting hexes		 :		h = " << h << "    /    " << grid_.ncells() << std::endl;

		std::array<vec3, 8> verts;
		FOR(hv, 8) verts[hv] = grid_.points[grid_.vert(h, hv)];

		std::array<int, 24> he_array;
		Marchinghex::extract_polyedra_dual(vert_is_in_[h], he_array);
		Marchinghex::Polyedra polyedra;
		Marchinghex::extract_polyedra(he_array, vert_is_in_[h], verts, bnd_intersections_[h], polyedra);
		Marchinghex::Triangulated_Polyedra polyedra_convex_hull;
		Marchinghex::compute_polyedra_convex_hull(polyedra, polyedra_convex_hull);
		Marchinghex::Hexadreized_Polyedra hexes;
		Marchinghex::extract_hexes_from_polyedra(polyedra_convex_hull, polyedra, hexes, wishes_[h]);

		int off_v = hex.points.create_points(hexes.verts_.size());
		FOR(v, hexes.verts_.size()) hex.points[off_v + v] = hexes.verts_[v];
		int off_c = hex.create_cells(hexes.hexaedra_.size());
		FOR(h, hexes.hexaedra_.size()) FOR(i, 8) hex.vert(off_c + h, i) = off_v + hexes.hexaedra_[h][i];
	}
	std::cerr << "Simpliying extra vertices..." << std::endl;
	double tol = grid_tolerance(hex);
	if (debug_) std::cerr << "tol = " << tol << std::endl;
	simplifying_vertices(hex, 1E-6);
}

inline std::array<bool, 8> inverse_cell(const std::array<bool, 8>& is_in) {
	std::array<bool, 8> res;
	FOR(i, 8) res[i] = !is_in[i];
	return res;
}

void  grid_pattern_extractor::extract(const int h, const std::array<bool, 8>& c_is_in, Marchinghex::Hexadreized_Polyedra& hexes, bool proj) {
	std::array<vec3, 8> verts;
	FOR(hv, 8) verts[hv] = grid_.points[grid_.vert(h, hv)];

	std::array<int, 24> he_array;
	Marchinghex::extract_polyedra_dual(c_is_in, he_array);
	Marchinghex::Polyedra polyedra;
	Marchinghex::Triangulated_Polyedra polyedra_convex_hull;

	if (proj) {
		Marchinghex::extract_polyedra(he_array, c_is_in, verts, bnd_intersections_[h], polyedra);
		Marchinghex::compute_polyedra_convex_hull(polyedra, polyedra_convex_hull);
		Marchinghex::extract_hexes_from_polyedra(polyedra_convex_hull, polyedra, hexes, wishes_[h]);
	}
	else {
		Marchinghex::extract_polyedra(he_array, c_is_in, verts, default_bnd_intersections_[h], polyedra);
		Marchinghex::compute_polyedra_convex_hull(polyedra, polyedra_convex_hull);
		Marchinghex::extract_hexes_from_polyedra(polyedra_convex_hull, polyedra, hexes, Marchinghex::DEFAULT_NO_WISH);
	}
}

void grid_pattern_extractor::extract_patterns(UM::Hexahedra& hex, std::vector<bool>& is_in, bool do_projection) {
	std::vector<vec3> vert_wish; std::vector<Marchinghex::VERT_TYPE> wish_type;
	extract_patterns(hex, is_in, vert_wish, wish_type);
	if (do_projection) {
		FOR(v, hex.nverts()) hex.points[v] = vert_wish[v];
	}
}


void  grid_pattern_extractor::extract_patterns(UM::Hexahedra& hex, std::vector<bool>& is_in, std::vector<vec3>& vert_wish, std::vector<Marchinghex::VERT_TYPE>& wish_type) {
	std::cerr << "Running marchinghex..." << std::endl;
	grid_pre_processing();
	hex.clear();
	int nb_cells_done = 0;
	#pragma omp parallel
	{
		std::vector<std::array<int, 8>> new_hexes;
		std::vector<bool> new_hex_is_in;

		std::vector<vec3> new_verts;
		std::vector<vec3> new_verts_wish;
		std::vector<Marchinghex::VERT_TYPE> new_verts_wish_type;
		#pragma omp for nowait
		FOR(h, grid_.ncells()) {
			//#pragma omp atomic
			//++nb_cells_done;
			//if (debug_) if (nb_cells_done % 1000 == 0) {
			//	#pragma omp critical
			//	std::cerr << " Extracting hexes		 :		c = " << (nb_cells_done / 1000) * 1000 << "    /    " << grid_.ncells() << std::endl;
			//}

			std::array<std::array<bool, 8>, 2> in_out = { vert_is_in_[h], inverse_cell(vert_is_in_[h]) };
			FOR(i, 2) {
				Marchinghex::Hexadreized_Polyedra proj_hexes, valid_hexes;
				extract(h, in_out[i], proj_hexes, true);
				extract(h, in_out[i], valid_hexes, false);

				int start_v = new_verts.size();
				FOR(v, valid_hexes.verts_.size()) new_verts.push_back(valid_hexes.verts_[v]);
				FOR(v, valid_hexes.verts_.size()) new_verts_wish.push_back(proj_hexes.verts_[v]);
				FOR(v, valid_hexes.verts_.size()) new_verts_wish_type.push_back(proj_hexes.vert_type_[v]);
				FOR(h, valid_hexes.hexaedra_.size()) {
					std::array<int, 8> hex;
					FOR(hv, 8) hex[hv] = start_v + valid_hexes.hexaedra_[h][hv];
					new_hexes.push_back(hex);
					new_hex_is_in.push_back(i == 0);

				}
			}

		}
		#pragma omp critical
		{
			int off_v = hex.points.create_points(new_verts.size());
			FOR(v, new_verts.size()) hex.points[off_v + v] = new_verts[v];
			vert_wish.resize(vert_wish.size() + new_verts.size());
			wish_type.resize(wish_type.size() + new_verts_wish_type.size());
			FOR(v, new_verts.size()) vert_wish[off_v + v] = new_verts_wish[v];
			FOR(v, new_verts.size()) wish_type[off_v + v] = new_verts_wish_type[v];

			int off_c = hex.create_cells(new_hexes.size());
			is_in.resize(is_in.size() + new_hexes.size());
			FOR(h, new_hexes.size()) FOR(i, 8) hex.vert(off_c + h, i) = off_v + new_hexes[h][i];
			FOR(h, new_hexes.size()) is_in[off_c + h]= new_hex_is_in[h];
		}
	}

	std::cerr << "Simpliying extra vertices..." << std::endl;
	double tol = grid_tolerance(grid_);
	if (debug_) std::cerr << "tol = " << tol << std::endl;

	//removing isolated verts
	std::vector<bool> tokill(hex.nverts(), true);
	FOR(h, hex.ncells()) FOR(hv, 8)  tokill[hex.vert(h, hv)] = false;
	std::vector<int> old2new;
	hex.points.delete_points(tokill, old2new);
	FOR(v, vert_wish.size()) if (old2new[v] >= 0) vert_wish[old2new[v]] = vert_wish[v];
	FOR(v, vert_wish.size()) if (old2new[v] >= 0) wish_type[old2new[v]] = wish_type[v];
	FOR(h, hex.ncells())  FOR(hc, 8) hex.vert(h, hc) = old2new[hex.vert(h, hc)];
	vert_wish.resize(hex.nverts());
	wish_type.resize(hex.nverts());

	// colocating vertices 
	std::vector<vec3> vecarray(hex.points.begin(), hex.points.end());
	colocate(vecarray, old2new, tol);
	DisjointSet ds(hex.nverts());
	FOR(v, hex.nverts()) ds.merge(v, old2new[v]);
	int nb_group = ds.get_sets_id(old2new);

	std::vector<vec3> new_verts(nb_group);
	FOR(v, hex.nverts()) new_verts[old2new[v]] = hex.points[v];
	hex.points.data->clear();
	hex.points.data->assign(new_verts.begin(), new_verts.end());

	FOR(v, vert_wish.size()) if (old2new[v] >= 0) vert_wish[old2new[v]] = vert_wish[v];
	FOR(v, vert_wish.size()) if (old2new[v] >= 0) wish_type[old2new[v]] = wish_type[v];
	FOR(h, hex.ncells())  FOR(hc, 8) hex.vert(h, hc) = old2new[hex.vert(h, hc)];
	vert_wish.resize(hex.nverts());
	wish_type.resize(hex.nverts());
}








void grid_pattern_extractor::grid_pre_processing() {
	std::vector<std::array<bool, 3>> split(grid_.ncells(), { 0 });

	constexpr int binary_facet_vertex[6][4] = { {0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, {0,1,2,3}, {4,5,6,7} };
	constexpr int rot_facet_vertex[6][4] = { {0,4,6,2}, {1,3,7,5}, {0,1,5,4}, {2,6,7,3}, {0,2,3,1}, {4,5,7,6} };

	CellAttribute<int> need_split(grid_, 0);
	std::vector<int> pile;
	// marking split needed
	FOR(h, grid_.ncells()) {
		FOR(d, 3) FOR(dd, 2) { 
			std::array<bool, 4> s;
			FOR(i, 4) s[i] = vert_is_in_[h][rot_facet_vertex[2 * d + dd][i]];
			// if face is 0-1-0-1 or 1-0-1-0 it needs splitting:
			if (s[0] && !s[1] && s[2] && !s[3] || !s[0] && s[1] && !s[2] && s[3]) split[h][d] = true;
		}
		if (split[h][0] || split[h][1] || split[h][2]) {
			pile.push_back(h);
			need_split[h] = 1;
		}
	}
	
	// propagating
	VolumeConnectivity vec(grid_);
	while (!pile.empty() ) {
		int h = pile.back(); pile.pop_back();
		FOR(d, 3) if (split[h][d]) {
			if (need_split[h] == 0) need_split[h] = 2;
			FOR(hf, 6) if ((hf / 2) != d) {

				int opp_ghf = vec.adjacent[6 * h + hf];
				if (opp_ghf == -1) continue;
				int opp_h = opp_ghf / 6;
				int opp_hf = opp_ghf % 6;
				int opp_d = opp_hf / 2;
				int opp_dd = opp_hf % 2;
				std::array<int, 4> q1, q2;
				std::array<int, 8> h1, h2;
				FOR(hfv, 4) {
					q1[hfv] = grid_.facet_vert(h, hf, hfv);
					q2[hfv] = grid_.facet_vert(opp_h, opp_hf, hfv);
				}
				FOR(hv, 8) {
					h1[hv] = grid_.vert(h, hv);
					h2[hv] = grid_.vert(opp_h, hv);
				}
				int opp_i0 = 0;
				FOR(opp_hfv, 4) if (grid_.vert(opp_h, binary_facet_vertex[opp_hf][opp_hfv]) == grid_.vert(h, binary_facet_vertex[hf][0]))
						opp_i0 = opp_hfv;


				int s = opp_i0 % 2;
				if (opp_dd == (hf % 2)) {
					s += 1;
				}

				int opp_split_d = (d + s + ((opp_d - (hf / 2)) % 3)) % 3;
				if (!split[opp_h][opp_split_d]) {
					split[opp_h][opp_split_d] = true;
					pile.push_back(opp_h);
				}
			}
		}

	}


	// splitting
	constexpr int d_jump[3][3] = {
		{1,2,4},
		{2,1,4},
		{4,1,2},
	};
	FOR(dim_cut, 3) {
		FOR(h, grid_.ncells()) if (split[h][dim_cut]) {
			split[h][dim_cut] = false;
			std::array<vec3, 8> hex;
			FOR(hv, 8) hex[hv] = grid_.points[grid_.vert(h, hv)];

			std::array<vec3, 8> hex1, hex2;
			std::array<bool, 8> config1, config2; 
			std::array<vec3, 24> cba1, cba2, dcba1, dcba2;

			FOR(i, 4) {
				int o_i_1 = binary_facet_vertex[2 * dim_cut][i];
				int o_i_2 = binary_facet_vertex[2 * dim_cut + 1][i];
				hex1[o_i_1] = hex[o_i_1];
				config1[o_i_1] = vert_is_in_[h][o_i_1];
				hex2[o_i_2] = hex[o_i_2];
				config2[o_i_2] = vert_is_in_[h][o_i_2];
				FOR(d, 3) cba1[3 * o_i_1 + d] = bnd_intersections_[h][3 * o_i_1 + d];
				FOR(d, 3) cba2[3 * o_i_2 + d] = bnd_intersections_[h][3 * o_i_2 + d];
			}

			FOR(d2, 2) FOR(d1, 2) {

				int i = d2 * d_jump[dim_cut][1] + d1 * d_jump[dim_cut][2];
				int fi = d1 * 2 + d2;
				int o_i_1 = binary_facet_vertex[2 * dim_cut][fi];
				int o_i_2 = binary_facet_vertex[2 * dim_cut + 1][fi];
				vec3 new_v = 0.5 * (hex[o_i_1] + hex[o_i_2]);
				if (vert_is_in_[h][o_i_2] && !vert_is_in_[h][o_i_1]) new_v = 0.5 * (hex[o_i_2] + bnd_intersections_[h][3 * o_i_2 + dim_cut]);
				if (!vert_is_in_[h][o_i_2] && vert_is_in_[h][o_i_1]) new_v = 0.5 * (hex[o_i_1] + bnd_intersections_[h][3 * o_i_1 + dim_cut]);
				hex1[d_jump[dim_cut][0] + i] = new_v;
				hex2[i] = new_v;
				bool has_in = vert_is_in_[h][o_i_2] || vert_is_in_[h][o_i_1];
				config1[d_jump[dim_cut][0] + i] = has_in;
				config2[i] = has_in;

				FOR(d, 3) cba1[3 * (d_jump[dim_cut][0] + i) + d] = vec3(1e100, 1e100, 1e100);
				FOR(d, 3) cba2[3 * (i)+d] = vec3(1e100, 1e100, 1e100);

				cba1[3 * (d_jump[dim_cut][0] + i) + dim_cut] = bnd_intersections_[h][3 * o_i_2 + dim_cut];
				cba2[3 * (i)+dim_cut] = bnd_intersections_[h][3 * o_i_1 + dim_cut];
			}
			FOR(z, 2) FOR(y, 2) FOR(x, 2) FOR(d, 3) {

				int s[3] = { 0,0,0 }; s[d] = 1;
				dcba1[3 * bin_to_int(x, y, z) + d] = 0.5 * (hex1[bin_to_int(x, y, z)]
					+ hex1[bin_to_int(x + s[0], y + s[1], z + s[2])]);
				dcba2[3 * bin_to_int(x, y, z) + d] = 0.5 * (hex2[bin_to_int(x, y, z)]
					+ hex2[bin_to_int(x + s[0], y + s[1], z + s[2])]);
				if (cba1[3 * bin_to_int(x, y, z) + d][0] == 1e100)
					cba1[3 * bin_to_int(x, y, z) + d] = 0.5 * (hex1[bin_to_int(x, y, z)]
						+ hex1[bin_to_int(x + s[0], y + s[1], z + s[2])]);
				if (cba2[3 * bin_to_int(x, y, z) + d][0] == 1e100)
					cba2[3 * bin_to_int(x, y, z) + d] = 0.5 * (hex2[bin_to_int(x, y, z)]
						+ hex2[bin_to_int(x + s[0], y + s[1], z + s[2])]);
			}

			int off_set_vert = grid_.points.create_points(16);
			FOR(hv, 8) {
				grid_.points[off_set_vert + hv] = hex1[hv];
				grid_.points[off_set_vert + 8 + hv] = hex2[hv];
			}
			int off_set_cells = grid_.create_cells(1);
			FOR(hv, 8) {
				grid_.vert(h, hv) = off_set_vert + hv;
				grid_.vert(off_set_cells, hv) = off_set_vert + 8 + hv;
			}
			vert_is_in_[h] = config1;
			vert_is_in_.push_back(config2);
			bnd_intersections_[h] = cba1;
			bnd_intersections_.push_back(cba2);
			wishes_[h] = Marchinghex::DEFAULT_NO_WISH;
			wishes_.push_back(Marchinghex::DEFAULT_NO_WISH);
			default_bnd_intersections_[h] = dcba1;
			default_bnd_intersections_.push_back(dcba2);
			split.push_back(split[h]);
		}
	}

	std::vector<bool> tokill(grid_.nverts(), true);
	FOR(h, grid_.ncells()) FOR(hv, 8)  tokill[grid_.vert(h, hv)] = false;
	std::vector<int> old2new;
	grid_.points.delete_points(tokill, old2new);
	FOR(h, grid_.ncells())  FOR(hc, 8) grid_.vert(h, hc) = old2new[grid_.vert(h, hc)];

}

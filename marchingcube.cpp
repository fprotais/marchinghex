#include "marchingcube.h"
#include "intersections.h"
#include <map>
#include <chrono>
#include <ultimaille/disjointset.h>
#include <ultimaille/colocate.h>

#define FOR(i, n) for(int i = 0; i < n; i++)
static constexpr int NOT_AN_ID = -1;
using namespace UM;
namespace marchingcube {
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
		if (Pfacet < 6) {
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

		std::vector<std::vector<int>> P_to_TP_facets(P.facets_.size());
		FOR(t, TPol.triangles_.size()) P_to_TP_facets[TPol.triangles_original_facets_[t]].push_back(t);
		std::vector<vec3> absolute_center(nb_of_composant);
		bool has_feature_corner = false;
		std::vector<int> nb_of_center(nb_of_composant, 0);
		FOR(f, P.facets_.size()) {
			if (P.facets_[f].size() < 3) continue;
			vec3 center = center_of_polyedra_facet_with_flattening(TPol, f);

			if (f > 5 && is_wish(wishes[0])) {
				if (nb_of_composant < 2) {
					center = wishes[0];
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
			std::map<std::pair<int, int>, vec3> edge_mid_vertices_wish;

			FOR(f, P.facets_.size()) FOR(fv, P.facets_[f].size()) {
				const int v1 = P.facets_[f][fv];
				const int v2 = P.facets_[f][(fv + 1) % P.facets_[f].size()];
				if (f < 6 && P.vertex_id_[v1] > 7 && P.vertex_id_[v2] > 7 && is_wish(wishes[f+1])) {
					edge_mid_vertices_wish[ordered_pair(v1, v2)] = wishes[f + 1];
				}
				edge_mid_vertices[ordered_pair(v1, v2)] = NOT_AN_ID;
			}
			for (auto& edge_mid_vertex : edge_mid_vertices) {
				edge_mid_vertex.second = (int)hexes.verts_.size() + (int)new_verts.size();
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
		}
		if (!has_feature_corner) FOR(f, P.facets_.size()) if (f > 5) {
			std::vector<vec3> set_vertices;
			FOR(fv, P.facets_[f].size()) {
				std::pair<int, int> p = ordered_pair(P.facets_[f][fv], P.facets_[f][(fv + 1) % P.facets_[f].size()]);
				if (vertex_is_hardedge[edge_mid_vertices[p]])
					set_vertices.push_back(hexes.verts_[edge_mid_vertices[p]]);
			}
			if (set_vertices.size() > 1) {
				vec3 new_center;
				FOR(i, set_vertices.size()) new_center += set_vertices[i];
				new_center = new_center / set_vertices.size();
				hexes.verts_[nb_pol_verts + nb_of_composant + f] = new_center;

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

	static void extract_nodes(const std::vector<vec3>& geometry, const std::vector<std::array<int, 4>>& tetraedras, std::vector<std::vector<std::vector<bool>>>& Vert_is_in, const int start_in_X, const int start_in_Y, const int start_in_Z) {
		FOR(tet, tetraedras.size()) {
			if (tet % 1000 == 0) std::cerr << "Extracting sources       :      " << tet << "  /   " << tetraedras.size() << std::endl;
			assert(tetraedras[tet][0] < geometry.size());
			assert(tetraedras[tet][1] < geometry.size());
			assert(tetraedras[tet][2] < geometry.size());
			assert(tetraedras[tet][3] < geometry.size());
			vec3 A = geometry[tetraedras[tet][0]];
			vec3 B = geometry[tetraedras[tet][1]];
			vec3 C = geometry[tetraedras[tet][2]];
			vec3 D = geometry[tetraedras[tet][3]];
			vec3 min(1E10, 1E10, 1E10), max(-1E10, -1E10, -1E10);
			FOR(cv, 4) FOR(d, 3) {
				if (geometry[tetraedras[tet][cv]][d] > max[d]) max[d] = geometry[tetraedras[tet][cv]][d];
				if (geometry[tetraedras[tet][cv]][d] < min[d]) min[d] = geometry[tetraedras[tet][cv]][d];
			}
			std::array<double, 4> l;
			for (int X = int(std::floor(min[0])); X <= std::ceil(max[0]); X++)
				for (int Y = int(std::floor(min[1])); Y <= std::ceil(max[1]); Y++)
					for (int Z = int(std::floor(min[2])); Z <= std::ceil(max[2]); Z++)
						if (intersections::point_is_in_tet(A, B, C, D, vec3(X, Y, Z), l)) {
							Vert_is_in[X - start_in_X][Y - start_in_Y][Z - start_in_Z] = true;
						}
		}
		
	}
	static void extracting_bnd_voxel_intersection(const std::vector<vec3>& geometry, const std::vector<std::array<int, 3>>& tris, std::vector<std::vector<std::vector<std::array<vec3, 3>>>>& Bnd_intersection,
		int start_in_X, int start_in_Y, int start_in_Z, int nb_in_X, int nb_in_Y, int nb_in_Z) {
		std::cerr << "matching boundary" << std::endl;

		for (int X = start_in_X; X < start_in_X + (int)nb_in_X; X++) {
			for (int Y = start_in_Y; Y < start_in_Y + (int)nb_in_Y; Y++) {
				for (int Z = start_in_Z; Z < start_in_Z + (int)nb_in_Z; Z++) {
					Bnd_intersection[X - start_in_X][Y - start_in_Y][Z - start_in_Z][0] = vec3(X + 0.5, Y, Z);
					Bnd_intersection[X - start_in_X][Y - start_in_Y][Z - start_in_Z][1] = vec3(X, Y + 0.5, Z);
					Bnd_intersection[X - start_in_X][Y - start_in_Y][Z - start_in_Z][2] = vec3(X, Y, Z + 0.5);
				}
			}
		}
		FOR(tri, tris.size()) {
			if (tri % 1000 == 0) std::cerr << "Matching boundary triangles with voxels :      " << tri << "  /   " << tris.size() << std::endl;
			assert(tris[tri][0] < geometry.size());
			assert(tris[tri][1] < geometry.size());
			assert(tris[tri][2] < geometry.size());

			vec3 A = geometry[tris[tri][0]];
			vec3 B = geometry[tris[tri][1]];
			vec3 C = geometry[tris[tri][2]];
			vec3 min(1E10, 1E10, 1E10), max(-1E10, -1E10, -1E10);
			FOR(tv, 3) FOR(d, 3) {
				if (geometry[tris[tri][tv]][d] > max[d]) max[d] = geometry[tris[tri][tv]][d];
				if (geometry[tris[tri][tv]][d] < min[d]) min[d] = geometry[tris[tri][tv]][d];
			}
			double l1 = 0, l2 = 0, l3 = 0;
			for (int X = int(std::floor(min[0])); X <= std::ceil(max[0]); X++)
				for (int Y = int(std::floor(min[1])); Y <= std::ceil(max[1]); Y++)
					for (int Z = int(std::floor(min[2])); Z <= std::ceil(max[2]); Z++) {
						vec3 Xi = vec3(X, Y, Z);
						FOR(d, 3) {
							vec3 Xj = Xi; Xj[d] += 1;
							if (intersections::inter_sec_triangle(A, B, C, Xi, Xj, l1, l2, l3)) {
								vec3 bndcut = l1 * A + l2 * B + l3 * C;
								if ((bndcut - Xi).norm() < 0.1) bndcut = 0.9 * Xi + 0.1 * Xj;
								if ((bndcut - Xj).norm() < 0.1) bndcut = 0.1 * Xi + 0.9 * Xj;
								Bnd_intersection[X - start_in_X][Y - start_in_Y][Z - start_in_Z][d] = bndcut;
							}

						}
					}

		}

	}



	static void extracting_voxel_hardedge_wishes(const std::vector<vec3>& geometry, const std::vector<std::array<int, 2>>& edges, std::vector<std::vector<std::vector<std::array<vec3, 4>>>>& voxel_wishes,
		int start_in_X, int start_in_Y, int start_in_Z, int nb_in_X, int nb_in_Y, int nb_in_Z) {
		for (int X = start_in_X; X < start_in_X + (int)nb_in_X; X++) {
			for (int Y = start_in_Y; Y < start_in_Y + (int)nb_in_Y; Y++) {
				for (int Z = start_in_Z; Z < start_in_Z + (int)nb_in_Z; Z++) {
					FOR(i, 4) voxel_wishes[X - start_in_X][Y - start_in_Y][Z - start_in_Z][i] = NO_WISH;
				}
			}
		}
		FOR(edge, edges.size()) FOR(ev, 2) assert(edges[edge][ev] < geometry.size());
		std::vector<int> vertex_hd(geometry.size(), 0);
		FOR(edge, edges.size()) FOR(ev, 2) vertex_hd[edges[edge][ev]]++;
		FOR(v, geometry.size()) if (vertex_hd[v] > 2) {
			voxel_wishes[(int)std::floor(geometry[v].x) - start_in_X][(int)std::floor(geometry[v].y) - start_in_Y][(int)std::floor(geometry[v].z) - start_in_Z][0] = geometry[v];
		}
		FOR(edge, edges.size()) {
			if (edge % 1000 == 0) std::cerr << "Intersecting hardedges with voxels :      " << edge << "  /   " << edges.size() << std::endl;
			
			vec3 A = geometry[edges[edge][0]];
			vec3 B = geometry[edges[edge][1]];
			vec3 min(1E10, 1E10, 1E10), max(-1E10, -1E10, -1E10);
			FOR(tv, 2) FOR(d, 3) {
				if (geometry[edges[edge][tv]][d] > max[d]) max[d] = geometry[edges[edge][tv]][d];
				if (geometry[edges[edge][tv]][d] < min[d]) min[d] = geometry[edges[edge][tv]][d];
			}
			std::vector<std::array<int, 3>> list_of_voxel_crossed;
			double l1 = 0, l2 = 0, l3 = 0;
			for (int X = int(std::floor(min[0])); X <= std::ceil(max[0]); X++)
				for (int Y = int(std::floor(min[1])); Y <= std::ceil(max[1]); Y++)
					for (int Z = int(std::floor(min[2])); Z <= std::ceil(max[2]); Z++) {
						vec3 Xi = vec3(X, Y, Z);

						FOR(dim, 3) {
							vec3 d1; d1[(dim + 1) % 3] = 1;
							vec3 d2; d2[(dim + 2) % 3] = 1;
							vec3 inter;
							if (intersections::inter_sec_quad(Xi, Xi + d1, Xi + d1 + d2, Xi + d2, A, B, inter)) {
								voxel_wishes[X - start_in_X][Y - start_in_Y][Z - start_in_Z][dim + 1] = inter;
							}
						}

					}

		}

	}
	static void simplifying_vertices(std::vector<vec3>& mc_verts, std::vector<int>& mc_cells) {
		//removing isolated verts
		{
			std::vector<bool> mark(mc_verts.size(), false);
			FOR(hv, mc_cells.size())  mark[mc_cells[hv]] = true;
			std::vector<vec3> new_verts;
			std::vector<int> old2new(mc_verts.size());
			int curr = 0;
			FOR(v, mc_verts.size()) if (mark[v]) {
				new_verts.push_back(mc_verts[v]);
				old2new[v] = curr++;
			}
			mc_verts.clear();
			mc_verts.assign(new_verts.begin(), new_verts.end());
			FOR(hv, mc_cells.size())  mc_cells[hv] = old2new[mc_cells[hv]];

		}

		// colocating vertices 
		{
			std::vector<int> old2new;
			colocate(mc_verts, old2new, 1E-6);
			DisjointSet ds(mc_verts.size());
			FOR(v, mc_verts.size()) ds.merge(v, old2new[v]);
			int nb_group = ds.get_sets_id(old2new);
			std::cerr << nb_group << std::endl;
			std::cerr << mc_verts.size() << std::endl;

			std::vector<vec3> new_verts(nb_group);
			FOR(v, mc_verts.size()) new_verts[old2new[v]] = mc_verts[v];
			mc_verts.clear();
			mc_verts.assign(new_verts.begin(), new_verts.end());
			FOR(hv, mc_cells.size())  mc_cells[hv] = old2new[mc_cells[hv]];
		}
	}


	void hexify(const std::vector<vec3>& verts, const std::vector<int>& tets, std::vector<vec3>& out_verts, std::vector<int>& out_hexes) {
		if (verts.empty() || tets.empty()) return;
		std::vector<std::array<int, 4>> tetraedras(tets.size() / 4);
		std::vector<vec3> geometry(verts.size());
		FOR(v, geometry.size()) geometry[v] = verts[v];
		FOR(f, tetraedras.size()) FOR(fv, 4) tetraedras[f][fv] = tets[4 * f + fv];


		vec3 bboxmax = geometry[0];
		vec3 bboxmin = geometry[0];

		FOR(v, geometry.size()) FOR(dim, 3) {
			if (bboxmax[dim] < geometry[v][dim]) bboxmax[dim] = geometry[v][dim];
			if (bboxmin[dim] > geometry[v][dim]) bboxmin[dim] = geometry[v][dim];
		}
		std::cerr << " bboxmax  : " << bboxmax << std::endl;
		std::cerr << " bboxmin  : " << bboxmin << std::endl;
		int nb_in_X = (int)std::ceil(bboxmax[0] - bboxmin[0]) + 1;
		int nb_in_Y = (int)std::ceil(bboxmax[1] - bboxmin[1]) + 1;
		int nb_in_Z = (int)std::ceil(bboxmax[2] - bboxmin[2]) + 1;
		int start_in_X = (int)std::floor(bboxmin[0]);
		int start_in_Y = (int)std::floor(bboxmin[1]);
		int start_in_Z = (int)std::floor(bboxmin[2]);
		std::vector<std::vector<std::vector<bool>>> Vert_is_in(nb_in_X, std::vector<std::vector<bool>>(nb_in_Y, std::vector<bool>(nb_in_Z, false)));
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		extract_nodes(geometry, tetraedras, Vert_is_in, start_in_X, start_in_Y, start_in_Z);
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cerr << "Source extraction: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000. << "sec." << std::endl;
		std::vector<vec3> mc_verts;
		std::vector<int> mc_cells;
		begin = std::chrono::steady_clock::now();
		for (int X = start_in_X; X < start_in_X + nb_in_X - 1; X++) {
			std::cerr << " Extracting cubes		 :		X = " << X << "    ->    " << start_in_X + nb_in_X - 1 << std::endl;
			for (int Y = start_in_Y; Y < start_in_Y + nb_in_Y - 1; Y++) {
				for (int Z = start_in_Z; Z < start_in_Z + nb_in_Z - 1; Z++)
				{
					std::array<bool, 8> corner_is_in;
					std::array<vec3, 8> corners;
					std::array<vec3, 24> corner_bnd_alternative;
					bool one_is_in = false;
					FOR(k, 2) FOR(j, 2) FOR(i, 2) {
						corner_is_in[4 * k + 2 * j + i] = Vert_is_in[X + i - start_in_X][Y + j - start_in_Y][Z + k - start_in_Z];
						one_is_in = one_is_in || corner_is_in[4 * k + 2 * j + i];
					}
					if (!one_is_in) continue;
					FOR(k, 2) FOR(j, 2) FOR(i, 2) {
						corners[4 * k + 2 * j + i] = vec3((double)X + (double)i, (double)Y + (double)j, (double)Z + (double)k);
					}
					FOR(k, 2) FOR(j, 2) FOR(i, 2) FOR(dim, 3) {
						int v = 4 * k + 2 * j + i;
						vec3 shift; shift[dim] = 0.5;
						if (dim == 0 && i == 1) shift[0] *= -1;
						if (dim == 1 && j == 1) shift[1] *= -1;
						if (dim == 2 && k == 1) shift[2] *= -1;
						corner_bnd_alternative[3 * v + dim] = corners[v] + shift;
					}

					std::array<int, 24> he_array;
					extract_polyedra_dual(corner_is_in, he_array);
					Polyedra polyedra;
					extract_polyedra(he_array, corner_is_in, corners, corner_bnd_alternative, polyedra);
					Triangulated_Polyedra polyedra_convex_hull;
					compute_polyedra_convex_hull(polyedra, polyedra_convex_hull);
					Hexadreized_Polyedra hexes;
					extract_hexes_from_polyedra(polyedra_convex_hull, polyedra, hexes);

					int vert_init_size = (int)mc_verts.size();
					mc_verts.resize(vert_init_size + hexes.verts_.size());
					FOR(v, hexes.verts_.size()) mc_verts[vert_init_size + v] = hexes.verts_[v];
					int hex_init_size = (int)mc_cells.size();
					mc_cells.resize(hex_init_size + 8 * hexes.hexaedra_.size());
					FOR(h, hexes.hexaedra_.size()) FOR(i, 8) mc_cells[hex_init_size + 8 * h + i] = vert_init_size + hexes.hexaedra_[h][i];
				}
			}
		}
		end = std::chrono::steady_clock::now();
		std::cerr << "Hexaedra extraction: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000. << "sec." << std::endl;
		std::cerr << "Simplifying unecessary vertices" << std::endl;
		begin = std::chrono::steady_clock::now();
		simplifying_vertices(mc_verts, mc_cells);
		end = std::chrono::steady_clock::now();
		std::cerr << "Colocating and vertex simplification: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000. << " sec." << std::endl;
		out_verts.assign(mc_verts.begin(), mc_verts.end());
		out_hexes.assign(mc_cells.begin(), mc_cells.end());
	}




	void hexify_with_bnd(const std::vector<UM::vec3>& verts, const std::vector<int>& tets, const std::vector<int> bnd_tri, const std::vector<int> hard_edge, std::vector<vec3>& out_verts, std::vector<int>& out_hexes){
		if (verts.empty() || tets.empty()) return;
		std::vector<std::array<int, 4>> tetraedras(tets.size() / 4);
		std::vector<std::array<int, 3>> tris(bnd_tri.size() / 3);
		std::vector<std::array<int, 2>> edges(hard_edge.size() / 2);
		std::vector<vec3> geometry(verts.size());
		FOR(v, geometry.size()) geometry[v] = verts[v];
		FOR(f, tetraedras.size()) FOR(fv, 4) tetraedras[f][fv] = tets[4 * f + fv];
		FOR(t, tris.size()) FOR(tv, 3) tris[t][tv] = bnd_tri[3 * t + tv];
		FOR(e, edges.size()) FOR(ev, 2) edges[e][ev] = hard_edge[2 * e + ev];


		vec3 bboxmax = geometry[0];
		vec3 bboxmin = geometry[0];

		FOR(v, geometry.size()) FOR(dim, 3) {
			if (bboxmax[dim] < geometry[v][dim]) bboxmax[dim] = geometry[v][dim];
			if (bboxmin[dim] > geometry[v][dim]) bboxmin[dim] = geometry[v][dim];
		}
		std::cerr << " bboxmax  : " << bboxmax << std::endl;
		std::cerr << " bboxmin  : " << bboxmin << std::endl;
		int nb_in_X = (int)std::ceil(bboxmax[0] - bboxmin[0]) + 2;
		int nb_in_Y = (int)std::ceil(bboxmax[1] - bboxmin[1]) + 2;
		int nb_in_Z = (int)std::ceil(bboxmax[2] - bboxmin[2]) + 2;
		int start_in_X = (int)std::floor(bboxmin[0]);
		int start_in_Y = (int)std::floor(bboxmin[1]);
		int start_in_Z = (int)std::floor(bboxmin[2]);
		std::cerr << " start_in_X  : " << start_in_X << std::endl;
		std::cerr << " start_in_Y  : " << start_in_Y << std::endl;
		std::cerr << " start_in_Z  : " << start_in_Z << std::endl;
		std::cerr << " nb_in_X  : " << nb_in_X << std::endl;
		std::cerr << " nb_in_Y  : " << nb_in_Y << std::endl;
		std::cerr << " nb_in_Z  : " << nb_in_Z << std::endl;

		std::vector<std::vector<std::vector<bool>>> Vert_is_in(nb_in_X, std::vector<std::vector<bool>>(nb_in_Y, std::vector<bool>(nb_in_Z, false)));
		std::vector<std::vector<std::vector<std::array<vec3, 3>>>> Bnd_intersection(nb_in_X, std::vector<std::vector<std::array<vec3, 3>>>(nb_in_Y, std::vector<std::array<vec3, 3>>(nb_in_Z)));
		std::vector<std::vector<std::vector<std::array<vec3, 4>>>> voxel_wishes(nb_in_X, std::vector<std::vector<std::array<vec3, 4>>>(nb_in_Y, std::vector<std::array<vec3, 4>>(nb_in_Z)));

		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		extract_nodes(geometry, tetraedras, Vert_is_in, start_in_X, start_in_Y, start_in_Z);
		extracting_bnd_voxel_intersection(geometry, tris, Bnd_intersection, start_in_X, start_in_Y, start_in_Z, nb_in_X, nb_in_Y, nb_in_Z);
		extracting_voxel_hardedge_wishes(geometry, edges, voxel_wishes, start_in_X, start_in_Y, start_in_Z, nb_in_X, nb_in_Y, nb_in_Z);
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cerr << "Data extraction: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000. << "sec." << std::endl;
		
		std::vector<vec3> mc_verts;
		std::vector<int> mc_cells;
		begin = std::chrono::steady_clock::now();
		for (int X = start_in_X; X < start_in_X + nb_in_X - 1; X++) {
			std::cerr << " Extracting cubes		 :		X = " << X << "    ->    " << start_in_X + nb_in_X - 1 << std::endl;
			for (int Y = start_in_Y; Y < start_in_Y + nb_in_Y - 1; Y++) {
				for (int Z = start_in_Z; Z < start_in_Z + nb_in_Z - 1; Z++)
				{
					std::array<bool, 8> corner_is_in;
					std::array<vec3, 8> corners;
					std::array<vec3, 24> corner_bnd_alternative;
					bool one_is_in = false;
					FOR(k, 2) FOR(j, 2) FOR(i, 2) {
						corner_is_in[4 * k + 2 * j + i] = Vert_is_in[X + i - start_in_X][Y + j - start_in_Y][Z + k - start_in_Z];
						one_is_in = one_is_in || corner_is_in[4 * k + 2 * j + i];
					}
					if (!one_is_in) continue;
					FOR(k, 2) FOR(j, 2) FOR(i, 2) {
						corners[4 * k + 2 * j + i] = vec3((double)X + (double)i, (double)Y + (double)j, (double)Z + (double)k);
					}
					FOR(k, 2) FOR(j, 2) FOR(i, 2) FOR(dim, 3) {
						int ip[3] = { i,j,k };
						if (ip[dim] == 1) ip[dim] = 0;
						corner_bnd_alternative[3 * (4 * k + 2 * j + i) + dim] = Bnd_intersection[X + ip[0] - start_in_X][Y + ip[1] - start_in_Y][Z + ip[2] - start_in_Z][dim];
					}

					std::array<vec3, 7> center_and_facet_wishes = { NO_WISH };
					center_and_facet_wishes[0] = voxel_wishes[X - start_in_X][Y - start_in_Y][Z - start_in_Z][0];
					center_and_facet_wishes[1] = voxel_wishes[X - start_in_X][Y - start_in_Y][Z - start_in_Z][1];
					center_and_facet_wishes[3] = voxel_wishes[X - start_in_X][Y - start_in_Y][Z - start_in_Z][2];
					center_and_facet_wishes[5] = voxel_wishes[X - start_in_X][Y - start_in_Y][Z - start_in_Z][3];
					center_and_facet_wishes[2] = voxel_wishes[X + 1 - start_in_X][Y - start_in_Y][Z - start_in_Z][1];
					center_and_facet_wishes[4] = voxel_wishes[X - start_in_X][Y + 1 - start_in_Y][Z - start_in_Z][2];
					center_and_facet_wishes[6] = voxel_wishes[X - start_in_X][Y - start_in_Y][Z + 1 - start_in_Z][3];

					std::array<int, 24> he_array;
					extract_polyedra_dual(corner_is_in, he_array);
					Polyedra polyedra;
					extract_polyedra(he_array, corner_is_in, corners, corner_bnd_alternative, polyedra);
					Triangulated_Polyedra polyedra_convex_hull;
					compute_polyedra_convex_hull(polyedra, polyedra_convex_hull);
					Hexadreized_Polyedra hexes;
					extract_hexes_from_polyedra(polyedra_convex_hull, polyedra, hexes, center_and_facet_wishes);

					int vert_init_size = (int)mc_verts.size();
					mc_verts.resize(vert_init_size + hexes.verts_.size());
					FOR(v, hexes.verts_.size()) mc_verts[vert_init_size + v] = hexes.verts_[v];
					int hex_init_size = (int)mc_cells.size();
					mc_cells.resize(hex_init_size + 8 * hexes.hexaedra_.size());
					FOR(h, hexes.hexaedra_.size()) FOR(i, 8) mc_cells[hex_init_size + 8 * h + i] = vert_init_size + hexes.hexaedra_[h][i];
				}
			}
		}
		end = std::chrono::steady_clock::now();
		std::cerr << "Hexaedra extraction: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000. << "sec." << std::endl;
		std::cerr << "Simplifying unecessary vertices" << std::endl;
		begin = std::chrono::steady_clock::now();
		simplifying_vertices(mc_verts, mc_cells);
		end = std::chrono::steady_clock::now();
		std::cerr << "Colocating and vertex simplification: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000. << " sec." << std::endl;
		out_verts.assign(mc_verts.begin(), mc_verts.end());
		out_hexes.assign(mc_cells.begin(), mc_cells.end());
		
	}

}


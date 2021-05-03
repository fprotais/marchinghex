#include "quality.h"
#include <limits>

#include "projection.h"

using namespace UM;


double hausdorff_dist(const std::vector<vec3>& P, const std::vector<std::array<vec3, 3>>& mesh, std::vector<double>& dist2mesh) {
	if (P.size() == 0 || mesh.size() == 0) return std::numeric_limits<double>::infinity();
	dist2mesh.clear();
	dist2mesh.resize(P.size(), 0);
	double max_dist = 0;

	std::vector<BBox3> boxes(mesh.size());
	for (int t : range(mesh.size())) for (int tv : range(3)) boxes[t].add(mesh[t][tv]);
	HBoxes looker(boxes);
	double avg_edge = 0;
	for (int t : range(mesh.size())) for (int tc : range(3)) avg_edge += (mesh[t][tc] - mesh[t][(tc + 1) % 3]).norm();
	avg_edge /= 3 * mesh.size();
	for (int v : range(P.size())) {
		double dilatation = 4 * avg_edge;
		BBox3 box;
		box.add(P[v]);
		std::vector<int> tris;
		while (tris.empty()) {
			box.dilate(dilatation);
			looker.intersect(box, tris);
			dilatation *= 2;
		}
		double loc_dist = std::numeric_limits<double>::max();
		for (int t : tris) {
			std::array<double, 3> l;

			vec3 X;
			double dist = point_triangle_squared_distance(P[v], mesh[t], X, l);
			loc_dist = std::min(dist, loc_dist);
		}
		dist2mesh[v] = loc_dist;
		max_dist = std::max(loc_dist, max_dist);


	}
	
	return max_dist;
}

inline void sampletriangle(const std::array<vec3, 3>& tri, int edge_sampling, std::vector<vec3>& points) {
	if (edge_sampling == 0) return;
	if (edge_sampling == 1) {
		points.push_back((tri[0] + tri[1] + tri[2]) / 3);
		return;
	}
	for (int i : range(edge_sampling))
		for (int j : range(edge_sampling - i)) {
			double l1 = (double)i / (edge_sampling - 1);
			double l2 = (double)j / (edge_sampling - 1);
			points.push_back(l1 * tri[1] + l2 * tri[2] + (1 - l1 - l2) * tri[0]);
	}
}

double hausdorff_dist(const UM::Hexahedra& hex, const UM::Tetrahedra& tet, int tri_sampling) {
	if (hex.nfacets() == 0 || tet.nfacets() == 0) return std::numeric_limits<double>::infinity();
	std::vector<vec3> tet_points;
	std::vector<std::array<vec3, 3>> tet_boundary;
	VolumeConnectivity tet_vec(tet);
	for (int t : range(tet.ncells())) for (int tf : range(4)) if (tet_vec.adjacent[4 * t + tf] == -1) {
		std::array<vec3, 3> tri;
		for (int tfv : range(3)) tri[tfv] = tet.points[tet.facet_vert(t, tf, tfv)];
		tet_boundary.push_back(tri);
		sampletriangle(tri, tri_sampling, tet_points);
	}

	std::vector<vec3> hex_points;
	std::vector<std::array<vec3, 3>> hex_boundary;
	VolumeConnectivity hex_vec(hex);
	for (int h : range(hex.ncells())) for (int hf : range(6)) if (hex_vec.adjacent[6 * h + hf] == -1) {
		std::array<vec3, 4> q;
		for (int hfv : range(4)) q[hfv] = hex.points[hex.facet_vert(h, hf, hfv)];
		vec3 center = 0.25 * (q[0] + q[1] + q[2] + q[3]);
		for (int hfv : range(4)){
			hex_boundary.push_back({ q[hfv], q[(hfv + 1) % 4], center });
			sampletriangle(hex_boundary.back(), tri_sampling, hex_points);
		}
	}
	std::vector<double> hex_dists_to_tet;
	std::vector<double> tet_dists_to_hex;
	double tet2hex = hausdorff_dist(tet_points, hex_boundary, tet_dists_to_hex);
	double hex2tet = hausdorff_dist(hex_points, tet_boundary, hex_dists_to_tet);

	{ // to dist distances:
		Triangles mesh_hex_bnd; mesh_hex_bnd.points.data->assign(hex_points.begin(), hex_points.end());
		Triangles mesh_tet_bnd; mesh_tet_bnd.points.data->assign(tet_points.begin(), tet_points.end());
		PointAttribute<double> att_hex_points(mesh_hex_bnd);
		for (int v : range(mesh_hex_bnd.nverts())) att_hex_points[v] = hex_dists_to_tet[v];
		PointAttribute<double> att_tet_points(mesh_tet_bnd);
		for (int v : range(mesh_tet_bnd.nverts())) att_tet_points[v] = tet_dists_to_hex[v];

		write_by_extension("dist_hex_mesh.mesh", hex);
		write_by_extension("dist_tet_mesh.mesh", tet);
		write_by_extension("dist_hex_sampling.geogram", mesh_hex_bnd, SurfaceAttributes{ {{"dist", att_hex_points.ptr}},{},{} });
		write_by_extension("dist_tet_sampling.geogram", mesh_tet_bnd, SurfaceAttributes{ {{"dist", att_tet_points.ptr}},{},{} });

	}



	return std::max(tet2hex, hex2tet);
}
constexpr int HEX_CORNER_SPLITTING[8][4] = {
	{0,1,2,4},
	{1,3,0,5},
	{2,0,3,6},
	{3,2,1,7},
	{4,6,5,0},
	{5,4,7,1},
	{6,7,4,2},
	{7,5,6,3},
};
double compute_scaled_jacobian(const UM::Hexahedra& hex, UM::CellAttribute<double>& cell_min_sj) {
	double glob_min = 1;
	for (int h : range(hex.ncells())) {
		double min_sj = 1;
		for (int hv : range(8)) {
			std::array<vec3, 4> v;
			for (int i : range(4)) v[i] = hex.points[hex.vert(h, HEX_CORNER_SPLITTING[hv][i])];
			vec3 n1 = v[1] - v[0]; n1.normalize();
			vec3 n2 = v[2] - v[0]; n2.normalize();
			vec3 n3 = v[3] - v[0]; n3.normalize();
			min_sj = std::min(min_sj, n3 * cross(n1, n2));
		}
		cell_min_sj[h] = min_sj;
		glob_min = std::min(glob_min, min_sj);
	}
	return glob_min;
}

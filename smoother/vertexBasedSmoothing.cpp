#include "vertexBasedSmoothing.h"
#include "projection.h"

#include <set>

using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)

#define DEBUG 0
#if DEBUG
#define DEBUGCODE(X) X 
#else
#define DEBUGCODE(X) 
#endif
#define DBGMSG(X) DEBUGCODE(std::cout << X << "\n";)
#define DBGVARIABLE(X) DEBUGCODE(std::cout << #X << ": " << X << "\n";)



vertex_smoother::vertex_smoother(UM::Hexahedra& m)
: _m(m)

, _pts(*m.points.data)
, _vert_data(_pts.size())
, _vert_order(_pts.size())

, _hex2tets(_m.ncells())
, _tets(_m.ncells() * 8)

{
	const vec3 hex_ref[8] = {
		{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
		{0,0,1}, {1,0,1}, {0,1,1}, {1,1,1},
	};
	std::array<UM::vec3, 4> tet_ref[8];
	FOR(hc, 8) {
		vec3 ref_pt[4]; 
		FOR(i,4) ref_pt[i]= hex_ref[HEX_TO_SJ_TET_SPLIT[hc][i]];
		UM::mat3x3 M = { ref_pt[1] - ref_pt[0], ref_pt[2] - ref_pt[0], ref_pt[3] - ref_pt[0] };
		UM::mat3x3 invM = M.invert();
		double detM = M.det();
		invM = invM.transpose();
		tet_ref[hc] = { -invM[0] - invM[1] - invM[2], invM[0], invM[1], invM[2] };
	}
	
	std::vector<std::set<int>> vert2vert(_pts.size());
	FOR(h, _m.ncells()) FOR(hc, 8) {
		int t = 8*h+hc;
		_hex2tets[h][hc] = t;
		FOR(tv, 4) _tets[t].verts[tv] = _m.vert(h, HEX_TO_SJ_TET_SPLIT[hc][tv]);

		_tets[t].pre_computed = tet_ref[hc];

		FOR(tv, 4) FOR(tv2, 3) vert2vert[_tets[t].verts[tv]].insert(_tets[t].verts[(tv + tv2 + 1) % 4]);
		FOR(tv, 4) _vert_data[_tets[t].verts[tv]].tets.push_back(t);
	}

	FOR(v, _pts.size()) {
		_vert_order[v] = v;
		_vert_data[v].verts.assign(vert2vert[v].begin(), vert2vert[v].end());
	}

	_proj.vert2triangles.resize(_pts.size());
 	_proj.vert2segments.resize(_pts.size());
 	_proj.vert2point.resize(_pts.size());
}

double vertex_smoother::compute_tet_SJ(int t) const {
	vec3 v0 = _pts[_tets[t].verts[1]] - _pts[_tets[t].verts[0]];
	vec3 v1 = _pts[_tets[t].verts[2]] - _pts[_tets[t].verts[0]];
	vec3 v2 = _pts[_tets[t].verts[3]] - _pts[_tets[t].verts[0]];
	double det = v2 * cross(v0, v1);
	double n0 = v0.norm();
	double n1 = v1.norm();
	double n2 = v2.norm();
	double dem = n0 * n1 * n2;
	if (dem < 1e-12) return 0;
	return det / dem;
}

double vertex_smoother::compute_hex_minSJ(int h) const {
	double minSJ = 1.;
	for (int t : _hex2tets[h]) {
		minSJ = std::min(minSJ, compute_tet_SJ(t));
	}
	return minSJ;
}

double vertex_smoother::compute_vert_minSJ(int v) const {
	double minSJ = 1.;
	for (int t : _vert_data[v].tets) {
		minSJ = std::min(minSJ, compute_tet_SJ(t));
	}
	return minSJ;
}

void vertex_smoother::set_locked_vertices(const std::vector<bool>& locks) {
	assert(_pts.size() == locks.size());
	FOR(v, _pts.size()) _vert_data[v].locked = locks[v];
}
void vertex_smoother::set_bnd_triangles(const UM::Triangles& bnd) {
	_proj.bnd.points.data->assign(bnd.points.data->begin(), bnd.points.data->end());
	_proj.bnd.facets.assign(bnd.facets.begin(), bnd.facets.end());
}
void vertex_smoother::set_features_segment(const UM::PolyLine& ft) {
	_proj.ft.points.data->assign(ft.points.data->begin(), ft.points.data->end());
	_proj.ft.segments.assign(ft.segments.begin(), ft.segments.end());
}
void vertex_smoother::set_vertex_triangles(const int i, const std::vector<int>& triangles) {
	_vert_data[i].type = 1;
	_proj.vert2triangles[i].assign(triangles.begin(), triangles.end());

}
void vertex_smoother::set_vertex_segments(const int i, const std::vector<int>& segments) {
	_vert_data[i].type = 2;
	_proj.vert2segments[i].assign(segments.begin(), segments.end());
}
void vertex_smoother::set_vertex_point(const int i, const UM::vec3 v) {
	_vert_data[i].type = 3;
	_proj.vert2point[i] = _proj.pts.size(); _proj.pts.push_back(v);
}

void vertex_smoother::update_order() {
	int cnt = 0;
	for (int t = 3; t >= 0; t--) {
		for (int i = 0; i < _pts.size(); i++) {
			if (_vert_data[i].type == t) _vert_order[cnt++] = i;
		}
	}
}

void vertex_smoother::udpdate_tet_quality(int t) {
	_tets[t].SJ = compute_tet_SJ(t);
	if (_tets[t].SJ <= _min_SJ_update && !_tets[t].is_bad) {
		_tets[t].is_bad = true;
		for (int v : _tets[t].verts) _vert_data[v].nb_bad_tets++;
	} 
	else if (_tets[t].SJ > _min_SJ_update && _tets[t].is_bad) {
		_tets[t].is_bad = false;
		for (int v : _tets[t].verts) _vert_data[v].nb_bad_tets--;
	}
}
void vertex_smoother::update_bad_elements() {
	FOR(t, _tets.size()) udpdate_tet_quality(t);
	FOR(v, _pts.size()) if (_vert_data[v].type == 0) _vert_data[v].is_bnd_compliant = true;
}

UM::vec3 vertex_smoother::projected_position(int v) const {
	vec3 proj = _pts[v];
	switch (_vert_data[v].type) {
	case 1:
		{
		double closet_dist = 1E100;
		for (int f : _proj.vert2triangles[v]) {
			std::array<vec3, 3> tri;
			for (int fv : range(3)) tri[fv] = _proj.bnd.points[_proj.bnd.vert(f, fv)];
			std::array<double, 3> l;
			vec3 X;
			double d = point_triangle_squared_distance(_pts[v], tri, X, l);
			if (d < closet_dist) {
				proj = X;
				closet_dist = d;
			}
		}
		}
		break;
	case 2:
		{
		double closet_dist = 1E100;
		for (int s : _proj.vert2segments[v]) {
			std::array<vec3, 2> seg;
			for (int sv : range(2)) seg[sv] = _proj.ft.points[_proj.ft.vert(s, sv)];
			std::array<double, 2> l;
			vec3 X;
			double d = point_segment_squared_distance(_pts[v], seg, X, l);
			if (d < closet_dist) {
				proj = X;
				closet_dist = d;
			}
		}
		}
		break;
	case 3:
		{
		proj = _proj.pts[_proj.vert2point[v]];
		}
		break;
	}
	return proj;
}


void vertex_smoother::run_iter() {
	auto maxSJ_line_search = [&](int v, vec3 dir) {
		constexpr double Tau[] = { 4,2,1,.5,.25, 0.125 };
		vec3 originalPt = _pts[v];
		vec3 bestLocation = _pts[v];
		double maxMinSJ = compute_vert_minSJ(v);
		DBGVARIABLE(maxMinSJ);
		for (double tau : Tau) {
			_pts[v] = originalPt + tau * dir; 
			double sj = compute_vert_minSJ(v);
			if (sj > maxMinSJ) {
				DBGMSG("new choice: " << tau << "| sj: " << sj);
				maxMinSJ = sj;
				bestLocation = _pts[v];
			}
		}
		return bestLocation;
	};
	auto furthest_line_search = [&](int v, vec3 dir) {
		constexpr double Tau[] = { 1,0.875,.75, 0.625, .5,.25, 0.125, 0.05 };
		vec3 originalPt = _pts[v];
		double minSJ = std::min(_min_SJ_project_, compute_vert_minSJ(v));
		DBGVARIABLE(minSJ);
		for (double tau : Tau) {
			_pts[v] = originalPt + tau * dir; 
			double sj = compute_vert_minSJ(v);
			if (sj > minSJ) {
				DBGMSG("new choice: " << tau << "| sj: " << sj);
				return _pts[v];
			}
		}
		return originalPt;
	};
	for (int v : _vert_order) {
		if (_vert_data[v].locked) continue;
		if (_vert_data[v].is_bnd_compliant && _vert_data[v].nb_bad_tets == 0) continue;
		DBGVARIABLE(v);
		DBGVARIABLE(_pts[v]);
		mat3x3 hess;
		vec3 dir = compute_vert_elliptic_grad_hess(v, hess);
		DBGVARIABLE(hess.det());
		DBGVARIABLE(hess);
		if (hess.det() >= 1e-10) {
			dir = compute_naive_laplacian_direction(v);
		}
		else {
			DBGMSG("lapl");
			dir = hess.invert() * dir;
			dir *= -1;
		}
		DBGVARIABLE(dir);
		_pts[v] = maxSJ_line_search(v, dir);
		DBGVARIABLE(_pts[v]);
		if (_vert_data[v].type > 0) {
			vec3 projPos = projected_position(v);
			DBGVARIABLE(projPos);
			_pts[v] = furthest_line_search(v, projPos - _pts[v]);
			if ((_pts[v] - projPos).norm() < 1e-8) _vert_data[v].is_bnd_compliant = true;
			DBGVARIABLE(_pts[v]);
			DBGVARIABLE(_vert_data[v].is_bnd_compliant);
		}
		for (int t : _vert_data[v].tets) udpdate_tet_quality(t);
	}
}
UM::vec3 vertex_smoother::compute_naive_laplacian_direction(int v) const {
	vec3 newPos = _pts[v];
	for(int o_v : _vert_data[v].verts) newPos += _pts[v];
	newPos /= _vert_data[v].verts.size() + 1;
	return newPos - _pts[v];
}


// UTILS :

inline UM::mat3x3 dual_basis(const UM::mat3x3& J) {
	return
	{
		{{
			 J[1].y * J[2].z - J[1].z * J[2].y,
			 J[1].z * J[2].x - J[1].x * J[2].z,
			 J[1].x * J[2].y - J[1].y * J[2].x
		 },
		{
			J[0].z * J[2].y - J[0].y * J[2].z,
			J[0].x * J[2].z - J[0].z * J[2].x,
			J[0].y * J[2].x - J[0].x * J[2].y
		},
		{
			J[0].y * J[1].z - J[0].z * J[1].y,
			J[0].z * J[1].x - J[0].x * J[1].z,
			J[0].x * J[1].y - J[0].y * J[1].x
		}}
	};
}

inline double chi(double eps, double det) {
	if (det > 0)
		return (det + std::sqrt(eps * eps + det * det)) * .5;
	return .5 * eps * eps / (std::sqrt(eps * eps + det * det) - det);
}

inline double chi_deriv(double eps, double det) {
	return .5 + det / (2. * std::sqrt(eps * eps + det * det));
}


// END UTILS

vec3 vertex_smoother::compute_vert_elliptic_grad_hess(int v, UM::mat3x3& hess) const {
	vec3 grad = { 0,0,0 };
	hess = { UM::vec3(0,0,0), UM::vec3(0,0,0),UM::vec3(0,0,0) };
	double E = 0;
	for (int t : _vert_data[v].tets) {

		mat3x3 J = { UM::vec3(0,0,0), UM::vec3(0,0,0),UM::vec3(0,0,0) };
		FOR(tv, 4) FOR(d, 3)
			J[d] += _tets[t].pre_computed[tv] * _pts[_tets[t].verts[tv]][d];
		mat3x3 K = dual_basis(J);
		double detJ = J.det();

		double c1 = chi(_eps, detJ);
		double c2 = std::pow(c1, 2. / 3.);
		double c3 = chi_deriv(_eps, detJ);
		double f = (J[0] * J[0] + J[1] * J[1] + J[2] * J[2]) / c2;
		double g = (1 + detJ * detJ) / c1;
		E += ((1 - _theta) * f + _theta * g);

		int loc_id = 0;
		FOR(tv, 4) if ( _tets[t].verts[tv] == v) loc_id = tv;
		FOR(d, 3) {
			UM::vec3 dfda = J[d] * (2. / c2) - K[d] * ((2. * f * c3) / (3. * c1));
			UM::vec3 dgda = K[d] * ((2 * detJ - g * c3) / c1);
			grad[d] += (dfda + (dgda - dfda) * _theta) * _tets[t].pre_computed[loc_id];
		}

		double ch1 = (1 - _theta) / c2;
		double ch2 = (1 - _theta) * (-4.) / 3. * c3 / std::pow(c1, 5. / 3.);

		double ch4 = (1 - _theta) * (2. / 3.) * (1 + 2. / 3.) * (f / (c1 * c1)) 
			+ _theta / c1 * (2 - 4 * detJ * c3 / c1 + 2 * (1 + detJ * detJ) * c3 * c3 / (c1 * c1));
		std::array<double, 9> a, b;
		FOR(d1, 3) FOR(d2, 3) {
			a[3 * d1 + d2] = J[d1][d2];
			b[3 * d1 + d2] = K[d1][d2];
		}
		std::array<double, 9> vch;
		FOR(k, 9) vch[k] = a[k] * ch2 + b[k] * ch4;
		mat<9, 9> M;
		FOR(k, 9) FOR(l, 9) {
			if (k == l) M[k][l] = ch1;
			else M[k][l] = 0;
			M[k][l] += b[l] * vch[k] + a[l] * b[k] * ch2;
		}

		FOR(d1, 3) {
			FOR(d2, 3) {
				mat3x3 locH;
				FOR(j, 3) FOR(k, 3) for (int k : range(3)) locH[j][k] = M[3 * d1 + j][3 * d2 + k];
				hess[d1][d2] += _tets[t].pre_computed[loc_id] * (locH * _tets[t].pre_computed[loc_id]);
			}
		}
	}
	return grad;
}

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

#define PICVALUEMSG(x, val, VAR) if (x == val) std::cout << #x << "= " << val << " -> " << #VAR << ": " << VAR << "\n";

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
		_vert_data[v].minSJ = compute_vert_minSJ(v);
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

void vertex_smoother::update_tet_quality(int t) {
	_tets[t].SJ = compute_tet_SJ(t);
	double eq = compute_tet_elliptic_quality(t);
	bool bad = _tets[t].SJ <= _min_SJ_target
			 || eq > 4; // magic value? 
	if (bad && !_tets[t].is_bad) {
		_tets[t].is_bad = true;
		for (int v : _tets[t].verts) _vert_data[v].nb_bad_tets++;
		return;
	} 
	if (!bad && _tets[t].is_bad) {
		_tets[t].is_bad = false;
		for (int v : _tets[t].verts) _vert_data[v].nb_bad_tets--;
	}
}

void vertex_smoother::update_bad_elements() {
	FOR(t, _tets.size()) update_tet_quality(t);
	FOR(v, _pts.size())  _vert_data[v].is_bnd_compliant = (_vert_data[v].type == 0);
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
		constexpr double Tau[] = { 1,.5,.25, 0.125, 0.0625 };
		vec3 originalPt = _pts[v];
		vec3 bestLocation = _pts[v];
		double best_E = compute_vert_elliptic_energy(v);
		double minSJ = std::min(_vert_data[v].minSJ, _min_SJ_update);
		DBGVARIABLE(best_E);
		for (double tau : Tau) {
			_pts[v] = originalPt + tau * dir; 
			double e = compute_vert_elliptic_energy(v);
			if (e < best_E) {
				double sj = compute_vert_minSJ(v);
				if (sj < minSJ) continue;
				DBGMSG("new choice: " << tau << "| sj: " << e);
				best_E = e;
				_vert_data[v].minSJ = sj;
				bestLocation = _pts[v];
			}
		}
		return bestLocation;
	};
	auto furthest_line_search = [&](int v, vec3 dir) {
		constexpr double Tau[] = { 1,0.875,.75, 0.625, .5,.25, 0.125, 0.05 };
		vec3 originalPt = _pts[v];
		double minSJ = std::min(_min_SJ_project_, compute_vert_minSJ(v));
		double maxeq = std::max(compute_vert_elliptic_max_quality(v), 10.);

		DBGVARIABLE(minSJ);
		for (double tau : Tau) {
			_pts[v] = originalPt + tau * dir; 
			double sj = compute_vert_minSJ(v);
			double eq = compute_vert_elliptic_max_quality(v);
			if (sj > minSJ && eq < maxeq) {
				DBGMSG("new choice: " << tau << "| sj: " << sj);
				_vert_data[v].minSJ = sj;
				return ;
			}
		}
		_pts[v] = originalPt;
		return;
	};

	for (int v : _vert_order) {
		if (_vert_data[v].locked) continue;
		if (_vert_data[v].is_bnd_compliant && _vert_data[v].nb_bad_tets == 0) continue;
		DBGVARIABLE(v);

		if (_vert_data[v].is_bnd_compliant) {
			DBGVARIABLE(_pts[v]);
			mat3x3 hess;
			vec3 dir = compute_vert_elliptic_grad_trunc_hess(v, hess);
			DBGVARIABLE(hess.det());
			DBGVARIABLE(hess);
			dir = hess.invert() * dir;
			dir *= -1;

			if (dir.norm() <= 1e-10 || dir.norm() >= 1e6) {
				dir = compute_naive_laplacian_direction(v);
				DBGMSG("lapl");
			}
			DBGVARIABLE(dir);
			_pts[v] = maxSJ_line_search(v, dir);
			DBGVARIABLE(_pts[v]);
		}
		if (_vert_data[v].type > 0) {
			DBGVARIABLE(_pts[v]);
			vec3 projPos = projected_position(v);
			DBGVARIABLE(projPos);
			furthest_line_search(v, projPos - _pts[v]);
			_vert_data[v].is_bnd_compliant = (_pts[v] - projPos).norm() < 1e-8;
			DBGVARIABLE(_pts[v]);
			DBGVARIABLE(_vert_data[v].is_bnd_compliant);
		}
		for (int t : _vert_data[v].tets) update_tet_quality(t);
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

double vertex_smoother::compute_tet_elliptic_quality(int t) const{
	mat3x3 J = { UM::vec3(0,0,0), UM::vec3(0,0,0),UM::vec3(0,0,0) };
	FOR(tv, 4) FOR(d, 3)
		J[d] += _tets[t].pre_computed[tv] * _pts[_tets[t].verts[tv]][d];
	double detJ = J.det();
	if (detJ < 1e-30) return 1e10;
	double c = std::pow(detJ, 2. / 3.);
	double f =  (J[0] * J[0] + J[1] * J[1] + J[2] * J[2]) / c;
	return std::min(1e10, f);
}
	
double vertex_smoother::compute_vert_elliptic_max_quality(int v) const {
	double maxEQ = 1.;
	for (int t : _vert_data[v].tets) {
		maxEQ = std::max(maxEQ, compute_tet_elliptic_quality(t));
	}
	return maxEQ;
}


double vertex_smoother::compute_tet_elliptic_energy(int t) const {
	mat3x3 J = { UM::vec3(0,0,0), UM::vec3(0,0,0),UM::vec3(0,0,0) };
	FOR(tv, 4) FOR(d, 3)
		J[d] += _tets[t].pre_computed[tv] * _pts[_tets[t].verts[tv]][d];
	double detJ = J.det();
	double c1 = chi(_eps, detJ);
	double c2 = std::pow(c1, 2. / 3.);
	return (J[0] * J[0] + J[1] * J[1] + J[2] * J[2]) / c2;
}

double vertex_smoother::compute_vert_elliptic_energy(int v) const {
	double E = 0;
	for (int t : _vert_data[v].tets) {
		E += compute_tet_elliptic_energy(t);
	}
	return E;
}

vec3 vertex_smoother::compute_vert_elliptic_grad_trunc_hess(int v, UM::mat3x3& hess) const {
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
		E +=  f;

		int loc_id = 0;
		FOR(tv, 4) if ( _tets[t].verts[tv] == v) loc_id = tv;

		double val1 =  2. / c2;
		double val2 =  (2. * f * c3) / (3. * c1);

		double c3dc1 = c3 / c1;
		double val4 = val1 * (2./3.) * c3dc1;//(4.*c3)/(3.*c2*c1);
		double val5 = val2 * (5./3.) * c3dc1;//(10.*f* c3 * c3)/(9.*c1 * c1);
		
		mat<3,3> Id = mat<3,3>::identity();

		FOR(d, 3) { 
			UM::vec3 dfda = J[d] * (val1) - K[d] * (val2);
			grad[d] += dfda * _tets[t].pre_computed[loc_id];
		
			//would looping over d1 and d2 here would work? is it more/less operations?
			mat<3,1> A = {{{J[d].x}, {J[d].y}, {J[d].z}}};
			mat<3,1> B = {{{K[d].x}, {K[d].y}, {K[d].z}}};

			mat<3,3> AB = A*B.transpose();
			mat<3,3> BB = B*B.transpose();

            mat<3,3> locH = Id*(val1) - (AB + AB.transpose())*(val4) + (BB)*(val5);
			hess[d][d] += _tets[t].pre_computed[loc_id] * (locH * _tets[t].pre_computed[loc_id]);

		}

	}
	return grad;
}

#include "locale_iterative_smoother.h"
#include "projection.h"

using namespace UM;

// UTILS :

const vec3 hex_ref[8] = {
	{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
	{0,0,1}, {1,0,1}, {0,1,1}, {1,1,1},
};

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

iterative_smoother::iterative_smoother(Hexahedra& m)
	: m_(m)
	, vert_one_ring_(m_.nverts())
	, vert_eps_(m_.nverts(), 0.01)
	, lock_(m_.nverts(), false)
	, vert_min_det_(m_.nverts(), 1.)
	, vert_bnd_projection_triangles_(m_.nverts())
	, vert_bnd_projection_segments_(m_.nverts())
	, vert_bnd_projection_point_(m_.nverts())
	, vert_type_(m_.nverts(), 0)
	, order_(m_.nverts(), 0)
{
	scale_ = 0;
	for (int c : range(m_.ncells())) for (int cf : range(6)) for (int cfv : range(4))
		scale_ += (m.points[m.facet_vert(c, cf, cfv)] - m.points[m.facet_vert(c, cf, (cfv + 1) % 4)]).norm();
	scale_ /= m.ncells() * 24;
	for (int v : range(m_.nverts())) m.points[v] /= scale_;
	for (int v : range(m_.nverts())) order_[v] = v;
}


void iterative_smoother::set_locks(const std::vector<bool>& locks) {
	assert(locks.size() == lock_.size());
	lock_.assign(locks.begin(), locks.end());
}
void iterative_smoother::set_bnd_triangles(const UM::Triangles& bnd) {
	bnd_.points.data->assign(bnd.points.data->begin(), bnd.points.data->end());
	bnd_.facets.assign(bnd.facets.begin(), bnd.facets.end());
	for (int v : range(bnd_.nverts())) bnd_.points[v] /= scale_;
}
void iterative_smoother::set_features_segment(const UM::PolyLine& ft) {
	ft_.points.data->assign(ft.points.data->begin(), ft.points.data->end());
	ft_.segments.assign(ft.segments.begin(), ft.segments.end());
	for (int v : range(ft_.nverts())) ft_.points[v] /= scale_;
}
void iterative_smoother::set_vertex_triangles(const int i, const std::vector<int>& triangles) {
	vert_type_[i] = 1;
	vert_bnd_projection_triangles_[i].assign(triangles.begin(), triangles.end());

}
void iterative_smoother::set_vertex_segments(const int i, const std::vector<int>& segments) {
	vert_type_[i] = 2;
	vert_bnd_projection_segments_[i].assign(segments.begin(), segments.end());
}
void iterative_smoother::set_vertex_point(const int i, const UM::vec3 v) {
	vert_type_[i] = 3;
	vert_bnd_projection_point_[i] = pts_.size(); pts_.push_back(v / scale_);
}
void iterative_smoother::scale_back() {
	for (int v : UM::range(m_.nverts())) m_.points[v] *= scale_;
}

void iterative_smoother::scale_up() {
	for (int v : UM::range(m_.nverts())) m_.points[v] /= scale_;
}

double iterative_smoother::compute_hextet_jacobian(const int i, const vec3& v, const int h, const int t, UM::mat3x3& J, UM::mat3x3& K) {
	J = { UM::vec3(0,0,0), UM::vec3(0,0,0),UM::vec3(0,0,0) };
	std::array<vec3, 4> V;
	for (int j : range(4)) {
		int index = m_.vert(h, optimized_tets_[t][j]);
		if (index == i) V[j] = v;
		else V[j] = m_.points[index];
	}
	for(int j : range(4)) for (int d : range(3))
		J[d] += hexrefcoef_[t][j] * V[j][d];
	K = dual_basis(J);
	return J.det();

}
double iterative_smoother::evaluate_vert_energy(const int i, const UM::vec3& v, double& minDetJ) {
	double E = 0;
	minDetJ = 1E100;
	for (int ht : vert_one_ring_[i]) {
		int h = ht / optimized_tets_.size();
		int t = ht % optimized_tets_.size();
		mat3x3 J, K;
		double detJ = compute_hextet_jacobian(i, v, h, t, J, K);
		minDetJ = std::min(minDetJ, detJ);
		double c = chi(vert_eps_[i], detJ);
		double f = (J[0] * J[0] + J[1] * J[1] + J[2] * J[2]) / std::pow(c, 2. / 3.);
		double g = (1 + detJ * detJ) / c;
		E += ((1 - theta_) * f + theta_ * g);
	}
	return E;
}


double iterative_smoother::compute_vert_grad_hess(const int i, UM::vec3& grad, UM::mat3x3& hess) {
	grad = { 0,0,0 };
	hess = { UM::vec3(0,0,0), UM::vec3(0,0,0),UM::vec3(0,0,0) };
	double E = 0;
	vert_min_det_[i] = 1E100;

	for (int ht : vert_one_ring_[i]) {
		int h = ht / optimized_tets_.size();
		int t = ht % optimized_tets_.size();
		mat3x3 J, K;
		double detJ = compute_hextet_jacobian(i, m_.points[i], h, t, J, K);
		vert_min_det_[i] = std::min(vert_min_det_[i], detJ);
		double c1 = chi(vert_eps_[i], detJ);
		double c2 = std::pow(c1, 2. / 3.);
		double c3 = chi_deriv(vert_eps_[i], detJ);
		double f = (J[0] * J[0] + J[1] * J[1] + J[2] * J[2]) / c2;
		double g = (1 + detJ * detJ) / c1;
		int i_ref = 0;
		E += ((1 - theta_) * f + theta_ * g);
		for (int j : range(4)) if (m_.vert(h, optimized_tets_[t][j]) == i) i_ref = j;
		for (int d : range(3)) {
			UM::vec3 dfda = J[d] * (2. / c2) - K[d] * ((2. * f * c3) / (3. * c1));
			UM::vec3 dgda = K[d] * ((2 * detJ - g * c3) / c1);
			grad[d] += (dfda * (1. - theta_) + dgda * theta_) * hexrefcoef_[t][i_ref];
		}
		double ch1 = (1 - theta_) / c2;
		double ch2 = (1 - theta_) * (-4.) / 3. * c3 / std::pow(c1, 5. / 3.);

		double ch4 = (1 - theta_) * (2. / 3.) * (1 + 2. / 3.) * (f / (c1 * c1)) 
			+ theta_ / c1 * (2 - 4 * detJ * c3 / c1 + 2 * (1 + detJ * detJ) * c3 * c3 / (c1 * c1));
		std::array<double, 9> a, b;
		for (int d1 : range(3)) for (int d2 : range(3)) {
			a[3 * d1 + d2] = J[d1][d2];
			b[3 * d1 + d2] = K[d1][d2];
		}
		std::array<double, 9> vch;
		for (int k : range(9)) vch[k] = a[k] * ch2 + b[k] * ch4;
		mat<9, 9> M;
		for (int k : range(9)) for (int l : range(9)) {
			if (k == l) M[k][l] = ch1;
			else M[k][l] = 0;
			M[k][l] += b[l] * vch[k] + a[l] * b[k] * ch2;
		}


		for (int d1 : range(3)) {
			for (int d2 : range(3)) {
				mat3x3 locH;
				for (int j : range(3)) for (int k : range(3)) locH[j][k] = M[3 * d1 + j][3 * d2 + k];
				hess[d1][d2] += hexrefcoef_[t][i_ref] * (locH * hexrefcoef_[t][i_ref]);
			}
		}
			 
	}
	return E;
}

double iterative_smoother::linesearch(const int i, const vec3& dir, const double prev_E, double& new_E) {
	double Tau[] = { 4,2,1,.5,.25 };
	new_E = prev_E;
	double tau_min = 0;
	for (double& tau : Tau) {
		double mindet;
		double E = evaluate_vert_energy(i, m_.points[i] - tau * dir, mindet);
		if (E < new_E) {
			new_E = E;
			tau_min = tau;
		}
	}
	return tau_min;
}

vec3 iterative_smoother::place_on_bnd(const int i, const vec3& v) {
	if (vert_type_[i] == 0) return v;
	//std::cerr << vert_type_[i] << std::endl;
	vec3 close_point;
	if (vert_type_[i] == 1) {
		if (vert_bnd_projection_triangles_[i].empty()) return v;
		double closet_dist = 1E100;
		for (int f : vert_bnd_projection_triangles_[i]) {
			std::array<vec3, 3> tri;
			for (int fv : range(3)) tri[fv] = bnd_.points[bnd_.vert(f, fv)];
			std::array<double, 3> l;
			vec3 X;
			double d = point_triangle_squared_distance(v, tri, X, l);
			if (d < closet_dist) {
				close_point = X;
				closet_dist = d;
			}
		}
	}
	if (vert_type_[i] == 2) {
		if (vert_bnd_projection_segments_[i].empty()) return v;
		double closet_dist = 1E100;
		for (int s : vert_bnd_projection_segments_[i]) {
			std::array<vec3, 2> seg;
			for (int sv : range(2)) seg[sv] = ft_.points[ft_.vert(s, sv)];
			std::array<double, 2> l;
			vec3 X;
			double d = point_segment_squared_distance(v, seg, X, l);
			if (d < closet_dist) {
				close_point = X;
				closet_dist = d;
			}
		}
	}
	if (vert_type_[i] == 3) {
		close_point = pts_[vert_bnd_projection_point_[i]];
	}

	vec3 dir = close_point - v;

	double Tau[] = { 1, 0.75, 0.5, 0.25, 0.1, 0.05, 0.01 };
	for (double& tau : Tau) {
		double minDet;
		double E = evaluate_vert_energy(i, v + tau * dir, minDet);
		if (minDet > min_det_project_) {
			vert_min_det_[i] = minDet;
			return v + tau * dir;
		}
	}

	return v;

}


void iterative_smoother::update_eps(const int i, const double prev_E, const double new_E) {
	double sigma = std::max(1. - new_E / prev_E, 1e-1);
	double mu = (1 - sigma) * chi(vert_eps_[i], vert_min_det_[i]);
	if (vert_min_det_[i] < mu)
		vert_eps_[i] = 2 * std::sqrt(mu * (mu - vert_min_det_[i]));
	else vert_eps_[i] = 0;//1e-10;

}

void iterative_smoother::update_order() {
	order_;
	int cnt = 0;
	for (int t = 3; t >= 0; t--) {
		for (int i = 0; i < m_.nverts(); i++) {
			if (vert_type_[i] == t) order_[cnt++] = i;
		}
	}
}

void iterative_smoother::update_connectivity() {
	hexrefcoef_.resize(optimized_tets_.size());
	for (int i : range(optimized_tets_.size())) {
		vec3 tet_ref[4];
		for (int j : range(4)) tet_ref[j] = hex_ref[optimized_tets_[i][j]];

		UM::mat3x3 M = { tet_ref[1] - tet_ref[0], tet_ref[2] - tet_ref[0], tet_ref[3] - tet_ref[0] };
		UM::mat3x3 invM = M.invert();
		double detM = M.det();
		invM = invM.transpose();
		hexrefcoef_[i] = { -invM[0] - invM[1] - invM[2], invM[0], invM[1], invM[2] };
	}
	vert_one_ring_.clear(); vert_one_ring_.resize(m_.nverts());
	for (int c : range(m_.ncells())) {
		for (int t : range(optimized_tets_.size())) for (int tv : range(4)) {
			vert_one_ring_[m_.vert(c, optimized_tets_[t][tv])].push_back(optimized_tets_.size() * c + t);
		}
	}
}


void iterative_smoother::run_iter() {
	update_order();
	update_connectivity();
	for (int ip = 0; ip < m_.nverts(); ip++) {
		int i = order_[ip];
		if (lock_[i]) continue;
		vec3 grad;
		mat3x3 hess;
		double actual_E = compute_vert_grad_hess(i, grad, hess);

		vec3 dir = hess.invert() * grad;
		double new_E;
		double l = linesearch(i, dir, actual_E, new_E);
		m_.points[i] = place_on_bnd(i, m_.points[i] - l * dir);

		//debug
		//{
		//	PointAttribute<int> prio(m_);
		//	PointAttribute<int> order(m_);
		//	for (int v : range(m_.nverts())) prio[v] = vert_type_[v];
		//	for (int v : range(m_.nverts())) order[v] = order_[v];
		//	write_by_extension("debug.geogram", m_, VolumeAttributes{ {{"type", prio.ptr},{"order", order.ptr}},{},{},{} });

		//}

		update_eps(i, actual_E, new_E);
	}
	double glob_min = 1E100;
	for (int i = 0; i < m_.nverts(); i++) if (!lock_[i]) if (glob_min > vert_min_det_[i]) glob_min = vert_min_det_[i];
	std::cerr << "Min det is: " << glob_min << std::endl;
}





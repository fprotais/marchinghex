#include <ultimaille/all.h>
#include <array>


class boundary_matcher {
public:
	inline boundary_matcher(const UM::Triangles& bnd, const UM::PolyLine& features)
		: bnd_(bnd)
		, ft_(features)
		, bnd_looker_()
		, ft_looker_()
	{
		{
			std::vector<UM::BBox3> boxes(bnd.nfacets());
			for (int t : UM::range(bnd_.nfacets())) for (int tv : UM::range(3)) boxes[t].add(bnd_.points[bnd_.vert(t, tv)]);
			bnd_looker_.init(boxes);
		}
		{
			std::vector<UM::BBox3> boxes(ft_.nsegments());
			for (int e : UM::range(ft_.nsegments())) for (int ev : UM::range(2)) boxes[e].add(ft_.points[ft_.vert(e, ev)]);
			ft_looker_.init(boxes);
		}
	}
	inline void get_vert_close_triangles(const UM::vec3 v, std::vector<int>& triangles, const double size_of_bbox = 1) {
		UM::BBox3 box;
		box.add(v);
		box.dilate(size_of_bbox / 2);
		bnd_looker_.intersect(box, triangles);
	}
	inline void get_vert_close_segments(const UM::vec3 v, std::vector<int>& edges, const double size_of_bbox = 1) {
		UM::BBox3 box;
		box.add(v);
		box.dilate(size_of_bbox / 2);
		ft_looker_.intersect(box, edges);
	}
private:
	const UM::Triangles& bnd_;
	const UM::PolyLine& ft_;
	UM::HBoxes bnd_looker_;
	UM::HBoxes ft_looker_;
};


class iterative_smoother {	
public:
	iterative_smoother(UM::Hexahedra& m);
	void run_iter();
	void set_locks(const std::vector<bool>& locks);
	void set_bnd_triangles(const UM::Triangles& bnd);
	void set_features_segment(const UM::PolyLine& ft);
	void set_vertex_triangles(const int i, const std::vector<int>& triangles);
	void set_vertex_segments(const int i, const std::vector<int>& segments);
	void set_vertex_point(const int i, const UM::vec3 v);
	void scale_back();
	void scale_up();

	UM::Hexahedra& m_;
	UM::Triangles bnd_;
	UM::PolyLine ft_;
	std::vector<UM::vec3> pts_;
	double min_det_project_ = 0.1;
	double theta_ = 0.5;



private:
	double compute_hextet_jacobian(const int i, const UM::vec3& v, const int h, const int t, UM::mat3x3& J, UM::mat3x3& K);
	double evaluate_vert_energy(const int i, const UM::vec3& v, double& minDetJ);
	double compute_vert_grad_hess(const int i, UM::vec3& grad, UM::mat3x3& hess);
	double linesearch(const int i, const UM::vec3& dir, const double prev_E, double& new_E);
	UM::vec3 place_on_bnd(const int i, const UM::vec3& v);

	void update_eps(const int i, const double prev_E, const double new_E);
	std::vector<bool> lock_;
	std::vector<double> vert_eps_;
	std::vector<double> vert_min_det_;
	std::vector<int> vert_type_;
	std::vector<std::vector<int>> vert_one_ring_;
	std::vector<std::vector<int>> vert_bnd_projection_triangles_;
	std::vector<std::vector<int>> vert_bnd_projection_segments_;
	std::vector<int> vert_bnd_projection_point_;
	std::array<std::array<UM::vec3, 4>, 24> hexrefcoef_;
	double scale_;


};
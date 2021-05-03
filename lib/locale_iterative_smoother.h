#include <ultimaille/all.h>
#include <array>

static const std::vector<std::array<int, 4>> TET_MATRIX_HEX_LEXICOGRAPHIC_24_SPLIT = { {
	{0,1,3,4},{0,1,2,5},{0,2,6,1},{0,2,4,3},
	{0,4,1,6},{0,4,5,2},{1,3,0,7},{1,3,2,5},
	{1,5,3,4},{1,5,7,0},{2,3,7,0},{2,3,6,1},
	{2,6,0,7},{2,6,4,3},{3,7,6,1},{3,7,2,5},
	{4,5,0,7},{4,5,1,6},{4,6,5,2},{4,6,7,0},
	{5,7,1,6},{5,7,3,4},{6,7,5,2},{6,7,4,3}
} };

static const std::vector<std::array<int, 4>> TET_MATRIX_HEX_LEXICOGRAPHIC_SJ_SPLIT = { {
	{0,1,2,4},
	{1,3,0,5},
	{2,0,3,6},
	{3,2,1,7},
	{4,6,5,0},
	{5,4,7,1},
	{6,7,4,2},
	{7,5,6,3},
} };

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

	std::vector<std::array<int, 4>> optimized_tets_ = TET_MATRIX_HEX_LEXICOGRAPHIC_SJ_SPLIT;

private:
	double compute_hextet_jacobian(const int i, const UM::vec3& v, const int h, const int t, UM::mat3x3& J, UM::mat3x3& K);
	double evaluate_vert_energy(const int i, const UM::vec3& v, double& minDetJ);
	double compute_vert_grad_hess(const int i, UM::vec3& grad, UM::mat3x3& hess);
	double linesearch(const int i, const UM::vec3& dir, const double prev_E, double& new_E);
	UM::vec3 place_on_bnd(const int i, const UM::vec3& v);

	void update_connectivity();
	void update_eps(const int i, const double prev_E, const double new_E);
	void update_order();

	std::vector<int> order_;
	std::vector<bool> lock_;
	std::vector<double> vert_eps_;
	std::vector<double> vert_min_det_;
	std::vector<int> vert_type_;
	std::vector<std::vector<int>> vert_one_ring_;
	std::vector<std::vector<int>> vert_bnd_projection_triangles_;
	std::vector<std::vector<int>> vert_bnd_projection_segments_;
	std::vector<int> vert_bnd_projection_point_;
	std::vector<std::array<UM::vec3, 4>> hexrefcoef_;
	double scale_;


};




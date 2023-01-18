#include <ultimaille/all.h>
#include <array>


static constexpr std::array<std::array<int, 4>, 8> HEX_TO_SJ_TET_SPLIT = { {
	{0,1,2,4},
	{1,3,0,5},
	{2,0,3,6},
	{3,2,1,7},
	{4,6,5,0},
	{5,4,7,1},
	{6,7,4,2},
	{7,5,6,3},
} };

class vertex_smoother {	
public:
	vertex_smoother(UM::Hexahedra& m);

	void set_locked_vertices(const std::vector<bool>& locks);

	double _min_SJ_project_ = 0.1; // won't project if quality is worse than that
	double _min_SJ_update = 0.7; // will try to improve elements with quality worse than that
	double _theta = 1e-3; // elliptic smoothing parameter
	double _eps = 1e-5; // elliptic smoothing parameter

	/* 
	to set projective data:
	1. first set  UM::Triangles bnd and  UM::PolyLine ft which contains all
		the data needed for projection.
	2. for each vertex of your mesh you want projected, set-up the type of projection
		with the entities the vertex can be projected to. 
	*/
	void set_bnd_triangles(const UM::Triangles& bnd);
	void set_features_segment(const UM::PolyLine& ft);

	void set_vertex_triangles(const int i, const std::vector<int>& triangles);
	void set_vertex_segments(const int i, const std::vector<int>& segments);
	void set_vertex_point(const int i, const UM::vec3 v);


	void update_bad_elements();
	void update_order(); // prioritise vertex based on feature projection
	void run_iter();

	void execute(unsigned nb_iter) {
		update_bad_elements();
		update_order();
		for (unsigned i = 0; i < nb_iter; i++){
			std::cout << "Smoothing iter " << i + 1 << " / " << nb_iter << std::endl; 
			run_iter();
		}	
	}


	double compute_hex_minSJ(int h) const;

private:
	double compute_tet_SJ(int t) const;
	double compute_vert_minSJ(int v) const;
	void udpdate_tet_quality(int t);

	UM::Hexahedra& _m;

	std::vector<UM::vec3>& _pts;
	struct VertData {
		bool locked = false;
		bool is_bnd_compliant = false;
		int type = 0; // 0 is vol, 1 bnd, 2 curve, 3 point 
		int nb_bad_tets = 0;
		std::vector<int> verts;
		std::vector<int> tets;
	};
	std::vector<VertData> _vert_data;
	std::vector<int> _vert_order;

	struct TetData {
		bool is_bad = false;
		double SJ;
		std::array<unsigned, 4> verts;
		std::array<UM::vec3, 4> pre_computed;
	};
	
	
	std::vector<std::array<int, 8>> _hex2tets;
	std::vector<TetData> _tets;

	struct ProjectData {
		UM::Triangles bnd;
		UM::PolyLine ft;
		std::vector<UM::vec3> pts;
		std::vector<std::vector<int>> vert2triangles;
		std::vector<std::vector<int>> vert2segments;
		std::vector<int> vert2point;
	} _proj;

	UM::vec3 projected_position(int v) const;

	UM::vec3 compute_vert_elliptic_grad_hess(int v, UM::mat3x3& hess) const;
	UM::vec3 compute_vert_elliptic_grad_trunc_hess(int v, UM::mat3x3& hess) const;
	UM::vec3 compute_naive_laplacian_direction(int v) const;


};




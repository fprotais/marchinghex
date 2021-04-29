#include <iostream>
#include <cstdlib>
#include <ultimaille/all.h>

#include "marchingcube.h"

#include "locale_iterative_smoother.h"
#include "quality.h"
using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)



inline void clean_outside(Hexahedra& hex, const std::vector<bool>& is_in, std::vector<Marchinghex::VERT_TYPE>& vert_type, std::vector<vec3>& wish) {

    {
        std::vector<bool> tokill(hex.ncells());
        FOR(i, hex.ncells()) tokill[i] = !is_in[i];
        hex.delete_cells(tokill);
    }


    std::vector<bool> tokill(hex.nverts(), true);
    FOR(h, hex.ncells()) FOR(hv, 8)  tokill[hex.vert(h, hv)] = false;
    std::vector<int> old2new;
    hex.points.delete_points(tokill, old2new);
    FOR(c, hex.ncells()) FOR(cv, 8) hex.vert(c,cv) = old2new[hex.vert(c, cv)];

    FOR(v, vert_type.size()) if (old2new[v] != -1) vert_type[old2new[v]] = vert_type[v];
    vert_type.resize(hex.nverts());
    FOR(v, wish.size()) if (old2new[v] != -1) wish[old2new[v]] = wish[v];
    wish.resize(hex.nverts());
}


// return number of unlocked vert
int lock_far_from_bnd(const UM::Hexahedra& m, iterative_smoother& smoother, const std::vector<Marchinghex::VERT_TYPE>& vert_type, int dist = 5) {
    std::vector<int> cell_dist(m.ncells(), 20000);
    for (int c : range(m.ncells())) for (int cv : range(8)) {
        if ((int)vert_type[m.vert(c, cv)] > 0) cell_dist[c] = 0;
    }
    VolumeConnectivity vec(m);
    for (int i : range(dist)) {
        for (int c : range(m.ncells())) if (cell_dist[c] == i) {
            for (int cf : range(6)) if (vec.adjacent[6 * c + cf] != -1) {
                int c2 = vec.adjacent[6 * c + cf] / 6;
                if (cell_dist[c2] > i) cell_dist[c2] = i + 1;
            }
        }
    }
    std::vector<bool> locked(m.nverts(), false);
    for (int c : range(m.ncells())) if (cell_dist[c] > dist) {
        for (int cv : range(8)) locked[m.vert(c, cv)] = true;
    }
    smoother.set_locks(locked);
    int cnt = 0;
    for (int v : range(m.nverts())) if (!locked[v]) cnt++;
    return cnt;
}


int main(int argc, char** argv) {

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " grid.ext domain.ext hexmesh.ext " << std::endl;
        std::cerr << "Input:" << std::endl;
        std::cerr << "- grid.ext must contain an hex mesh" << std::endl;
        std::cerr << "- domain.ext must contain a tet mesh, and can contain triangles and edges." << std::endl;
        std::cerr << "  > triangle will be taken as domain boundary, and edges as its features." << std::endl;
        std::cerr << "Output: [OPTIONAL]" << std::endl;
        std::cerr << "- hexmesh.ext is the results hexahedral mesh (default: hexmesh.mesh)" << std::endl;

        std::cerr << "ext formats are ultimaille supported volumic formats (geogram, medit -.mesh- and vtk)." << std::endl;
        std::cerr << "For more details, see README." << std::endl;
        std::cerr << "contact: francois.protais@inria.fr" << std::endl;
        return 1;
    }

    std::string gridname = argv[1];
    std::string domainname = argv[2];
    std::string hexmeshname = "hexmesh.mesh";

    if (argc > 3) hexmeshname = argv[3];



    Hexahedra grid; read_by_extension(gridname, grid);

    Tetrahedra domain; read_by_extension(domainname, domain);
    Triangles boundary; read_by_extension(domainname, boundary);
    PolyLine features; read_by_extension(domainname, features);

    std::cerr << "Starting marchinghex pattern extraction." << std::endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    grid_pattern_extractor extractor(grid);
    extractor.set_domain(domain);
    extractor.set_boundary(boundary);
    extractor.set_hard_edges(features);

    Hexahedra hex;
    std::vector<bool> is_in;
    std::vector<Marchinghex::VERT_TYPE> vert_type;
    std::vector<vec3> wish;
    extractor.extract_patterns(hex, is_in, wish, vert_type);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cerr << "Marching Hex Run Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000. << "sec." << std::endl;


    clean_outside(hex, is_in, vert_type, wish);
    CellAttribute<double> hex_sj(hex);
    double min_sj = compute_scaled_jacobian(hex, hex_sj);
    write_by_extension(hexmeshname, hex, VolumeAttributes{ {},{{"sj", hex_sj.ptr}},{},{} });

    std::cerr << "Grid has: " << grid.nverts() << " verts and " << grid.ncells() << " hexahedra." << std::endl;
    std::cerr << "Hexmesh has: " << hex.nverts() << " verts and " << hex.ncells() << " hexahedra." << std::endl;
    std::cerr << "Min Scaled jacobian: " << min_sj << std::endl;
    std::cerr << "Saving marchinghex result..." << std::endl;

    std::cerr << "Beguinning smoothing of the result. It is an iterative process, you can stop it when it is not saving." << std::endl;
    std::cerr << "The smoothing time is linear to the side of the hexmesh (5000 nodes/sec for each iteration on my computer). Speeding it up is a work in progress." << std::endl;
    write_by_extension("mh_result.mesh", hex);

    begin = std::chrono::steady_clock::now();
    iterative_smoother smoother(hex);
    int nb_unlocked = lock_far_from_bnd(hex, smoother, vert_type);
    std::cerr << "Number of vertices optimized: " << nb_unlocked << std::endl;
    smoother.set_bnd_triangles(boundary);
    smoother.set_features_segment(features);
    boundary_matcher matcher(boundary, features);
    double grid_size_edge = 0;
    FOR(h, grid.ncells()) FOR(hf, 6) FOR(hfv, 4)
        grid_size_edge = std::max(grid_size_edge,(grid.points[grid.facet_vert(h, hf, hfv)] - grid.points[grid.facet_vert(h, hf, (hfv + 1) % 3)]).norm());

    FOR(i, hex.nverts()) {
        if (vert_type[i] == Marchinghex::VERT_IS_INSIDE) continue;
        if (vert_type[i] == Marchinghex::VERT_IS_ON_BOUNDARY) {
            std::vector<int> tri;
            matcher.get_vert_close_triangles(wish[i], tri, 8*grid_size_edge);
            smoother.set_vertex_triangles(i, tri);
        }
        if (vert_type[i] == Marchinghex::VERT_IS_ON_FEATURE) {
            std::vector<int> seg;
            matcher.get_vert_close_segments(wish[i], seg, 2*grid_size_edge);
            smoother.set_vertex_segments(i, seg);
        }
        if (vert_type[i] == Marchinghex::VERT_IS_FEATURE_POINT) {
            smoother.set_vertex_point(i, wish[i]);
        }
    }

    FOR(i, 20) {
        std::cerr << "Iter - " << i << std::endl;
        begin = std::chrono::steady_clock::now();
        smoother.run_iter();
        smoother.scale_back();
        min_sj = compute_scaled_jacobian(hex, hex_sj);
        std::cerr << "Min Scaled jacobian: " << min_sj << std::endl;
        std::cerr << "Saving, do not quit...";
        write_by_extension(hexmeshname, hex, VolumeAttributes{ {},{{"sj", hex_sj.ptr}},{},{} });
        write_by_extension("iter_"+std::to_string(i)+".mesh", hex);
        std::cerr << "Done." << std::endl;
        smoother.scale_up();
        end = std::chrono::steady_clock::now();
        std::cerr << "Iter time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000. << "sec." << std::endl;

    }
    smoother.scale_back();

    std::cerr << "FINISHED." << std::endl;
    double dist = hausdorff_dist(hex, domain);
    std::cerr << "Hausdorff dist: " << dist << std::endl;
    if (domain.nverts() == 0) return 0;
    vec3 min = domain.points[0], max = domain.points[0];
    FOR(v, domain.nverts()) FOR(d, 3) {
        min[d] = std::min(min[d], domain.points[v][d]);
        max[d] = std::max(max[d], domain.points[v][d]);
    }
    double bboxdiagsize = (max - min).norm();
    std::cerr << "HR = Hausdorff dist / BBOX diagonal =  " << 100 * dist / bboxdiagsize << " %" << std::endl;



}

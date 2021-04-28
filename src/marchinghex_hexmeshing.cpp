#include <iostream>
#include <cstdlib>
#include <ultimaille/all.h>

#include "marchingcube.h"

#include "locale_iterative_smoother.h"

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
    write_by_extension(hexmeshname, hex);
    std::cerr << "Grid has: " << grid.nverts() << " verts and " << grid.ncells() << " hexahedra." << std::endl;
    std::cerr << "Hexmesh has: " << hex.nverts() << " verts and " << hex.ncells() << " hexahedra." << std::endl;
    std::cerr << "Saving marchinghex result..." << std::endl;

    std::cerr << "Beguinning smoothing of the result. It is an iterative process, you can stop it when it is not saving." << std::endl;
    std::cerr << "The smoothing time is linear to the side of the hexmesh (5000 nodes/sec for each iteration on my computer). Speeding it up is a work in progress." << std::endl;

    begin = std::chrono::steady_clock::now();
    iterative_smoother smoother(hex);
    smoother.set_bnd_triangles(boundary);
    smoother.set_features_segment(features);
    boundary_matcher matcher(boundary, features);
    double grid_size_edge = 0;
    FOR(h, domain.ncells()) FOR(hf, 6) FOR(hfv, 4)
        grid_size_edge += (domain.points[domain.facet_vert(h, hf, hfv)] - domain.points[domain.facet_vert(h, hf, (hfv + 1) % 3)]).norm();
    grid_size_edge /= domain.ncells() * 24;

    FOR(i, hex.nverts()) {
        if (vert_type[i] == Marchinghex::VERT_IS_INSIDE) continue;
        if (vert_type[i] == Marchinghex::VERT_IS_ON_BOUNDARY) {
            std::vector<int> tri;
            matcher.get_vert_close_triangles(wish[i], tri, grid_size_edge);
            smoother.set_vertex_triangles(i, tri);
        }
        if (vert_type[i] == Marchinghex::VERT_IS_ON_FEATURE) {
            std::vector<int> seg;
            matcher.get_vert_close_segments(wish[i], seg, grid_size_edge);
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
        std::cerr << "Saving, do not quit...";
        write_by_extension(hexmeshname, hex);
        std::cerr << "Done." << std::endl;
        smoother.scale_up();
        end = std::chrono::steady_clock::now();
        std::cerr << "Iter time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000. << "sec." << std::endl;

    }

    std::cerr << "FINISHED." << std::endl;

}

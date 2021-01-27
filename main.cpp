#include <iostream>
#include <cstdlib>
#include <ultimaille/io/fromsuffix.h>

#include <ultimaille/geometry.h>

#include <chrono>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "marchingcube.h"
#include "smoother.h"
using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)

void get_ultimaille_mesh(const std::vector<vec3>&verts_, const std::vector<int>&hexes_, Hexahedra & m) {
    m.points.create_points(verts_.size());
    FOR(v, verts_.size()) m.points[v] = verts_[v];
    m.create_cells(hexes_.size() / 8);
    FOR(h, m.ncells()) FOR(hv, 8) m.vert(h, hv) = hexes_[8 * h + hv];
}

void reading_with_ultimaille(const std::string& filename, std::vector<vec3>& verts, std::vector<int>& tets, std::vector<int>& tris, std::vector<int>& edges) {
    Tetrahedra m;
    Triangles bnd;
    PolyLine hd;
    VolumeAttributes tmp = read_fromsuffix(filename, m);
    SurfaceAttributes tmp1 = read_fromsuffix(filename, bnd);
    PolyLineAttributes tmp2 = read_fromsuffix(filename, hd);

    write_fromsuffix("input.mesh", m);
    write_fromsuffix("boundary.mesh", bnd);
    write_fromsuffix("hardedges.mesh", hd);
    write_fromsuffix("input.vtk", m);
    write_fromsuffix("boundary.vtk", bnd);
    write_fromsuffix("hardedges.vtk", hd);
    verts.resize(m.nverts());
    FOR(v, m.nverts()) verts[v] = m.points[v];
    tets.resize(4 * m.ncells());
    FOR(t, m.ncells()) FOR(i, 4) tets[4 * t + i] = m.vert(t, i);
    tris.resize(3 * bnd.nfacets());
    FOR(t, bnd.nfacets()) FOR(i, 3) tris[3 * t + i] = bnd.vert(t, i);
    edges.resize(2 * hd.nsegments());
    FOR(e, hd.nsegments()) FOR(i, 2) edges[2 * e + i] = hd.vert(e, i);
}

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.mesh" << std::endl;
        std::cerr << "Alternatively, model.geogram" << std::endl;
        std::cerr << "Model must be a tet mesh." << std::endl;
        std::cerr << "Opt : Arg2 is scale (double), default: 1." << std::endl;
        std::cerr << "Opt : Arg3 is output name, default: res.mesh" << std::endl;
        std::cerr << "Full use: " << argv[0] << " model.mesh 1 res.mesh" << std::endl;
        return 1;
    }
    double scale = 1;
    std::string outfile = "res.vtk";
    if (argc > 2) scale = std::stod(argv[2]);
    if (argc > 3) outfile = argv[3];
    std::cerr << "scale: " << scale << std::endl;

    std::cerr << "Loading File:" << argv[1] << std::endl;

    std::vector<vec3> verts;
    std::vector<int> tets;
    std::vector<int> tris;
    std::vector<int> edges;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    reading_with_ultimaille(argv[1], verts, tets, tris, edges);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cerr << "File loaded: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()/1000. << "sec." << std::endl;
    std::cerr << tets.size() << std::endl;
    std::cerr << "Nb vertices: " << verts.size() << std::endl;
    std::cerr << "Nb tets: " << tets.size() / 4 << std::endl;

    std::cerr << "Rescaling" << std::endl;

    FOR(v, verts.size()) verts[v] += vec3(0.27, 0.27, 0.27);
    double average_edge_size = 0;
    FOR(t, tets.size() / 4) average_edge_size += (verts[tets[4 * t + 0]] - verts[tets[4 * t + 1]]).norm();
    average_edge_size /= tets.size() / 4;
    FOR(v, verts.size()) verts[v] = scale * verts[v] / average_edge_size;


    std::vector<int> mc_hexes;
    std::vector<vec3> mc_verts;
    marchingcube::hexify_with_bnd(verts, tets, tris, edges, mc_verts, mc_hexes);

    FOR(v, mc_verts.size()) mc_verts[v] = mc_verts[v] * average_edge_size / scale;
    FOR(v, mc_verts.size()) mc_verts[v] += -vec3(0.27, 0.27, 0.27);
    std::cerr << "Nb vertices: " << mc_verts.size() << std::endl;
    std::cerr << "Nb Hexes: " << mc_hexes.size() / 8 << std::endl;
    std::cerr << "Writting " << outfile << std::endl;
    begin = std::chrono::steady_clock::now();
    Hexahedra m;
    get_ultimaille_mesh(mc_verts, mc_hexes, m);
    write_fromsuffix(outfile, m);

    elliptic_smoother(m);
    write_fromsuffix("smoothed.mesh", m);

    end = std::chrono::steady_clock::now();
    std::cerr << "File writen: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000. << "sec." << std::endl;

    return 0;
}

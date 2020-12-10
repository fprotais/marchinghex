#include <iostream>
#include <cstdlib>
#include <ultimaille/io/medit.h>
#include <ultimaille/io/geogram.h>
#include <ultimaille/geometry.h>

#include <chrono>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "marchingcube.h"

using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)


void write_medit_format(const std::string& filename, std::vector<vec3> verts_, std::vector<int> hexes_) {
    Hexahedra m;
    m.points.create_points(verts_.size());
    FOR(v, verts_.size()) m.points[v] = verts_[v];
    m.create_cells(hexes_.size() / 8);
    FOR(h, m.ncells()) FOR(hv, 8) m.vert(h, hv) = hexes_[8 * h + hv];
    write_medit(filename, m);
}

void mcreadfile(const std::string& filename, std::vector<vec3>& verts, std::vector<int>& tets) {
    Tetrahedra m;
    if (std::string(filename.end() - 5, filename.end()) == std::string(".mesh")) {
        read_medit(filename, m);
    }
    else if (std::string(filename.end() - 8, filename.end()) == std::string(".geogram")) {
        VolumeAttributes tmp = read_geogram(filename, m);
    }
    else {
        std::cerr << "File format suffix not supported, use .geogram or .mesh" << std::endl;
        exit(1);
    }
    verts.resize(m.nverts());
    FOR(v, m.nverts()) verts[v] = m.points[v];
    tets.resize(4 * m.ncells());
    FOR(t, m.ncells()) FOR(i, 4) tets[4 * t + i] = m.vert(t, i);
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
    double scale = 0.75;
    std::string outfile = "res.mesh";
    if (argc > 2) scale = std::stod(argv[2]);
    if (argc > 3) outfile = argv[3];
    std::cerr << "scale: " << scale << std::endl;

    std::cerr << "Loading File:" << argv[1] << std::endl;

    std::vector<vec3> verts;
    std::vector<int> tets;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    mcreadfile(argv[1], verts, tets);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cerr << "File loaded: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()/1000. << "sec." << std::endl;
    std::cerr << tets.size() << std::endl;
    std::cerr << "Nb vertices: " << verts.size() << std::endl;
    std::cerr << "Nb tets: " << tets.size() / 4 << std::endl;

    std::cerr << "Rescaling" << std::endl;
    double average_edge_size = 0;
    FOR(t, tets.size() / 4) average_edge_size += (verts[tets[4 * t + 0]] - verts[tets[4 * t + 1]]).norm();
    average_edge_size /= tets.size() / 4;
    FOR(v, verts.size()) verts[v] = scale * verts[v] / average_edge_size;


    std::vector<int> mc_hexes;
    std::vector<vec3> mc_verts;
    marchingcube::hexify(verts, tets, mc_verts, mc_hexes);


    std::cerr << "Nb vertices: " << mc_verts.size() << std::endl;
    std::cerr << "Nb Hexes: " << mc_hexes.size() / 8 << std::endl;
    std::cerr << "Writting " << outfile << std::endl;
    begin = std::chrono::steady_clock::now();
    write_medit_format(outfile, mc_verts, mc_hexes);
    end = std::chrono::steady_clock::now();
    std::cerr << "File writen: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000. << "sec." << std::endl;

    return 0;
}

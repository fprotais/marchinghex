#include <iostream>
#include <cstdlib>
#include <ultimaille/mesh_io.h>
#include <ultimaille/surface.h>
#include <ultimaille/polyline.h>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "marchingcube.h"

using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)

inline void file_must_no_be_at_end(std::ifstream & f, const std::string& reason = " should'nt") {
    if (f.eof()) {
        f.close();
        std::cout << "File ended to soon while : " << reason << std::endl;
        exit(1);
    }
}
inline static bool string_start(const std::string& string, const std::string& start_of_string) {
    size_t start = 0;
    FOR(i, string.size()) if (string[i] != ' ' && string[i] != '\t') {
        start = (size_t)i;
        break;
    }
    std::string copy_without_space(string.begin() + start, string.end());
    if (copy_without_space.size() < start_of_string.size()) return false;
    return (std::string(copy_without_space.begin(), copy_without_space.begin() + (long int)start_of_string.size()) == start_of_string);
}
void read_medit_format(const std::string& filename, std::vector<vec3>& verts_,  std::vector<int>& tets_) {
    //std::vector<vec3> verts_;
    std::vector<int> edges_;
    std::vector<int> tris_;
    std::vector<int> quads_;
    //std::vector<int> tets_;
    std::vector<int> hexes_;

    std::ifstream in;
    in.open(filename, std::ifstream::in);
    if (in.fail()) {
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }

    std::string firstline;

    while (!in.eof()) {
        std::getline(in, firstline);
        if (string_start(firstline, "Vertices")) {
            std::string line;
            int nb_of_vertices = 0;
            {
                file_must_no_be_at_end(in, "parsing vertices");
                std::getline(in, line);
                std::istringstream iss(line.c_str());
                iss >> nb_of_vertices;
            }
            verts_.resize(nb_of_vertices);
            FOR(v, nb_of_vertices) {
                file_must_no_be_at_end(in, "parsing vertices");
                std::getline(in, line);
                std::istringstream iss(line.c_str());
                FOR(i, 3 )  iss >> verts_[v][i];
            }
        }
        if (string_start(firstline, "Edges")) {
            std::string line;
            int nb_of_edges = 0;
            {
                file_must_no_be_at_end(in, "parsing Edges");
                std::getline(in, line);
                std::istringstream iss(line.c_str());
                iss >> nb_of_edges;
            }
            edges_.resize(2 * nb_of_edges);
            FOR(e, nb_of_edges) {
                file_must_no_be_at_end(in, "parsing Edges");
                std::getline(in, line);
                std::istringstream iss(line.c_str());
                FOR(i, 2) {
                    int a = 0;
                    iss >> a; 
                    edges_[2 * e + i] = a - 1;
                }
            }
        }
        if (string_start(firstline, "Triangles")) {
            std::string line;
            int nb_of_tri = 0;
            {
                file_must_no_be_at_end(in, "parsing Triangles");
                std::getline(in, line);
                std::istringstream iss(line.c_str());
                iss >> nb_of_tri;
            }
            tris_.resize(3 * nb_of_tri);
            FOR(t, nb_of_tri) {
                file_must_no_be_at_end(in, "parsing Triangles");
                std::getline(in, line);
                std::istringstream iss(line.c_str());
                FOR(i, 3) {
                    int a = 0;
                    iss >> a;
                    tris_[3 * t + i] = a - 1;
                }
            }
        }
        if (string_start(firstline, "Quadrilaterals")) {
            std::string line;
            int nb_of_quads = 0;
            {
                file_must_no_be_at_end(in, "parsing Quadrilaterals");
                std::getline(in, line);
                std::istringstream iss(line.c_str());
                iss >> nb_of_quads;
            }
            quads_.resize(4 * nb_of_quads);
            FOR(q, nb_of_quads) {
                file_must_no_be_at_end(in, "parsing Quadrilaterals");
                std::getline(in, line);
                std::istringstream iss(line.c_str());
                FOR(i, 4) {
                    int a = 0;
                    iss >> a;
                    quads_[4 * q + i] = a - 1;
                }
            }
        }
        if (string_start(firstline, "Tetrahedra")) {
            std::string line;
            int nb_of_tets = 0;
            {
                file_must_no_be_at_end(in, "parsing Tetrahedra");
                std::getline(in, line);
                std::istringstream iss(line.c_str());
                iss >> nb_of_tets;
            }
            tets_.resize(4 * nb_of_tets);
            FOR(t, nb_of_tets) {
                file_must_no_be_at_end(in, "parsing Tetrahedra");
                std::getline(in, line);
                std::istringstream iss(line.c_str());
                FOR(i, 4) {
                    int a = 0;
                    iss >> a;
                    tets_[4 * t + i] = a - 1;
                }
            }
        }
        if (string_start(firstline, "Hexahedra")) {
            std::string line;
            int nb_of_hexs = 0;
            {
                file_must_no_be_at_end(in, "parsing Hexahedra");
                std::getline(in, line);
                std::istringstream iss(line.c_str());
                iss >> nb_of_hexs;
            }
            hexes_.resize(8 * nb_of_hexs);
            FOR(h, nb_of_hexs) {
                file_must_no_be_at_end(in, "parsing Hexahedra");
                std::getline(in, line);
                std::istringstream iss(line.c_str());
                FOR(i, 8) {
                    int a = 0;
                    iss >> a;
                    hexes_[8 * h + i] = a - 1;
                }
            }
        }

    }

}

void write_medit_format(const std::string& filename, std::vector<vec3> verts_, std::vector<int> hexes_) {
    std::ofstream out_f;
    out_f.open(filename, std::ifstream::out);
    if (out_f.fail()) {
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }
    std::stringstream out;
    out << std::fixed << std::setprecision(4);
    out << "MeshVersionFormatted 2" << std::endl << std::endl;
    out << "Dimension" << std::endl << "3" <<  std::endl << std::endl;

    out << "Vertices" << std::endl;
    out << verts_.size() << std::endl;
    FOR(v, verts_.size()) {
        FOR(d, 3) out << verts_[v][d] << " ";
        out << "1" << std::endl;
    }
    out << std::endl;
    out << "Hexahedra" << std::endl;
    out << hexes_.size() / 8 << std::endl;
    FOR(h, hexes_.size() / 8) {
        out << hexes_[8 * h + 1] + 1 << " ";
        out << hexes_[8 * h + 0] + 1 << " ";
        out << hexes_[8 * h + 2] + 1 << " ";
        out << hexes_[8 * h + 3] + 1 << " ";
        out << hexes_[8 * h + 5] + 1 << " ";
        out << hexes_[8 * h + 4] + 1 << " ";
        out << hexes_[8 * h + 6] + 1 << " ";
        out << hexes_[8 * h + 7] + 1 << " ";
        //FOR(i, 8) out << hexes_[8*h + i]+1 << " ";
        out << "1" << std::endl;
    }
    out << "end" << std::endl;

    out_f << out.rdbuf();
    out_f.close();
}

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.mesh" << std::endl;
        return 1;
    }
    std::vector<vec3> verts;
    std::vector<int> tets;
    read_medit_format(argv[1], verts, tets);
    double average_edge_size = 0;
    FOR(t, tets.size() / 4) average_edge_size += (verts[tets[4 * t + 0]] - verts[tets[4 * t + 1]]).norm();
    average_edge_size /= tets.size() / 4;
    FOR(v, verts.size()) verts[v] = verts[v] / average_edge_size;

    std::cerr << "File loaded" << std::endl;
    std::cerr << tets.size() << std::endl;
    std::vector<int> mc_hexes;
    std::vector<vec3> mc_verts;
    marchingcube::hexify(verts, tets, mc_verts, mc_hexes);
    std::cerr << "Writting res.mesh" << std::endl;

    write_medit_format("res.mesh", mc_verts, mc_hexes);
    return 0;
}

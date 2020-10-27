#include <iostream>
#include <cstdlib>
#include <ultimaille/mesh_io.h>
#include <ultimaille/surface.h>
#include <ultimaille/polyline.h>

using namespace UM;

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }
    {
        PolyLine pl;
        PolyLineAttributes attributes = read_geogram(argv[1], pl);
        write_geogram("pl_write_test.geogram", pl, attributes);
    }
    return 0;
}

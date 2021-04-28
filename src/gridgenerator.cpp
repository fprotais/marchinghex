#include <ultimaille/all.h>

using namespace UM;
#define FOR(i, n) for(int i = 0; i < n; i++)





int main(int argc, char** argv) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " domain.ext grid.ext scale" << std::endl;
        std::cerr << "Input:" << std::endl;
        std::cerr << "- domain.ext must contain a tet mesh" << std::endl;
        std::cerr << "Output: [OPTIONAL]" << std::endl;
        std::cerr << "- grid.ext is an hex mesh containing domain, with cube size similar to tet size. (default: grid.mesh)" << std::endl;
        std::cerr << "- scale is a double. 5 is very fine grid, 0.1 is very coarse. (default: 1.)" << std::endl;

        std::cerr << "ext formats are ultimaille supported volumic formats (geogram, medit -.mesh- and vtk)." << std::endl;
        std::cerr << "For more details, see README." << std::endl;
        std::cerr << "contact: francois.protais@inria.fr" << std::endl;
        return 1;
    }

    std::string domainname = argv[1];
    std::string gridname = "grid.mesh";
    double scale = 1;

    if (argc > 2) gridname = argv[2];
    if (argc > 3) scale = std::stod(argv[3]);
    
    Tetrahedra domain; read_by_extension(domainname, domain);
    UM::Hexahedra grid;

    if (scale <= 0) {
        write_by_extension(gridname, grid);
        return 0;
    }


    double avg_edge = 0;
    FOR(tet, domain.ncells()) FOR(tet_f, 4) FOR(tfv, 3) 
        avg_edge += (domain.points[domain.facet_vert(tet, tet_f, tfv)] - domain.points[domain.facet_vert(tet, tet_f, (tfv+1)%3)]).norm();
    avg_edge /= domain.ncells() * 12;
    scale /= avg_edge;

    vec3 max = { -1E10,-1E10,-1E10 };
    vec3 min = { 1E10,1E10,1E10 };
    FOR(v, domain.nverts()) FOR(d, 3) {
        if (domain.points[v][d] > max[d]) max[d] = domain.points[v][d];
        if (domain.points[v][d] < min[d]) min[d] = domain.points[v][d];
    }
    max = max * scale;
    min = min * scale;

    int nb_in_X = (int)std::ceil(max[0] - min[0]) + 3;
    int nb_in_Y = (int)std::ceil(max[1] - min[1]) + 3;
    int nb_in_Z = (int)std::ceil(max[2] - min[2]) + 3;
    int start_in_X = (int)std::floor(min[0]) - 1;
    int start_in_Y = (int)std::floor(min[1]) - 1;
    int start_in_Z = (int)std::floor(min[2]) - 1;
    grid.points.create_points(nb_in_X * nb_in_Y * nb_in_Z);
    for (int Z = 0; Z < nb_in_Z; Z++)
        for (int Y = 0; Y < nb_in_Y; Y++)
            for (int X = 0; X < nb_in_X; X++)
                grid.points[Z * (nb_in_X * nb_in_Y) + Y * nb_in_X + X] = vec3(X + start_in_X, Y + start_in_Y, Z + start_in_Z);


    grid.create_cells((nb_in_X - 1) * (nb_in_Y - 1) * (nb_in_Z - 1));

    for (int Z = 0; Z < nb_in_Z - 1; Z++)
        for (int Y = 0; Y < nb_in_Y - 1; Y++)
            for (int X = 0; X < nb_in_X - 1; X++)
                FOR(Xp, 2) FOR(Yp, 2) FOR(Zp, 2)
                grid.vert(Z * ((nb_in_X - 1) * (nb_in_Y - 1)) + Y * (nb_in_X - 1) + X, 4 * Zp + 2 * Yp + Xp) = (Z + Zp) * (nb_in_X * nb_in_Y) + (Y + Yp) * nb_in_X + X + Xp;


    FOR(v, grid.nverts()) FOR(d, 3) grid.points[v][d] = grid.points[v][d] + 0.271;
    FOR(v, grid.nverts()) grid.points[v] = grid.points[v] / scale;

    write_by_extension(gridname, grid);
    return 0;
}
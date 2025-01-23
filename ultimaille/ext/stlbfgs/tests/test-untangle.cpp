#include <catch2/catch.hpp>

#include <iostream>
#include <fstream>
#include <cstring>
#include <limits>
#define _USE_MATH_DEFINES
#include <cmath>
#undef NDEBUG
#include <cassert>

#include "stlbfgs.h"

struct QuadMesh {
    QuadMesh() { // build a Belinsky's Z quad mesh
        int size = 8;
        boundary = std::vector<bool>(size*size);
        for (int j=0; j<size; j++) {
            for (int i=0; i<size; i++) {
                double x =   i/(double)size;
                double y = 2*j/(double)size;
                if (j>=size/2) {
                    x += 3./5;
                    y -= 3./5;
                }
                points.push_back(x);
                points.push_back(y);
                boundary[i+j*size] = (i==0 || i==size-1 || j==0 || j==size-1);
            }
        }
        for (int j=0; j<size-1; j++) {
            for (int i=0; i<size-1; i++) {
                quads.push_back(i  +(j  )*size);
                quads.push_back(i+1+(j  )*size);
                quads.push_back(i+1+(j+1)*size);
                quads.push_back(i  +(j+1)*size);
            }
        }
    }

    int nverts() const { return static_cast<int>(points.size()/2); }
    int nquads() const { return static_cast<int>(quads.size()/4);  }

    std::vector<int> quads;
    std::vector<double> points;
    std::vector<bool> boundary;
};

std::ostream& operator<<(std::ostream& out, const QuadMesh &m) {
    for (int v=0; v<m.nverts(); v++)
        out << "v " << m.points[v*2] << " " << m.points[v*2+1] << std::endl;
    for (int q=0; q<m.nquads(); q++) {
        out << "f ";
        for (int lv=0; lv<4; lv++) {
            out << (m.quads[q*4+lv]+1) << " ";
        }
        out << std::endl;
    }
    return out;
}

constexpr double Q[4][4][2] = { // quadratures for every quad corner
    {{-1,-1},{1,0},{0,0},{0,1}},
    {{-1,0},{1,-1},{0,1},{0,0}},
    {{0,0},{0,-1},{1,1},{-1,0}},
    {{0,-1},{0,0},{1,0},{-1,1}}
};

// evaluate the Jacobian matrix at the given quadrature point, return the Jacobian determinant
double jacobian(const std::vector<int> &quads, const std::vector<double>& x, const int q, const int qc, double J[2][2]) {
     double A[2][4] = {
         {x[quads[q*4+0]*2+0], x[quads[q*4+1]*2+0], x[quads[q*4+2]*2+0], x[quads[q*4+3]*2+0]},
         {x[quads[q*4+0]*2+1], x[quads[q*4+1]*2+1], x[quads[q*4+2]*2+1], x[quads[q*4+3]*2+1]}
     };
     const double (&B)[4][2] = Q[qc];
     for (int j=0; j<2; j++) // J = A*B
         for (int i=0; i<2; i++) {
            J[j][i] = 0;
            for (int k=0; k<4; k++)
                J[j][i] += A[j][k]*B[k][i];
         }
     return J[0][0]*J[1][1] - J[0][1]*J[1][0];
}

TEST_CASE("Quad mesh untangling", "[L-BFGS]") {
    QuadMesh m;
    int nvars = 2*m.nverts();

    double mindet = 0.;
    for (int iter=0; iter<10; iter++) {
        mindet = std::numeric_limits<double>::max();
        for (int q=0; q<m.nquads(); q++) // compute min det over every quad corner
            for (int qc=0; qc<4; qc++) {
                double J[2][2];
                double det = jacobian(m.quads, m.points, q, qc, J);
                mindet = std::min(mindet, det);
            }

        double eps = std::sqrt(pow(1e-6, 2) + .04*pow(std::min(mindet, 0.), 2)); // the regularization parameter e

        // compute the energy and its gradient
        const STLBFGS::Optimizer::func_grad_eval energy = [&](const std::vector<double>& x, double& f, std::vector<double>& g) {
            f = 0;
            g = std::vector<double>(nvars, 0);

            for (int q=0; q<m.nquads(); q++) { // sum over all quads
                for (int qc=0; qc<4; qc++) {   // evaluate the Jacobian matrix for every quad corner
                    double J[2][2];
                    double det = jacobian(m.quads, x, q, qc, J);

                    double chi = det/2. + std::sqrt(eps*eps + det*det)/2.;    // the penalty function
                    double chip = .5 + det/(2.*std::sqrt(eps*eps + det*det)); // its derivative
                    double F = (J[0][0]*J[0][0] + J[1][1]*J[1][1] + J[1][0]*J[1][0] + J[0][1]*J[0][1])/chi; // quad corner shape quality
                    f += F;

                    double dfdjT[2][2] = { // gradient of the energy w.r.t. the Jacobian matrix entries
                        {
                            (2*J[0][0] - J[1][1]*F*chip)/chi,
                            (2*J[1][0] + J[0][1]*F*chip)/chi
                        },
                        {
                            (2*J[0][1] + J[1][0]*F*chip)/chi,
                            (2*J[1][1] - J[0][0]*F*chip)/chi
                        }
                    };
                    double dfdu[4][2] = {{0,0},{0,0},{0,0},{0,0}}; // chain rule for the actual variables
                    for (int j=0; j<4; j++) // dfdu = Q[qc]*dfdjT
                        for (int i=0; i<2; i++)
                            for (int k=0; k<2; k++)
                                dfdu[j][i] += Q[qc][j][k]*dfdjT[k][i];

                    for (int i=0; i<4; i++) {
                        if (m.boundary[m.quads[q*4+i]]) continue; // the boundary verts are locked
                        g[(m.quads[q*4+i])*2+0] += dfdu[i][0];
                        g[(m.quads[q*4+i])*2+1] += dfdu[i][1];
                    }
                }
            }
        };
        STLBFGS::Optimizer opt{energy};
        opt.run(m.points); // inner L-BFGS loop
    }

    std::ofstream ofs("Z.obj", std::ios::binary);
    ofs << m;
    ofs.close();

    CHECK( std::abs(mindet-0.0010) < 1e-4 );
}


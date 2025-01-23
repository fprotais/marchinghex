#include <catch2/catch.hpp>

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include "stlbfgs.h"

using namespace STLBFGS;
static const double xtol = 1e-3;

TEST_CASE("Powell badly scaled function", "[L-BFGS]") {
    const Optimizer::func_grad_eval fcn = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        double f0 = pow(10, 4)*x[0]*x[1]-1;
        double f1 = exp(-x[0])+exp(-x[1])-1.0001;
        f = f0*f0+f1*f1;

        double dx0f0 = 2*pow(10, 4)*f0*x[1];
        double dx0f1 = -2*f1*exp(-x[0]);
        double dx0f = dx0f0+dx0f1;
        double dx1f0 = 2*pow(10, 4)*f0*x[0];
        double dx1f1 = -2*f1*exp(-x[1]);
        double dx1f = dx1f0+dx1f1;
        g =  {dx0f, dx1f};
    };

    std::vector<double> x = {0, 1};
    Optimizer opt{fcn};
    opt.ftol = 1e-15;
    opt.run(x);

    CHECK( std::abs(x[0]-1.098e-5) < xtol );
    CHECK( std::abs(x[1]-9.106) < xtol );
}


TEST_CASE("Powell singular function", "[L-BFGS]") {
    const Optimizer::func_grad_eval fcn = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        double f0 = x[0]+10*x[1];
        double f1 = pow(5, 0.5)*(x[2]-x[3]);
        double f2 = pow(x[1]-2*x[2], 2);
        double f3 = pow(10, 0.5)*pow(x[0]-x[3], 2);
        f = f0*f0+f1*f1+f2*f2+f3*f3;

        double dx0f = 2*(x[0]+10*x[1])+ 4*10*pow(x[0]-x[3], 3);
        double dx1f = 20*(x[0]+10*x[1])+ 4*pow(x[1]-2*x[2], 3);
        double dx2f = 2*5*(x[2]-x[3])+ -8*pow(x[1]-2*x[2], 3);
        double dx3f = -2*5*(x[2]-x[3])+ -4*10*pow(x[0]-x[3], 3);
        g = {dx0f, dx1f, dx2f, dx3f};
    };

    std::vector<double> x = {3, -1, 0, 1};
    Optimizer opt{fcn};
    opt.ftol = 1e-12;

    opt.run(x);

    for (int i=0; i<4; i++) {
        CHECK( std::abs(x[i]) < xtol );
    }
}


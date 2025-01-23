#define _USE_MATH_DEFINES
#include <cmath>
#include <catch2/catch.hpp>

#include <iostream>
#include "stlbfgs.h"

using namespace STLBFGS;
static const double xtol = 1e-3;

TEST_CASE("Helical valley function", "[L-BFGS]") {
    const Optimizer::func_grad_eval fcn = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        const auto phi = [](double x0, double x1) {
            double l = 1/(2*M_PI)*atan(x1/x0);
            if (x0>0) return l;
            else return l+.5;
        };

        double f0 = 10*(x[2]-10*phi(x[0], x[1]));
        double f1 = 10*(pow(x[0]*x[0]+x[1]*x[1], .5) - 1);
        double f2 = x[2];
        f = f0*f0+f1*f1+f2*f2;


        double dx0phi = 1/(2*M_PI)*1/(1+pow(x[1]/x[0] ,2))* x[1]/pow(x[0], 2)*(-1);
        double dx0f = 2*10*(x[2]-10*phi(x[0], x[1]))* (-100)*dx0phi + 2*10*(pow(x[0]*x[0]+x[1]*x[1], .5) - 1)* 5*pow(pow(x[0], 2)+pow(x[1], 2), .5)* 2*x[0];

        double dx1phi = 1/(2*M_PI)*1/(1+pow(x[1]/x[0] ,2))*1/x[0];
        double dx1f = 2*10*(x[2]-10*phi(x[0], x[1]))* (-100)*dx1phi + 2*10*(pow(x[0]*x[0]+x[1]*x[1], .5) - 1)* 5*pow(pow(x[0], 2)+pow(x[1], 2), .5)*2*x[1];

        double dx2f = 2*10*(x[2]-10*phi(x[0], x[1]))*10 + 2*x[2];

        g = {dx0f, dx1f, dx2f};
    };

    std::vector<double> x = {-1, 0, 0};
    Optimizer opt{fcn};
    opt.run(x);

    CHECK( std::abs(x[0]-1) < xtol );
    CHECK( std::abs(x[1]-0) < xtol );
    CHECK( std::abs(x[2]-0) < xtol );
}




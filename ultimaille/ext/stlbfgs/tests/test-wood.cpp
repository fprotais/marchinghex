#include <catch2/catch.hpp>

#define _USE_MATH_DEFINES
#include <cmath>
#include "stlbfgs.h"

using namespace STLBFGS;
static const double xtol = 1e-3;

TEST_CASE("Wood function", "[L-BFGS]") {
    const Optimizer::func_grad_eval fcn = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        double f0 = 10*(x[1]-pow(x[0], 2));
        double f1 = 1-x[0];
        double f2 = std::sqrt(90)*(x[3]-pow(x[2], 2));
        double f3 = 1-x[2];
        double f4 = std::sqrt(10)*(x[1]+x[3]-2);
        double f5 = (x[1]-x[3])/std::sqrt(10.);
        f = f0*f0+f1*f1+f2*f2+f3*f3+f4*f4+f5*f5;

        double dx0f = -400*(x[1]-pow(x[0], 2))*x[0]-2*(1-x[0]);
        double dx1f = 200*(x[1]-pow(x[0], 2))+20*(x[1]+x[3]-2)+2/10.0*(x[1]-x[3]);
        double dx2f = -4*90*(x[3]-pow(x[2], 2))*x[2]-2*(1-x[2]);
        double dx3f = 2*90*(x[3]-pow(x[2], 2))+ 20*(x[1]+x[3]-2)-2/10.0*(x[1]-x[3]);

        g = {dx0f, dx1f, dx2f, dx3f};
    };

    std::vector<double> x = {-3, -1, -3, -1};
    Optimizer opt{fcn};
    opt.ftol = 1e-10;
    opt.run(x);

    for (int i=0; i<4; i++)
        CHECK(std::abs(x[i]-1.) < xtol);
}


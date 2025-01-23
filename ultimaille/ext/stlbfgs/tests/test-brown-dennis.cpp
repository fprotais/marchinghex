#include <catch2/catch.hpp>
#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>
#include "stlbfgs.h"

using namespace STLBFGS;
static const double xtol = 1e-5;

TEST_CASE("Brown and Dennis function", "[L-BFGS]") {
    const Optimizer::func_grad_eval fcn = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        f = 0;
        g = {0, 0, 0, 0};

        for (int i = 1; i <= 20; i++) {
            double c = i/5.;
            double f1 = x[0] + c*x[1] - std::exp(c);
            double f2 = x[2] + x[3]*std::sin(c) - std::cos(c);
            double t = f1*f1 + f2*f2;

            f += t*t;
            g[0] += 4 * t * f1;
            g[1] += 4 * t * f1 * c;
            g[2] += 4 * t * f2;
            g[3] += 4 * t * f2 * std::sin(c);
        }
    };

    std::vector<double> x = {25, 5, -5, -1};
    Optimizer opt{fcn};
    opt.ftol = 1e-8;
    opt.run(x);

    double f;
    std::vector<double> g;
    fcn(x, f, g);

    CHECK( std::abs(f-85822.2)<0.1 );

    std::vector<double> globmin = {-11.594439904762162, 13.203630051207202, -0.4034394881768612, 0.2367787744557347};

    for (int i : {0,1,2,3})
        CHECK( std::abs(x[i]-globmin[i])<xtol );
}


#include <catch2/catch.hpp>

#define _USE_MATH_DEFINES
#include <cmath>
#include "stlbfgs.h"

using namespace STLBFGS;
static const double xtol = 1e-3;

TEST_CASE("Beale function", "[L-BFGS]") {
    const Optimizer::func_grad_eval fcn = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        double t1 = 1 - x[1];
        double t2 = 1 - x[1]*x[1];
        double t3 = 1 - x[1]*x[1]*x[1];

        double f1 = 1.5   - x[0]*t1;
        double f2 = 2.25  - x[0]*t2;
        double f3 = 2.625 - x[0]*t3;

        f = f1*f1 + f2*f2 + f3*f3;
        g = {
            -2. * (f1*t1 + f2*t2 + f3*t3),
            +2. * x[0] * (f1 + 2*f2*x[1] + 3*f3*x[1]*x[1])
        };
    };

    std::vector<double> x = {1., 1.};
    Optimizer opt{fcn};
    opt.run(x);

    CHECK( std::abs(x[0]-3.0) < xtol );
    CHECK( std::abs(x[1]-0.5) < xtol );
}


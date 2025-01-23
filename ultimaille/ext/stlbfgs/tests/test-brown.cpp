#include <catch2/catch.hpp>

#define _USE_MATH_DEFINES
#include <cmath>
#include "stlbfgs.h"

using namespace STLBFGS;
static const double xtol = 1e-5;

TEST_CASE("Brown's badly scaled function", "[L-BFGS]") {
    const Optimizer::func_grad_eval fcn = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        double f1 = x[0] - 1e6;
        double f2 = x[1] - 2e-6;
        double f3 = static_cast<double>(x[0]*x[1]) - 2; // Prevent fused multiply subtract.
        f = f1*f1 + f2*f2 + f3*f3;
        g = {
            2*f1 + 2*f3*x[1],
            2*f2 + 2*f3*x[0]
        };
    };

    std::vector<double> x = {1., 1.};
    Optimizer opt{fcn};
    opt.run(x);

    double f;
    std::vector<double> g;
    fcn(x, f, g);

    CHECK( std::abs(f-0)<xtol );
    CHECK( std::abs(x[0]-1e+6)<xtol );
    CHECK( std::abs(x[1]-2e-6)<xtol );
}


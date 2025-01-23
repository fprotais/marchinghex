#include <catch2/catch.hpp>

#define _USE_MATH_DEFINES
#include <cmath>
#include "stlbfgs.h"

using namespace STLBFGS;
static const double xtol = 1e-5;

TEST_CASE("Box' three-dimensional function", "[L-BFGS]") {
    const Optimizer::func_grad_eval fcn = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        f = 0;
        g = {0, 0, 0};

        for (int i=1; i<=10; i++) {
            double c = -i/10.;
            double y = std::exp(c) - std::exp(10*c);
            double t = std::exp(c*x[0]) - std::exp(c*x[1]) - x[2]*y;

            f += t*t;
            g[0] +=  2 * t * c * std::exp(c*x[0]);
            g[1] += -2 * t * c * std::exp(c*x[1]);
            g[2] += -2 * t * y;
        }
    };

    std::vector<double> x = {0., 10., 20.};
    Optimizer opt{fcn};
    opt.ftol = 1e-11;
    opt.run(x);

    double f;
    std::vector<double> g;
    fcn(x, f, g);

    CHECK( std::abs(f-0)<xtol );
    // three global minima
    CHECK( (
                (
                 std::abs(x[0] -  1.) < xtol &&
                 std::abs(x[1] - 10.) < xtol &&
                 std::abs(x[2] -  1.) < xtol
                ) ||
                (
                 std::abs(x[0] - 10.) < xtol &&
                 std::abs(x[1] -  1.) < xtol &&
                 std::abs(x[2] +  1.) < xtol
                ) ||
                (
                 std::abs(x[0] - x[1]) < xtol &&
                 std::abs(x[2] -  0.)  < xtol
                )
           )
         );
}


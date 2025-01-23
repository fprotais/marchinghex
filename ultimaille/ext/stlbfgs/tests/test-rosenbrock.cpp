#include <catch2/catch.hpp>
#include <iostream>
#include "stlbfgs.h"

using namespace STLBFGS;

static const double xtol = 1e-3;

TEST_CASE("Rosenbrock function", "[L-BFGS]") {
    std::vector<double> sol = {-1.2, 1.0};

    const Optimizer::func_grad_eval rosenbrock = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        f = (1. - x[0])*(1. - x[0]) + 100.*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]);
        g[0] = 2.*(200.*x[0]*x[0]*x[0] - 200.*x[0]*x[1] + x[0] - 1.);
        g[1] = 200.*(x[1] - x[0]*x[0]);

    };

    Optimizer opt{rosenbrock};
    opt.run(sol);

    CHECK(std::abs(sol[0]-1)<xtol);
    CHECK(std::abs(sol[1]-1)<xtol);
}


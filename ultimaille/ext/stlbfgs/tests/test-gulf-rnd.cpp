#include <catch2/catch.hpp>
#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>
#include "stlbfgs.h"

using namespace STLBFGS;
static const double xtol = 1e-5;

TEST_CASE("Gulf research and development function", "[L-BFGS]") {
    const Optimizer::func_grad_eval fcn = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        f = 0;
        g = {0, 0, 0};

        for (int i=1; i<=99; i++) {
            double arg = i/100.;
            double r = std::pow(-50*std::log(arg), 2.0/3.0) + 25 - x[1];

            double t1 = std::pow(std::abs(r), x[2]) / x[0];
            double t2 = std::exp(-t1);
            double t = t2 - arg;
            double s1 = t1 * t2 * t;

            f += t*t;
            g[0] += s1;
            g[1] += s1 / r;
            g[2] -= s1 * std::log(std::abs(r));
        }
        g[0] *= 2. / x[0];
        g[1] *= 2. * x[2];
        g[2] *= 2.;

    };

    std::vector<double> x = {5., 2.5, 0.15};
    Optimizer opt{fcn};
    opt.ftol = 1e-10;
    opt.run(x);

    double f;
    std::vector<double> g;
    fcn(x, f, g);

    CHECK( std::abs(f-0)<xtol );

    std::vector<double> gnd_truth = {50., 25., 1.5};
    for (int i : {0,1,2})
        CHECK( std::abs(x[i]-gnd_truth[i])<xtol );
}


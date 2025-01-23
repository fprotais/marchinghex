#include <catch2/catch.hpp>
#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>
#include "stlbfgs.h"

using namespace STLBFGS;
static const double xtol = 1e-5;

TEST_CASE("Gaussian function", "[L-BFGS]") {
    const Optimizer::func_grad_eval fcn = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        f = 0;
        g = {0, 0, 0};

        for (int i=1; i<=15; i++) {
            double c = 0.5 * (8-i);
            double b = c - x[2];
            double d = b*b;
            double e = std::exp(-0.5 * x[1] * d);

            double yi = 0;
            switch (i) {
                case 1:
                case 15: yi = 0.0009; break;
                case 2:
                case 14: yi = 0.0044; break;
                case 3:
                case 13: yi = 0.0175; break;
                case 4:
                case 12: yi = 0.0540; break;
                case 5:
                case 11: yi = 0.1295; break;
                case 6:
                case 10: yi = 0.2420; break;
                case 7:
                case 9:  yi = 0.3521; break;
                case 8:  yi = 0.3989;
            }

            double t = x[0]*e - yi;

            f += t*t;
            g[0] += 2 * t * e;
            g[1] -= t * e * d * x[0];
            g[2] += 2 * t * e * x[0] * x[1] * b;
        }
    };

    std::vector<double> x = {.4, 1., 0.};
    Optimizer opt{fcn};
    opt.ftol = 1e-11;
    opt.run(x);

    double f;
    std::vector<double> g;
    fcn(x, f, g);

    CHECK( std::abs(f-1.12793276961912e-08)<xtol );

    std::vector<double> gnd_truth = {0.398956137837997, 1.0000190844805048, 0};

    for (int i : {0,1,2})
        CHECK( std::abs(x[i]-gnd_truth[i])<xtol );
}


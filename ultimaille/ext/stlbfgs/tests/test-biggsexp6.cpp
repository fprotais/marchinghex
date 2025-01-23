#include <catch2/catch.hpp>

#define _USE_MATH_DEFINES
#include <cmath>
#include "stlbfgs.h"

using namespace STLBFGS;
static const double xtol = 1e-5;

TEST_CASE("Biggs EXP6 function", "[L-BFGS]") {
    const Optimizer::func_grad_eval fcn = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        f = 0;
        g = {0, 0, 0, 0, 0, 0};

        for (int i=1; i<=13; i++) {
            double z = i/10.;
            double y = std::exp(-z) - 5*std::exp(-10*z) + 3*std::exp(-4*z);
            double t = x[2]*std::exp(-x[0]*z) - x[3]*std::exp(-x[1]*z) + x[5]*std::exp(-x[4]*z) - y;

            double dfdx0 = -z * x[2] * std::exp(-x[0]*z);
            double dfdx1 =  z * x[3] * std::exp(-x[1]*z);
            double dfdx2 =  std::exp(-x[0] * z);
            double dfdx3 = -std::exp(-x[1] * z);
            double dfdx4 = -z * x[5] * std::exp(-x[4]*z);
            double dfdx5 = std::exp(-x[4] * z);

            f += t*t;
            g[0] += 2*t*dfdx0;
            g[1] += 2*t*dfdx1;
            g[2] += 2*t*dfdx2;
            g[3] += 2*t*dfdx3;
            g[4] += 2*t*dfdx4;
            g[5] += 2*t*dfdx5;
        }
    };

    std::vector<double> x = {1., 2., 1., 1., 1., 1.};
    Optimizer opt{fcn};
    opt.ftol = 1e-12;
    opt.run(x);

    double f;
    std::vector<double> g;
    fcn(x, f, g);

    CHECK( (std::abs(f-5.65565e-3)<xtol || std::abs(f-0)<xtol) );
    if (std::abs(f-0)<xtol) { // global minimum
        CHECK( std::abs(x[0]- 1.) < xtol );
        CHECK( std::abs(x[1]-10.) < xtol );
        CHECK( std::abs(x[2]- 1.) < xtol );
        CHECK( std::abs(x[3]- 5.) < xtol );
        CHECK( std::abs(x[4]- 4.) < xtol );
        CHECK( std::abs(x[5]- 3.) < xtol );
    } else { // local minimum
        CHECK( std::abs(x[0] -  1.711416) < xtol );
        CHECK( std::abs(x[1] - 17.683198) < xtol );
        CHECK( std::abs(x[2] -  1.163144) < xtol );
        CHECK( std::abs(x[3] -  5.186562) < xtol );
        CHECK( std::abs(x[4] -  1.711416) < xtol );
        CHECK( std::abs(x[5] -  1.163144) < xtol );
    }
}


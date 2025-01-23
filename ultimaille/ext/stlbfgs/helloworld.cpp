#include <iostream>
#include "stlbfgs.h"

int main() {
    const STLBFGS::Optimizer::func_grad_eval func = [](const std::vector<double> &x, double &f, std::vector<double> &g) {
        f = (x[0] - 7)*(x[0] - 7) +
            (x[1] - 1)*(x[1] - 1);
        g[0] = 2*(x[0] - 7);
        g[1] = 2*(x[1] - 1);
    };

    STLBFGS::Optimizer opt{func};
    std::vector<double> x = {10, 10};
    opt.run(x);

    std::cout << "Result: x=" << x[0] << ", y=" << x[1] << std::endl;

    if (std::abs(x[0]-7)<1e-3 && std::abs(x[1]-1)<1e-3) {
        std::cout << "Optimization succeeded" << std::endl;
        return 0;
    } else {
        std::cout << "Optimization failed" << std::endl;
        return 1;
    }
}


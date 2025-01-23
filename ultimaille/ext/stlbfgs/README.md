# C++ L-BFGS implementation using plain STL

[![Build](https://github.com/ssloy/stlbfgs/actions/workflows/continuous.yml/badge.svg)](https://github.com/ssloy/stlbfgs/actions/workflows/continuous.yml) [![Nightly](https://github.com/ssloy/stlbfgs/actions/workflows/nightly.yml/badge.svg)](https://github.com/ssloy/stlbfgs/actions/workflows/nightly.yml)

L-BFGS is a quasi-Newton optimization algorithm for solving large nonlinear optimization problems [1,2]. It employs function value and gradient information to search for the local optimum. It uses (as the name suggests) the BGFS (Broyden-Goldfarb-Fletcher-Shanno) algorithm to approximate the inverse Hessian matrix. The size of the memory available to store the approximation of the inverse Hessian is limited (hence L- in the name): in fact, we do not need to store the approximated matrix directly, but rather we need to be able to multiply it by a vector (most often the gradient), and this can be done efficiently by using the history of past updates.

The project has zero external dependencies, no Eigen, nothing, plain standard library. It is a from-scratch implementation and not a wrapper/translation of the original Fortran/Matlab codes by Jorge Nocedal / Dianne O'Leary.

This implementation uses Moré-Thuente line search algorithm [3]. The preconditioning is performed using the M1QN3 strategy [5,6].
Moré-Thuente line search routine is tested against data [3] (Tables 1-6), and L-BFGS routine is tested on problems taken from [4].

# Hello world
Here is a minimal usage example: it shows a minimization of a 2D quadratic function f(x, y) = (x-7)^2 + (y-1)^2, whose minimizer is (x, y) = (7, 1):

```cpp
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

    return std::abs(x[0]-7)>1e-3 || std::abs(x[1]-1)>1e-3;
}

```

# Build and run unit tests:
```sh
git clone https://github.com/ssloy/stlbfgs.git &&
cd stlbfgs &&
mkdir build &&
cd build &&
cmake -DSTLBFGS_UNIT_TESTS:BOOL=ON .. &&
cmake --build . -j &&
cd tests &&
ctest
```

# References

**[1]** J. Nocedal (1980). Updating quasi-Newton matrices with limited storage. Mathematics of Computation, 35/151, 773-782. 3

**[2]** Dong C. Liu, Jorge Nocedal. On the limited memory BFGS method for large scale optimization. Mathematical Programming (1989).

**[3]** Jorge J. Moré, David J. Thuente. Line search algorithms with guaranteed sufficient decrease. ACM Transactions on Mathematical Software  (1994)

**[4]** Jorge J. Moré, Burton S. Garbow, Kenneth E. Hillstrom, "Testing Unconstrained Optimization Software", ACM Transactions on Mathematical Software  (1981)

**[5]** Gilbert JC, Lemaréchal C. The module M1QN3. INRIA Rep., version. 2006;3:21.

**[6]** Gilbert JC, Lemaréchal C. Some numerical experiments with variable-storage quasi-Newton algorithms. Mathematical Programming 45, pp. 407-435, 1989.


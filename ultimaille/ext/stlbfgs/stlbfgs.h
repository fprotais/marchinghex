#ifndef __STLBFGS_H__
#define __STLBFGS_H__

#include <functional>
#include <vector>
#include <deque>
#include "linesearch.h"

#define M1QN3_PRECOND 1
namespace STLBFGS {
    typedef std::vector<double> vector;

    struct Optimizer {
        typedef std::function<void(const vector& x, double &f, vector& g)> func_grad_eval;
        struct IHessian { // L-BFGS approximates inverse Hessian matrix by storing a limited history of past updates
            void mult(const vector &g, vector &result) const; // matrix-vector multiplication
            void add_correction(const vector&s, const vector& y);

            int history_depth = 10;
            typedef std::deque<vector> history;
            history S = {};
            history Y = {};
#if M1QN3_PRECOND
            vector diag = {};
#else
            double gamma = 1.; // TODO remove non-preconditioned version
#endif
        };

        bool run(vector& sol); // actual optimization loop

        const func_grad_eval func_grad;
        IHessian invH = {};  // current inverse Hessian approximation

        // L-BFGS user parameters
        int maxiter = 10000; // maximum number of quasi-Newton updates
        double ftol = 1e-6;  // the iteration stops when (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= ftol
        double gtol = 1e-14; // the iteration stops when ||g||/max(1,||x||) <= gtol
        double gmax = 1e-14; // the iteration stops when max{|g_i|, i = 1, ..., n} <= gmax

        // Line search user parameters: the step size must satisfy Wolfe conditions with these parameters
        double mu  = 1e-4; // sufficient decrease constant (Armijo rule)
        double eta = 9e-1; // curvature condition constant
//      int lsmaxfev = 16;  // TODO move all line search parameters here

        bool verbose = true;
    };
}

#endif //__STLBFGS_H__


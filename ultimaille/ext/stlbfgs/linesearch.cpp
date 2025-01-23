#include <tuple>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>

#include "linesearch.h"

namespace STLBFGS {
    bool sufficient_decrease(Sample phi0, Sample phia, double mu) {
        return phia.f <= phi0.f + mu*phia.a*phi0.d;
    }

    bool curvature_condition(Sample phi0, Sample phia, double eta) {
        return std::abs(phia.d) <= eta*std::abs(phi0.d);
    }

    // minimizer of the cubic function that interpolates f(a), f'(a), f(b), f'(b) within the given interval
    double find_cubic_minimizer(double a, double fa, double ga, double b, double fb, double gb) {
        if (a>b) {
            std::swap( a,  b);
            std::swap(fa, fb);
            std::swap(ga, gb);
        }
        double z = 3.*(fa - fb)/(b - a) + ga + gb;
        double D = z*z - ga*gb;
        if (D<=0) return std::numeric_limits<double>::max(); // no minumum in the interval, +inf here because of the linesearch nature
        double w = std::sqrt(D); // this code assumes a<b; negate this value if b<a.
        return b - ((b - a)*(gb + w - z))/(gb - ga + 2.*w);
    }

    // minimizer of the quadratic function that interpolates f(a), f'(a), f(b) within the given interval
    double find_quadratic_minimizer(double a, double fa, double ga, double b, double fb) {
        return a + (((b - a)*(b - a))*ga)/(2.*(fa - fb + (b - a)*ga));
    }

    // minimizer of the quadratic function that interpolates f'(a), f'(b) within the given interval
    // N.B. the function itself is undetermined since we miss information like f(a) or f(b); however the minimizer is well-defined
    double find_quadratic_minimizer(double a, double ga, double b, double gb) {
        return b + ((b - a)*gb)/(ga - gb);
    }

    std::tuple<double,int> trial_value(const Sample &l, const Sample &t, const Sample &u, const bool bracketed) {
        double ac = find_cubic_minimizer(l.a, l.f, l.d, t.a, t.f, t.d);
        if (t.f > l.f) { // Case 1: a higher function value. The minimum is bracketed.
            double aq = find_quadratic_minimizer(l.a, l.f, l.d, t.a, t.f);
            double res = (std::abs(ac - l.a) < std::abs(aq - l.a)) ? ac : (aq + ac)/2.;
            return std::make_tuple(res, 1);
        }

        double as = find_quadratic_minimizer(l.a, l.d, t.a, t.d);
        if ((l.d>0 && t.d<0) || (l.d<0 && t.d>0)) { // Case 2: A lower function value and derivatives of opposite sign. The minimum is bracketed.
            double res = (std::abs(ac - t.a) >= std::abs(as - t.a)) ? ac : as;
            return std::make_tuple(res, 2);
        }

        if (std::abs(t.d) <= std::abs(l.d)) { // Case 3: A lower function value, derivatives of the same sign, and the magnitude of the derivative decreases.
            if ( // the cubic function may not have a minimizer; moreover, even if it exists, it can be in the wrong direction; fix it
                    (l.a<u.a && ac<=t.a) ||
                    (l.a>u.a && ac>=t.a)
               ) ac = u.a;
            double res = bracketed ?
                (std::abs(ac - t.a) < std::abs(as - t.a) ? ac : as) :
                (std::abs(ac - t.a) > std::abs(as - t.a) ? ac : as);
            return std::make_tuple(res, 3);
        }

        // Case 4: A lower function value, derivatives of the same sign, and the magnitude of the derivative does not decrease.
        double res = bracketed ? find_cubic_minimizer(t.a, t.f, t.d, u.a, u.f, u.d) : u.a;
        return std::make_tuple(res, 4);
    }

    bool line_search_more_thuente(const linesearch_function phi, const Sample phi0, double &at, const double mu, const double eta, const int lsmaxfev) {
        double init_step = at;
        bool stage1 = true;  // use function psi instead if phi
        bool bracketed = false;
        Sample phil = phi0;
        Sample phiu = { 0, 0, 0 }; // N.B: .f and .d members are not used until the minimum is bracketed
        double width_prev = std::numeric_limits<double>::max();

        for (int nfev=0; nfev<lsmaxfev; nfev++) {
            if (!bracketed)
                phiu.a = at + 4.*(at - phil.a);
            Sample phit = phi(at);

            if (sufficient_decrease(phi0, phit, mu) && curvature_condition(phi0, phit, eta))
                return true;

            // TODO error handling
            if (
                    !std::isfinite(phil.a) || !std::isfinite(phil.f) || !std::isfinite(phil.d) ||
                    !std::isfinite(phit.a) || !std::isfinite(phit.f) || !std::isfinite(phit.d) ||
                    !std::isfinite(phiu.a) || !std::isfinite(phiu.f) || !std::isfinite(phiu.d) ||
                    (
                     (phil.a>=phiu.a || phil.a>=phit.a || phit.a>=phiu.a || phil.d>=0) &&
                     (phil.a<=phiu.a || phil.a<=phit.a || phit.a<=phiu.a || phil.d<=0)
                    )
               ) break;

            auto psi = [phi0, mu](const Sample &phia) {
                return Sample{ phia.a, phia.f - phi0.f - mu*phi0.d*phia.a, phia.d - mu*phi0.d };
            };

            // Decide if we want to switch to using a "Modified Updating Algorithm" (shown after theorem 3.2 in the paper)
            // by switching from using function psi to using function phi.
            // The decision follows the logic in the paragraph right before theorem 3.3 in the paper.
            stage1 = stage1 && (psi(phit).f>0 /*|| phit.d<=0*/ || phit.d < std::min(mu, eta)*phi0.d);
            int caseno;
            std::tie(at, caseno) = stage1 && phit.f<= phil.f && psi(phit).f>0 ?
                trial_value(psi(phil), psi(phit), psi(phiu), bracketed) :
                trial_value(    phil,      phit,      phiu,  bracketed);

            bracketed = bracketed || (caseno==1 || caseno==2);
            double width = std::abs(phiu.a - phil.a);

            // update the interval of uncertainty; note that the update does not depend on the new trial value
            if (1==caseno)
                phiu = phit;
            else if (2==caseno) {
                phiu = phil;
                phil = phit;
            } else
                phil = phit;

            if (bracketed) {
                if (caseno==1 || caseno==3) {
                    if (std::abs(phiu.a-phil.a) >= .66*width_prev) // force a sufficient decrease in the size of the interval of uncertainty.
                        at = (phil.a+phiu.a)/2.;
                    else { // safeguard the trial value
                        double safeguard1 = phil.a + .66*(phiu.a-phil.a); // the magic constant is used in the Moré-Thuente paper (Section 4, Case 3).
                        at = phil.a<phiu.a ?
                            std::min(safeguard1, at) :
                            std::max(safeguard1, at);
                        double safeguard2 = phil.a + .001*(phiu.a-phil.a);
                        at = phil.a>phiu.a ?
                            std::min(safeguard2, at) :
                            std::max(safeguard2, at);
                    }
                }
                width_prev = width;
            }

            at = phil.a<phiu.a ? // force the step to be within the interval bounds
                std::max(phil.a, std::min(phiu.a, at)) :
                std::min(phil.a, std::max(phiu.a, at));
        }

        at = init_step;
        return false;
    }

    // a version of line search by Nocedal and Wright, Numerical optimization (2006)
    bool line_search_backtracking(const linesearch_function phi, const Sample phi0, double &at, const double mu, const double eta, const int lsmaxfev) {
        double init_step = at;
        Sample phil = phi0;
        Sample phiu = { 0, 0, 0 };

        int nfev = 0;
        for (; nfev<lsmaxfev; nfev++) { // bracketing phase, see Numerical optimization, Algorithm 3.5 (Line Search Algorithm)
            phiu = phi(at);
            // TODO error handling
            if (!sufficient_decrease(phi0, phiu, mu) || (nfev && phiu.f>=phil.f))
                break;
            if (curvature_condition(phi0, phiu, eta))
                return true;
            if (phiu.d>=0) {
                std::swap(phil, phiu);
                break;
            }
            at += 2.*(at - phil.a);
            phil = phiu;
        }

        for (; nfev<lsmaxfev; nfev++) { // zoom phase (bracketing succeeded), see Numerical optimization, Algorithm 3.6 (Zoom)
            if (
                    (phil.d<0 || phiu.a>=phil.a) &&
                    (phil.d>0 || phiu.a<=phil.a)
               ) break;
            at = (phil.a+phiu.a)/2.;
            Sample phit = phi(at);
            if (!sufficient_decrease(phi0, phit, mu) || phit.f>=phil.f)
                phiu = phit;
            else {
                if (curvature_condition(phi0, phit, eta))
                    return true;
                if (
                        (phit.d>=0 && phiu.a>phil.a) ||
                        (phit.d<=0 && phiu.a<phil.a)
                   )
                    phiu = phil;
                phil = phit;
            }
        }

        at = init_step;
        return false; // zoom failed
    }
}


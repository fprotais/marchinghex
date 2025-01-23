#define _USE_MATH_DEFINES
#include <cmath>
#include <catch2/catch.hpp>
#include <iostream>

#include "linesearch.h"

using namespace STLBFGS;

namespace STLBFGS {
    double find_cubic_minimizer(double a, double fa, double ga, double b, double fb, double gb);
    double find_quadratic_minimizer(double a, double fa, double ga, double b, double fb);
    double find_quadratic_minimizer(double a, double ga, double b, double gb);
}

TEST_CASE("cubic 1", "[interpolation]") {
    auto f = [](const double a) { return a*a*a - a; };
    auto g = [](const double a) { return 3*a*a - 1; };
    double a = -.4;
    double b = 19.;
    double minimizer = find_cubic_minimizer(a, f(a), g(a), b, f(b), g(b));
    REQUIRE( std::abs(minimizer-1./std::sqrt(3.))<1e-14 );
}

TEST_CASE("cubic 2", "[interpolation]") {
    auto f = [](const double a) { return a*a*a - 3*a*a - 144*a + 432; };
    auto g = [](const double a) { return 3*a*a - 6*a - 144; };
    double a = -5.;
    double b = 113.;
    double minimizer = find_cubic_minimizer(a, f(a), g(a), b, f(b), g(b));
    REQUIRE( std::abs(minimizer-8.)<1e-14 );
}

TEST_CASE("quadratic", "[interpolation]") {
    for (int i=0; i<8; i++) {
        double A = rand()/(double)RAND_MAX + .2134;
        double B = rand()/(double)RAND_MAX;
        double C = rand()/(double)RAND_MAX;
        double M = -B/(2.*A);
        auto f = [&](const double a) { return A*a*a + B*a + C; };
        auto g = [&](const double a) { return 2*A*a + B; };
        double a = M - rand()/(double)RAND_MAX;
        double b = M + rand()/(double)RAND_MAX;
        double minimizer1 = find_quadratic_minimizer(a, f(a), g(a), b, f(b));
        REQUIRE( std::abs(minimizer1 - M)<1e-14 );
        double minimizer2 = find_quadratic_minimizer(a, g(a), b, g(b));
        REQUIRE( std::abs(minimizer2 - M)<1e-14 );
    }
}

template <typename T> auto square(const T &number) { return number * number; }

TEST_CASE("Table 1", "[Mor\\'e-Thuente]") {
    int nfev = 0;
    const linesearch_function func = [&nfev](const double alpha) -> Sample {
        nfev++;
        constexpr double beta = 2.;
        double f = -alpha/(alpha*alpha + beta);
        double g = (square(alpha)-beta)/square(beta+square(alpha));
        return { alpha, f, g };
    };
    double alpha = 1e-3;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-1) );
    CHECK( std::abs(alpha-1.3650)< 1e-4 );
    CHECK( nfev == 6+1 );
    alpha = 1e-1;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-1) );
    CHECK( std::abs(alpha-1.4414)< 1e-4 );
    CHECK( nfev == 3+1 );
    alpha = 1e+1;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-1) );
    CHECK( std::abs(alpha-10)< 1e-4 );
    CHECK( nfev == 1+1 );
    alpha = 1e+3;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-1) );
    CHECK( std::abs(alpha-36.8876)< 1e-4 );
    CHECK( nfev == 4+1 );
}

TEST_CASE("Table 2", "[Mor\\'e-Thuente]") {
    int nfev = 0;
    const linesearch_function func = [&nfev](const double alpha) -> Sample {
        nfev++;
        constexpr double beta = .004;
        double f = pow(alpha+beta, 5.) - 2.*pow(alpha+beta, 4.);
        double g = 5.*pow(alpha+beta, 4.) - 8.*pow(alpha+beta, 3.);
        return { alpha, f, g };
    };
    double alpha = 1e-3;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-1, 1e-1) );
    CHECK( std::abs(alpha-1.5960)< 1e-4 );
    CHECK( nfev == 12+1 );
    alpha = 1e-1;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-1, 1e-1) );
    CHECK( std::abs(alpha-1.5960)< 1e-4 );
    CHECK( nfev == 8+1 );
    alpha = 1e+1;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-1, 1e-1) );
    CHECK( std::abs(alpha-1.5960)< 1e-4 );
    CHECK( nfev == 8+1 );
    alpha = 1e+3;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-1, 1e-1) );
    CHECK( std::abs(alpha-1.5960)< 1e-4 );
    CHECK( nfev == 11+1 );
}

TEST_CASE("Table 3", "[Mor\\'e-Thuente]") {
    int nfev = 0;
    const linesearch_function func = [&nfev](const double alpha) -> Sample {
        nfev++;
        constexpr double beta = .01;
        constexpr double l = 39;

        double f = 2.*(1.-beta)/(l*M_PI)*sin(l*M_PI_2*alpha);
        double g = (1.-beta)*cos(l*M_PI_2*alpha);
        if (alpha<=1.-beta) {
            f += 1.-alpha;
            g += -1.;
        } else if (alpha>=1+beta) {
            f += alpha-1.;
            g += 1.;
        } else {
            f += square(alpha-1.)/(2.*beta) + beta/2.;
            g += (alpha-1.)/beta;
        }

        return { alpha, f, g };
    };
    double alpha = 1e-3;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-1, 1e-1) );
    CHECK( std::abs(alpha-1.)< 1e-4 );
    CHECK( nfev == 12+1 );
    alpha = 1e-1;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-1, 1e-1) );
    CHECK( std::abs(alpha-1.)< 1e-4 );
    CHECK( nfev == 12+1 );
    alpha = 1e+1;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-1, 1e-1) );
    CHECK( std::abs(alpha-1.)< 1e-4 );
    CHECK( nfev == 10+1 );
    alpha = 1e+3;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-1, 1e-1) );
    CHECK( std::abs(alpha-1.)< 1e-4 );
    CHECK( nfev == 13+1 );
}

double gamma(const double beta) {
    return std::sqrt(1.+square(beta)) - beta;
}

double phi(const double beta1, const double beta2, const double alpha) {
    return gamma(beta1)*std::sqrt(square(1.-alpha) + square(beta2)) + gamma(beta2)*std::sqrt(square(alpha) + square(beta1));
}

double phi_diff(const double beta1, const double beta2, const double alpha) {
    return gamma(beta1)*(alpha-1.)/std::sqrt(square(1.-alpha) + square(beta2)) + gamma(beta2)*alpha/std::sqrt(square(alpha) + square(beta1));
}

TEST_CASE("Table 4", "[Mor\\'e-Thuente]") {
    int nfev = 0;
    const linesearch_function func = [&nfev](const double alpha) -> Sample {
        nfev++;
        constexpr double beta1 = 1e-3;
        constexpr double beta2 = 1e-3;
        return { alpha, phi(beta1, beta2, alpha), phi_diff(beta1, beta2, alpha) };
    };

    double alpha = 1e-3;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-3) );
    CHECK( std::abs(alpha-0.0850)< 1e-4 );
    CHECK( nfev == 4+1 );

    alpha = 1e-1;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-3) );
    CHECK( std::abs(alpha-0.10)< 1e-4 );
    CHECK( nfev == 1+1 );

    alpha = 1e+1;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-3) );
    CHECK( std::abs(alpha-0.3491)< 1e-4 );
    CHECK( nfev == 3+1 );

    alpha = 1e+3;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-3) );
    CHECK( std::abs(alpha-0.8294)< 1e-4 );
    CHECK( nfev == 4+1 );

    // The tests from the paper do not enter the code path where the function is modified much. This test does run that part.

    alpha = 1.;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 0.1, 0.9) );
    CHECK( std::abs(alpha-0.0039)< 1e-4 );
    CHECK( nfev == 6+1 );
}

TEST_CASE("Table 5", "[Mor\\'e-Thuente]") {
    int nfev = 0;
    const linesearch_function func = [&nfev](const double alpha) -> Sample {
        nfev++;
        constexpr double beta1 = 1e-2;
        constexpr double beta2 = 1e-3;
        return { alpha, phi(beta1, beta2, alpha), phi_diff(beta1, beta2, alpha) };
    };
    double alpha = 1e-3;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-3) );
    CHECK( std::abs(alpha-0.0750)< 1e-4 );
    CHECK( nfev == 6+1 );

    alpha = 1e-1;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-3) );
    CHECK( std::abs(alpha-0.0775)< 1e-4 );
    CHECK( nfev == 3+1 );

    alpha = 1e+1;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-3) );
    CHECK( std::abs(alpha-0.0731)< 1e-4 );
    CHECK( nfev == 7+1 );

    alpha = 1e+3;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-3) );
    CHECK( std::abs(alpha-0.0762)< 1e-4 );
    CHECK( nfev == 8+1 );

    // The tests from the paper do not enter the code path where the function is modified much. This test does run that part.

    alpha = 1.;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 0.1, 0.9) );
    CHECK( std::abs(alpha-0.0434)< 1e-4 );
    CHECK( nfev == 4+1 );
}

TEST_CASE("Table 6", "[Mor\\'e-Thuente]") {
    int nfev = 0;
    const linesearch_function func = [&nfev](const double alpha) -> Sample {
        nfev++;
        constexpr double beta1 = 1e-3;
        constexpr double beta2 = 1e-2;
        return { alpha, phi(beta1, beta2, alpha), phi_diff(beta1, beta2, alpha) };
    };
    double alpha = 1e-3;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-3) );
    CHECK( std::abs(alpha-0.9279)< 1e-4 );
    CHECK( nfev == 13+1 );

    alpha = 1e-1;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-3) );
    CHECK( std::abs(alpha-0.9262)< 1e-4 );
    CHECK( nfev == 11+1 );

    alpha = 1e+1;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-3) );
    CHECK( std::abs(alpha-0.9248)< 1e-4 );
    CHECK( nfev == 8+1 );

    alpha = 1e+3;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 1e-3, 1e-3) );

// low safeguarding affects this
#if 0
    CHECK( std::abs(alpha-0.9244)< 1e-4 );
    CHECK( nfev == 11+1 );
#else
    CHECK( std::abs(alpha-0.9250)< 1e-4 );
    CHECK( nfev == 10+1 );
#endif
    // The tests from the paper do not enter the code path where the function is modified much. This test does run that part.

    alpha = 1.;
    nfev = 0;
    CHECK( line_search_more_thuente(func, func(0), alpha, 0.1, 0.9) );
    CHECK( std::abs(alpha-0.0040)< 1e-4 );
    CHECK( nfev == 6+1 );
}


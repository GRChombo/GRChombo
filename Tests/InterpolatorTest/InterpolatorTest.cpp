/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to LICENSE, in Chombo's root directory.
 */
#endif

// Chombo includes
#include "parstream.H" //Gives us pout()

// General includes:
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>

#include "Lagrange.hpp"
#include "QuinticConvolution.hpp"
#include "SetupFunctions.hpp"
#include "SimpleArrayBox.hpp"
#include "SimpleInterpSource.hpp"

// Chombo namespace
#include "UsingNamespace.H"

// cell centered grid
// double get_dx(double L, int num_points_u) { return L / num_points_u; }
// double get_x(int i, double m_dx) { return m_dx * (i + 0.5); }
// double get_idx(double x, double m_dx) { return x / m_dx - 0.5; }

// node centered grid
double get_dx(double L, int num_points_u) { return L / (num_points_u - 1.); }
double get_x(int i, double m_dx) { return m_dx * i; }
double get_idx(double x, double m_dx) { return x / m_dx; }

double func(double x, double L)
{
    // cubic
    // maximum is 1 at x = L / 6 * ( 3 +- sqrt(3) )
    // return x * (x - L / 2.) * (x - L) / (L * L * L) * 12. * sqrt(3.);

    // quadratic
    // maximum is 1 at x = L / 2
    // return x * (x - L) / (L * L) * 4.;

    // non-linear
    return sin(2. * M_PI * x / L);
}

int runInterpolatorTest(int argc, char *argv[])
{

    const int num_points_u = 200;
    const double L = 10.;
    const double m_dx = get_dx(L, num_points_u);

    std::vector<double> f(num_points_u);

    for (int i = 0; i < num_points_u; ++i)
    {
        double x = get_x(i, m_dx);
        f[i] = func(x, L);
    }

    // Note that the 2nd argument of .setup 'dx' does not matter for 0th order
    // interpolation -> set to 0
    SimpleInterpSource<1> source({num_points_u}, {0.});
    SimpleArrayBox<1> box({num_points_u}, f);

    double test_point =
        M_PI; // just because it's irrational, hence for sure not in the grid
    double test_index = get_idx(test_point, m_dx);

    bool verbosity = true;
    Lagrange<4, 1> interpolator1(source, verbosity);
    interpolator1.setup({0}, {test_index});
    double val1 = interpolator1.interpData(box);

    QuinticConvolution<1> interpolator2(source, verbosity);
    interpolator2.setup({0}, {test_index});
    double val2 = interpolator2.interpData(box);

    double exact = func(test_point, L);

    double err1 = abs(val1 / exact - 1.) * 100.;
    double err2 = abs(val2 / exact - 1.) * 100.;

    pout() << std::setprecision(9);
    pout() << "Exact value is: " << exact << std::endl;
    pout() << "Lagrange Interpolator gave: " << val1 << " (" << err1
           << "% error)" << std::endl;
    pout() << "Quintic Convolution Interpolator gave: " << val2 << " (" << err2
           << "% error)" << std::endl;

    bool wrong = (err1 + err2 > 1e-5);

    return wrong;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runInterpolatorTest(argc, argv);

    if (status == 0)
        std::cout << "Interpolator test passed." << endl;
    else
        std::cout << "Interpolator test failed with return code " << status
                  << endl;

    mainFinalize();
    return status;
}

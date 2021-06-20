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
#include <iomanip>
#include <iostream>
#include <sys/time.h>

using std::endl;
#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"

#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "AMRInterpolator.hpp"
#include "Lagrange.hpp"
#include "SphericalExtraction.hpp"
#include "SphericalExtractionTestLevel.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

// Chombo namespace
#include "UsingNamespace.H"

int runSphericalExtractionTest(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    std::string in_string = argv[argc - 1];
    pout() << in_string << std::endl;
    char const *in_file = argv[argc - 1];
    GRParmParse pp(0, argv + argc, NULL, in_file);
    SimulationParameters sim_params(pp);

    GRAMR gr_amr;
    DefaultLevelFactory<SphericalExtractionTestLevel>
        surface_extraction_test_level_fact(gr_amr, sim_params);
    // the initial data for the two variables is the spherical harmonic
    // specified by params
    setupAMRObject(gr_amr, surface_extraction_test_level_fact);

    AMRInterpolator<Lagrange<4>> interpolator(
        gr_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params);

    // low resolution spherical extraction
    SphericalExtraction spherical_extraction_lo(
        sim_params.extraction_params_lo,
        sim_params.coarsest_dx * sim_params.dt_multiplier, 0.0, true, 0.0);
    spherical_extraction_lo.add_var(c_phi_Re);
    spherical_extraction_lo.add_var(c_phi_Im);
    spherical_extraction_lo.extract(&interpolator);
    spherical_extraction_lo.write_extraction("ExtractionOutLo_");

    // high resolution spherical extraction
    spherical_extraction_params_t extraction_params_hi =
        sim_params.extraction_params_lo;
    // we are only checking the converence in theta integration
    // extraction_params_hi.num_points_phi *= 2;
    extraction_params_hi.num_points_theta *= 2;
    // need to subtract a point as it's the number of subintervals we want to
    // double for theta
    extraction_params_hi.num_points_theta -= 1;
    SphericalExtraction spherical_extraction_hi(
        extraction_params_hi, sim_params.coarsest_dx * sim_params.dt_multiplier,
        0.0, true, 0.0);
    spherical_extraction_hi.add_var(c_phi_Re);
    spherical_extraction_hi.add_var(c_phi_Im);
    spherical_extraction_hi.extract(&interpolator);
    spherical_extraction_hi.write_extraction("ExtractionOutHi_");

    // real part is the zeroth componenent and imaginary part is first component
    auto extracted_harmonic = [](std::vector<double> &data, double, double,
                                 double) {
        return std::make_pair(data[0], data[1]);
    };

    // add the spherical harmonic mode integrands for each resolution and for
    // the trapezium rule, Simpson's rule, Simpson's 3/8 rule and Boole's rule
    // Always use trapezium rule in phi as this is periodic
    bool broadcast_integral = true;
    std::pair<std::vector<double>, std::vector<double>> integral_lo_trapezium,
        integral_hi_trapezium;
    spherical_extraction_lo.add_mode_integrand(
        sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
        integral_lo_trapezium, IntegrationMethod::trapezium,
        IntegrationMethod::simpson, broadcast_integral);
    spherical_extraction_hi.add_mode_integrand(
        sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
        integral_hi_trapezium, IntegrationMethod::trapezium,
        IntegrationMethod::simpson, broadcast_integral);
    std::pair<std::vector<double>, std::vector<double>> integral_lo_simpson,
        integral_hi_simpson;
    spherical_extraction_lo.add_mode_integrand(
        sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
        integral_lo_simpson, IntegrationMethod::simpson,
        IntegrationMethod::simpson, broadcast_integral);
    spherical_extraction_hi.add_mode_integrand(
        sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
        integral_hi_simpson, IntegrationMethod::simpson,
        IntegrationMethod::simpson, broadcast_integral);
    std::pair<std::vector<double>, std::vector<double>> integral_lo_simpson38,
        integral_hi_simpson38;
    spherical_extraction_lo.add_mode_integrand(
        sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
        integral_lo_simpson38, IntegrationMethod::simpson38,
        IntegrationMethod::simpson, broadcast_integral);
    spherical_extraction_hi.add_mode_integrand(
        sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
        integral_hi_simpson38, IntegrationMethod::simpson38,
        IntegrationMethod::simpson, broadcast_integral);
    std::pair<std::vector<double>, std::vector<double>> integral_lo_boole,
        integral_hi_boole;
    spherical_extraction_lo.add_mode_integrand(
        sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
        integral_lo_boole, IntegrationMethod::boole, IntegrationMethod::simpson,
        broadcast_integral);
    spherical_extraction_hi.add_mode_integrand(
        sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
        integral_hi_boole, IntegrationMethod::boole, IntegrationMethod::simpson,
        broadcast_integral);

    // do the surface integration
    spherical_extraction_lo.integrate();
    spherical_extraction_hi.integrate();

    int status = 0;
    pout() << std::setprecision(10);

    for (int iradius = 0;
         iradius < sim_params.extraction_params_lo.num_extraction_radii;
         ++iradius)
    {
        double r = sim_params.extraction_params_lo.extraction_radii[iradius];
        double integral_re_lo_trapezium =
            (integral_lo_trapezium.first)[iradius];
        double integral_re_hi_trapezium =
            (integral_hi_trapezium.first)[iradius];
        double integral_re_lo_simpson = (integral_lo_simpson.first)[iradius];
        double integral_re_hi_simpson = (integral_hi_simpson.first)[iradius];
        double integral_re_lo_simpson38 =
            (integral_lo_simpson38.first)[iradius];
        double integral_re_hi_simpson38 =
            (integral_hi_simpson38.first)[iradius];
        double integral_re_lo_boole = (integral_lo_boole.first)[iradius];
        double integral_re_hi_boole = (integral_hi_boole.first)[iradius];
        double analytic_integral = 1.0;

        double convergence_factor_trapezium =
            std::abs((integral_re_lo_trapezium - analytic_integral) /
                     (integral_re_hi_trapezium - analytic_integral));
        double convergence_factor_simpson =
            std::abs((integral_re_lo_simpson - analytic_integral) /
                     (integral_re_hi_simpson - analytic_integral));
        double convergence_factor_simpson38 =
            std::abs((integral_re_lo_simpson38 - analytic_integral) /
                     (integral_re_hi_simpson38 - analytic_integral));
        double convergence_factor_boole =
            std::abs((integral_re_lo_boole - analytic_integral) /
                     (integral_re_hi_boole - analytic_integral));

        double convergence_order_trapezium =
            std::log2(convergence_factor_trapezium);
        double convergence_order_simpson =
            std::log2(convergence_factor_simpson);
        double convergence_order_simpson38 =
            std::log2(convergence_factor_simpson38);
        double convergence_order_boole = std::log2(convergence_factor_boole);

        // trapezium rule should have second order convergence
        status |= (convergence_order_trapezium < 1.5);
        // Simpson's rule should have fourth order convergence
        status |= (convergence_order_simpson < 3.5);
        // Simpson 3/8's rule should have fourth order convergence
        status |= (convergence_order_simpson38 < 3.5);
        // Boole's rule should have sixth order convergence
        status |= (convergence_order_boole < 5.5);

        pout() << "At r = " << r << ":\n";
        pout() << "analytic_integral = " << analytic_integral << "\n";
        pout() << "integral_re_lo_trapezium = " << integral_re_lo_trapezium
               << "\n";
        pout() << "integral_re_hi_trapezium = " << integral_re_hi_trapezium
               << "\n";
        pout() << "integral_re_lo_simpson = " << integral_re_lo_simpson << "\n";
        pout() << "integral_re_hi_simpson = " << integral_re_hi_simpson << "\n";
        pout() << "integral_re_lo_simpson38 = " << integral_re_lo_simpson38
               << "\n";
        pout() << "integral_re_hi_simpson38 = " << integral_re_hi_simpson38
               << "\n";
        pout() << "integral_re_lo_boole = " << integral_re_lo_boole << "\n";
        pout() << "integral_re_hi_boole = " << integral_re_hi_boole << "\n";
        pout() << "convergence_order_trapezium = "
               << convergence_order_trapezium << "\n";
        pout() << "convergence_order_simpson = " << convergence_order_simpson
               << "\n";
        pout() << "convergence_order_simpson38 = "
               << convergence_order_simpson38 << "\n";
        pout() << "convergence_order_boole = " << convergence_order_boole
               << "\n"
               << endl;
    }

    return status;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runSphericalExtractionTest(argc, argv);

    if (status == 0)
        pout() << "SphericalExtractionTest test passed." << endl;
    else
        pout() << "SphericalExtractionTest test failed with return code "
               << status << endl;

    mainFinalize();
    return status;
}

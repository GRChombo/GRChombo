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
#include "SphericalExtractionUniformTestLevel.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

// Chombo namespace
#include "UsingNamespace.H"

void surface_integration(
    const spherical_extraction_params_t &extraction_params_lo,
    const SimulationParameters &sim_params,
    AMRInterpolator<Lagrange<4>> &interpolator,
    std::pair<std::vector<double>, std::vector<double>> &integral_lo,
    std::pair<std::vector<double>, std::vector<double>> &integral_hi,
    const IntegrationMethod &method)
{
    // low resolution spherical extraction
    SphericalExtractionUniform spherical_extraction_lo(
        extraction_params_lo, sim_params.coarsest_dx * sim_params.dt_multiplier,
        0.0, true, 0.0);
    spherical_extraction_lo.add_var(c_phi_Re);
    spherical_extraction_lo.add_var(c_phi_Im);
    spherical_extraction_lo.extract(&interpolator);
    spherical_extraction_lo.write_extraction("ExtractionOutLo_");

    // high resolution spherical extraction
    spherical_extraction_params_t extraction_params_hi = extraction_params_lo;
    // we are only checking the converence in theta integration
    // extraction_params_hi.num_points_phi *= 2;
    extraction_params_hi.num_points_theta *= 2;
    SphericalExtractionUniform spherical_extraction_hi(
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

    // add the spherical harmonic mode integrands for each resolution
    // Always use trapezium rule in phi as this is periodic
    bool broadcast_integral = true;
    spherical_extraction_lo.add_mode_integrand(
        sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
        integral_lo, method, IntegrationMethod::trapezium, broadcast_integral);
    spherical_extraction_hi.add_mode_integrand(
        sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
        integral_hi, method, IntegrationMethod::trapezium, broadcast_integral);

    // do the surface integration
    spherical_extraction_lo.integrate();
    spherical_extraction_hi.integrate();
}

int runSphericalExtractionUniformTest(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    std::string in_string = argv[argc - 1];
    pout() << in_string << std::endl;
    char const *in_file = argv[argc - 1];
    GRParmParse pp(0, argv + argc, NULL, in_file);
    SimulationParameters sim_params(pp);

    GRAMR gr_amr;
    DefaultLevelFactory<SphericalExtractionUniformTestLevel>
        surface_extraction_test_level_fact(gr_amr, sim_params);
    // the initial data for the two variables is the spherical harmonic
    // specified by params
    setupAMRObject(gr_amr, surface_extraction_test_level_fact);

    AMRInterpolator<Lagrange<4>> interpolator(
        gr_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params);

    std::pair<std::vector<double>, std::vector<double>> integral_lo_midpoint,
        integral_hi_midpoint;
    spherical_extraction_params_t extraction_params_lo_midpoint =
        sim_params.extraction_params_lo;
    // must be multiple of 2
    extraction_params_lo_midpoint.num_points_theta =
        std::ceil(extraction_params_lo_midpoint.num_points_theta / 2) * 2;
    surface_integration(extraction_params_lo_midpoint, sim_params, interpolator,
                        integral_lo_midpoint, integral_hi_midpoint,
                        IntegrationMethod::midpoint);

    std::pair<std::vector<double>, std::vector<double>>
        integral_lo_open_2nd_order, integral_hi_open_2nd_order;
    spherical_extraction_params_t extraction_params_lo_open_2nd_order =
        sim_params.extraction_params_lo;
    // must be multiple of 3
    extraction_params_lo_open_2nd_order.num_points_theta =
        std::ceil(extraction_params_lo_open_2nd_order.num_points_theta / 3) * 3;
    surface_integration(extraction_params_lo_open_2nd_order, sim_params,
                        interpolator, integral_lo_open_2nd_order,
                        integral_hi_open_2nd_order,
                        IntegrationMethod::open_2nd_order);

    std::pair<std::vector<double>, std::vector<double>>
        integral_lo_open_3rd_order, integral_hi_open_3rd_order;
    spherical_extraction_params_t extraction_params_lo_open_3rd_order =
        sim_params.extraction_params_lo;
    // must be multiple of 4
    extraction_params_lo_open_3rd_order.num_points_theta =
        std::ceil(extraction_params_lo_open_3rd_order.num_points_theta / 4) * 4;
    surface_integration(extraction_params_lo_open_3rd_order, sim_params,
                        interpolator, integral_lo_open_3rd_order,
                        integral_hi_open_3rd_order,
                        IntegrationMethod::open_3rd_order);

    std::pair<std::vector<double>, std::vector<double>>
        integral_lo_open_4th_order, integral_hi_open_4th_order;
    spherical_extraction_params_t extraction_params_lo_open_4th_order =
        sim_params.extraction_params_lo;
    // must be multiple of 5
    extraction_params_lo_open_4th_order.num_points_theta =
        std::floor(extraction_params_lo_open_4th_order.num_points_theta / 5) *
        5;
    surface_integration(extraction_params_lo_open_4th_order, sim_params,
                        interpolator, integral_lo_open_4th_order,
                        integral_hi_open_4th_order,
                        IntegrationMethod::open_4th_order);

    int status = 0;
    pout() << std::setprecision(10);

    for (int iradius = 0;
         iradius < sim_params.extraction_params_lo.num_extraction_radii;
         ++iradius)
    {
        double r = sim_params.extraction_params_lo.extraction_radii[iradius];
        double integral_re_lo_midpoint = (integral_lo_midpoint.first)[iradius];
        double integral_re_hi_midpoint = (integral_hi_midpoint.first)[iradius];
        double integral_re_lo_open_2nd_order =
            (integral_lo_open_2nd_order.first)[iradius];
        double integral_re_hi_open_2nd_order =
            (integral_hi_open_2nd_order.first)[iradius];
        double integral_re_lo_open_3rd_order =
            (integral_lo_open_3rd_order.first)[iradius];
        double integral_re_hi_open_3rd_order =
            (integral_hi_open_3rd_order.first)[iradius];
        double integral_re_lo_open_4th_order =
            (integral_lo_open_4th_order.first)[iradius];
        double integral_re_hi_open_4th_order =
            (integral_hi_open_4th_order.first)[iradius];
        double analytic_integral = 1.0;

        double convergence_factor_midpoint =
            std::abs((integral_re_lo_midpoint - analytic_integral) /
                     (integral_re_hi_midpoint - analytic_integral));
        double convergence_factor_open_2nd_order =
            std::abs((integral_re_lo_open_2nd_order - analytic_integral) /
                     (integral_re_hi_open_2nd_order - analytic_integral));
        double convergence_factor_open_3rd_order =
            std::abs((integral_re_lo_open_3rd_order - analytic_integral) /
                     (integral_re_hi_open_3rd_order - analytic_integral));
        double convergence_factor_open_4th_order =
            std::abs((integral_re_lo_open_4th_order - analytic_integral) /
                     (integral_re_hi_open_4th_order - analytic_integral));

        double convergence_order_midpoint =
            std::log2(convergence_factor_midpoint);
        double convergence_order_open_2nd_order =
            std::log2(convergence_factor_open_2nd_order);
        double convergence_order_open_3rd_order =
            std::log2(convergence_factor_open_3rd_order);
        double convergence_order_open_4th_order =
            std::log2(convergence_factor_open_4th_order);

        // midpoint rule should have second order convergence
        status |= (convergence_order_midpoint < 1.5);
        // Milne's adapted rule should have fourth order convergence
        status |= (convergence_order_open_2nd_order < 3.5);
        // Open 3rd order rule should have fourth order convergence
        status |= (convergence_order_open_3rd_order < 3.5);
        // Open 4th order rule should have sixth order convergence
        status |= (convergence_order_open_4th_order < 5.5);

        pout() << "At r = " << r << ":\n";
        pout() << "analytic_integral = " << analytic_integral << "\n";
        pout() << "integral_re_lo_midpoint = " << integral_re_lo_midpoint
               << "\n";
        pout() << "integral_re_hi_midpoint = " << integral_re_hi_midpoint
               << "\n";
        pout() << "integral_re_lo_open_2nd_order = "
               << integral_re_lo_open_2nd_order << "\n";
        pout() << "integral_re_hi_open_2nd_order = "
               << integral_re_hi_open_2nd_order << "\n";
        pout() << "integral_re_lo_open_3rd_order = "
               << integral_re_lo_open_3rd_order << "\n";
        pout() << "integral_re_hi_open_3rd_order = "
               << integral_re_hi_open_3rd_order << "\n";
        pout() << "integral_re_lo_open_4th_order = "
               << integral_re_lo_open_4th_order << "\n";
        pout() << "integral_re_hi_open_4th_order = "
               << integral_re_hi_open_4th_order << "\n";
        pout() << "convergence_order_midpoint = " << convergence_order_midpoint
               << "\n";
        pout() << "convergence_order_open_2nd_order = "
               << convergence_order_open_2nd_order << "\n";
        pout() << "convergence_order_open_3rd_order = "
               << convergence_order_open_3rd_order << "\n";
        pout() << "convergence_order_open_4th_order = "
               << convergence_order_open_4th_order << "\n"
               << endl;
    }

    return status;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runSphericalExtractionUniformTest(argc, argv);

    if (status == 0)
        pout() << "SphericalExtractionUniformTest test passed." << endl;
    else
        pout() << "SphericalExtractionUniformTest test failed with return code "
               << status << endl;

    mainFinalize();
    return status;
}

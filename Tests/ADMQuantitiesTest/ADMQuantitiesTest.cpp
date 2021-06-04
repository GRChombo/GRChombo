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

#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"

#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "ADMQuantitiesExtraction.hpp"
#include "ADMQuantitiesTestLevel.hpp"
#include "AMRInterpolator.hpp"
#include "Lagrange.hpp"

// Chombo namespace
#include "UsingNamespace.H"

int runADMQuantitiesTest(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    std::string in_string = argv[argc - 1];
    char const *in_file = argv[argc - 1];
    GRParmParse pp(0, argv + argc, NULL, in_file);
    SimulationParameters sim_params(pp);

    GRAMR gr_amr;
    DefaultLevelFactory<ADMQuantitiesTestLevel> test_level_fact(gr_amr,
                                                                sim_params);
    setupAMRObject(gr_amr, test_level_fact);

    AMRInterpolator<Lagrange<4>> interpolator(
        gr_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);

    // one integral for ADM mass and one for momentum
    std::vector<std::vector<double>> out_integrals(2);

    // perform extraction
    ADMQuantitiesExtraction adm_extraction(sim_params.extraction_params,
                                           sim_params.coarsest_dt, 0.0, true,
                                           0.0, c_Madm, c_Jadm);
    adm_extraction.extract(&interpolator);
    bool broadcast_integral = true;
    adm_extraction.add_var_integrand(
        0, out_integrals[0], IntegrationMethod::simpson,
        IntegrationMethod::trapezium, broadcast_integral);
    adm_extraction.add_var_integrand(
        1, out_integrals[1], IntegrationMethod::simpson,
        IntegrationMethod::trapezium, broadcast_integral);
    adm_extraction.integrate();

    // extrapolate results to infinity
    std::vector<double> extrapolations =
        adm_extraction.richardson_extrapolation(out_integrals);

    CH_assert(extrapolations.size() == 2);

    int status = 0;

    // spin = J / M
    double spin = extrapolations[1] / extrapolations[0];
    double mass_rel_error_percentage =
        abs(extrapolations[0] / sim_params.kerr_params.mass - 1.) * 100.;
    double spin_rel_error_percentage =
        abs(spin / sim_params.kerr_params.spin - 1.) * 100.;

    pout() << std::setprecision(10);
    pout() << "\nExtrapolation to infinity results:" << std::endl;
    pout() << "ADM Mass = " << extrapolations[0] << std::endl;
    pout() << "ADM Spin = " << spin << std::endl;
    pout() << "Expected results:" << std::endl;
    pout() << "ADM Mass = " << sim_params.kerr_params.mass << std::endl;
    pout() << "ADM Spin = " << sim_params.kerr_params.spin << std::endl;
    pout() << std::setprecision(3);
    pout() << "Relative Error:" << std::endl;
    pout() << "ADM Mass = " << mass_rel_error_percentage << "%" << std::endl;
    pout() << "ADM Spin = " << spin_rel_error_percentage << "%" << std::endl;

    status |= (mass_rel_error_percentage > 1.);
    status |= (spin_rel_error_percentage > 1.);

    return status;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runADMQuantitiesTest(argc, argv);

    if (status == 0)
        pout() << "ADMQuantitiesTest test passed." << endl;
    else
        pout() << "ADMQuantitiesTest test failed with return code " << status
               << endl;

    mainFinalize();
    return status;
}

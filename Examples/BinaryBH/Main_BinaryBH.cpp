/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CH_Timer.H"
#include "parstream.H" //Gives us pout()
#include <chrono>
#include <iostream>
#include <unistd.h>

#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "BinaryBHLevel.hpp"

int runGRChombo(int argc, char *argv[])
{
    do
    {
        // Load the parameter file and construct the SimulationParameter class
        // To add more parameters edit the SimulationParameters file.
        char *in_file = argv[1];
        GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
        SimulationParameters sim_params(pp);

        // The line below selects the problem that is simulated
        // (To simulate a different problem, define a new child of AMRLevel
        // and an associated LevelFactory)
        BHAMR bh_amr;
        DefaultLevelFactory<BinaryBHLevel> binary_bh_level_fact(bh_amr,
                                                                sim_params);
        setupAMRObject(bh_amr, binary_bh_level_fact);

        // call this after amr object setup so grids known
        // and need it to stay in scope throughout run
        AMRInterpolator<Lagrange<4>> interpolator(
            bh_amr, sim_params.origin, sim_params.dx, sim_params.verbosity);
        bh_amr.set_interpolator(&interpolator);

        using Clock = std::chrono::steady_clock;
        using Minutes = std::chrono::duration<double, std::ratio<60, 1>>;

        std::chrono::time_point<Clock> start_time = Clock::now();

        bh_amr.run(sim_params.stop_time, sim_params.max_steps);

        auto now = Clock::now();
        auto duration = std::chrono::duration_cast<Minutes>(now - start_time);
        pout() << "Total simulation time (mins): " << duration.count() << ".\n";

        bh_amr.conclude();

        CH_TIMER_REPORT(); // Report results when running with Chombo timers.
    } while (UserRestart::activate());
    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRChombo(argc, argv);

    if (status == 0)
        pout() << "GRChombo finished." << std::endl;
    else
        pout() << "GRChombo failed with return code " << status << std::endl;

    mainFinalize();
    return status;
}

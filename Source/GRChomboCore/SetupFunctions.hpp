/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETUP_FUNCTIONS_HPP_
#define SETUP_FUNCTIONS_HPP_
// This file incldues several functions that need to be called to
// set up the runs but aren't very interesting for the normal user.

#include "parstream.H" //Gives us pout()
#include <iostream>
using std::cerr;
using std::endl;
#include "AMRLevelFactory.H"
#include "ChomboParameters.hpp"
#include "DerivativeSetup.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "UserRestart.hpp"

#include "simd.hpp"

#ifdef EQUATION_DEBUG_MODE
#include "DebuggingTools.hpp"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/// This function calls MPI_Init, makes sure a parameter file is supplied etc...
void mainSetup(int argc, char *argv[]);

/// This function calls all finalisations
void mainFinalize();

/// Sets up the grid parameters, problem domain and AMR object
void setupAMRObject(AMR &gr_amr, AMRLevelFactory &a_factory);

void mainSetup(int argc, char *argv[])
{
#ifdef CH_MPI
    // Start MPI
    MPI_Init(&argc, &argv);
#ifdef CH_AIX
    H5dont_atexit();
#endif
// setChomboMPIErrorHandler();
#endif

    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank = 0;
    number_procs = 1;
#endif

#ifdef EQUATION_DEBUG_MODE
    EquationDebugging::check_no_omp();
    MayDay::Warning("GRChombo is running in equation debug mode. This mode is "
                    "intended only for debugging and leads to significantly "
                    "worse performance.");
#endif

    if (rank == 0)
    {
        pout() << " number_procs = " << number_procs << endl;
#ifdef _OPENMP
        pout() << " threads = " << omp_get_max_threads() << endl;
#endif
        pout() << " simd width (doubles) = " << simd_traits<double>::simd_len
               << endl;
    }

    const int required_argc = 2;
    if (argc < required_argc)
    {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
    }
}

void mainFinalize()
{
#ifdef CH_MPI
    // Exit MPI
#ifdef CH_USE_MEMORY_TRACKING
    dumpmemoryatexit();
#endif
    MPI_Finalize();
#endif
}

void setupAMRObject(GRAMR &gr_amr, AMRLevelFactory &a_factory)
{
    // Reread the params - just the base ones
    // Note that we could have passed these through the function
    // but this way preserves backwards compatibility
    GRParmParse pp;
    ChomboParameters chombo_params(pp);

    // remove restart trigger file
    if (chombo_params.allow_user_restart && UserRestart::activate())
    {
        UserRestart::remove_trigger(chombo_params.restart_trigger_file);
    }

    // set size of box
    Box problem_domain(IntVect::Zero, chombo_params.ivN);
    ProblemDomain physdomain(problem_domain);

    // set periodicity
    for (int dir = 0; dir < SpaceDim; dir++)
    {
        physdomain.setPeriodic(dir, chombo_params.isPeriodic[dir]);
    }

    // Define the AMR object
    gr_amr.define(chombo_params.max_level, chombo_params.ref_ratios, physdomain,
                  &a_factory);

    // Need to be able to cleanly stop evolutions if the user wants to restart
    gr_amr.allowEvolutionStop(chombo_params.allow_user_restart);

    // To preserve proper nesting we need to know the maximum ref_ratio, this
    // is now hard coded to 2 in the base params
    // The buffer is width of ghost cells + additional_grid_buffer
    // and defines the minimum number of level l cells there have to be
    // between level l+1 and level l-1
    const int max_ref_ratio = 2;
    const int additional_grid_buffer = 3;
    int grid_buffer_size =
        std::ceil(((double)chombo_params.num_ghosts) / (double)max_ref_ratio) +
        additional_grid_buffer;
    gr_amr.gridBufferSize(grid_buffer_size);

    // set checkpoint and plot intervals and prefixes
    gr_amr.checkpointInterval(chombo_params.checkpoint_interval);
    gr_amr.checkpointPrefix(chombo_params.checkpoint_prefix);
    if (chombo_params.plot_interval != 0)
    {
        gr_amr.plotInterval(chombo_params.plot_interval);
        gr_amr.plotPrefix(chombo_params.plot_prefix);
    }

    // Number of coarse time steps from one regridding to the next
    gr_amr.regridIntervals(chombo_params.regrid_interval);

    // max and min box sizes, fill ratio determining accuracy of regrid
    gr_amr.maxGridSize(chombo_params.max_grid_size);
    gr_amr.blockFactor(chombo_params.block_factor);
    gr_amr.fillRatio(chombo_params.fill_ratio);

    // Set verbosity
    gr_amr.verbosity(chombo_params.verbosity);

    // Set up input files
    if (!pp.contains("restart_file"))
    {
        gr_amr.setupForNewAMRRun();
    }
    else
    {
        std::string restart_file;
        pp.query("restart_file", restart_file);

#ifdef CH_USE_HDF5
        HDF5Handle handle(restart_file, HDF5Handle::OPEN_RDONLY);
        // read from checkpoint file
        gr_amr.setupForRestart(handle);
        handle.close();
#else
        MayDay::Error("GRChombo restart only defined with hdf5");
#endif
    }
}

#endif /* SETUP_FUNCTIONS_HPP_ */

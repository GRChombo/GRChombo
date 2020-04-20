/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHOMBOPARAMETERS_HPP_
#define CHOMBOPARAMETERS_HPP_

// General includes
#include "BoundaryConditions.hpp"
#include "GRParmParse.hpp"

class ChomboParameters
{
  public:
    ChomboParameters(GRParmParse &pp) { read_params(pp); }

    void read_params(GRParmParse &pp)
    {
        pp.load("verbosity", verbosity, 0);
        // Grid setup
        pp.load("L", L, 1.0);
        pp.load("regrid_threshold", regrid_threshold, 0.5);
        pp.load("num_ghosts", num_ghosts, 3);
        pp.load("tag_buffer_size", tag_buffer_size, 3);
        pp.load("dt_multiplier", dt_multiplier, 0.25);
        pp.load("fill_ratio", fill_ratio, 0.7);

        // Setup the grid size
        ivN = IntVect::Unit;
        int max_N = 0;
        for (int dir = 0; dir < CH_SPACEDIM; ++dir)
        {
            char dir_str[20];
            sprintf(dir_str, "N%d", dir + 1);
            int N;
            pp.load(dir_str, N);
            N_vect.push_back(N);
            ivN[dir] = N - 1;
            max_N = max(N, max_N);
        }
        coarsest_dx = L / max_N;

        pp.load("center", center,
                {0.5 * N_vect[0] * coarsest_dx, 0.5 * N_vect[1] * coarsest_dx,
                 0.5 * N_vect[2] * coarsest_dx}); // default to center

        // Periodicity and boundaries
        pp.load("isPeriodic", isPeriodic, {true, true, true});
        int bc = BoundaryConditions::STATIC_BC;
        pp.load("hi_boundary", boundary_params.hi_boundary, {bc, bc, bc});
        pp.load("lo_boundary", boundary_params.lo_boundary, {bc, bc, bc});
        // set defaults, then override them where appropriate
        boundary_params.vars_parity.fill(BoundaryConditions::EVEN);
        boundary_params.vars_asymptotic_values.fill(0.0);
        boundary_params.is_periodic.fill(true);
        nonperiodic_boundaries_exist = false;
        symmetric_boundaries_exist = false;
        FOR1(idir)
        {
            if (isPeriodic[idir] == false)
            {
                nonperiodic_boundaries_exist = true;
                boundary_params.is_periodic[idir] = false;

                // read in relevent params - note that no defaults are set so as
                // to force the user to specify them where the relevant BCs are
                // selected
                if ((boundary_params.hi_boundary[idir] ==
                     BoundaryConditions::REFLECTIVE_BC) ||
                    (boundary_params.lo_boundary[idir] ==
                     BoundaryConditions::REFLECTIVE_BC))
                {
                    symmetric_boundaries_exist = true;
                    pp.load("vars_parity", boundary_params.vars_parity);

                    if ((boundary_params.hi_boundary[idir] ==
                         BoundaryConditions::REFLECTIVE_BC) &&
                        (boundary_params.lo_boundary[idir] !=
                         BoundaryConditions::REFLECTIVE_BC))
                    {
                        center[idir] = N_vect[idir] * coarsest_dx;
                    }

                    if ((boundary_params.lo_boundary[idir] ==
                         BoundaryConditions::REFLECTIVE_BC) &&
                        (boundary_params.hi_boundary[idir] !=
                         BoundaryConditions::REFLECTIVE_BC))
                    {
                        center[idir] = 0;
                    }
                }
                if ((boundary_params.hi_boundary[idir] ==
                     BoundaryConditions::SOMMERFELD_BC) ||
                    (boundary_params.lo_boundary[idir] ==
                     BoundaryConditions::SOMMERFELD_BC))
                {
                    pp.load("vars_asymptotic_values",
                            boundary_params.vars_asymptotic_values);
                }
            }
        }
        pout() << "Center has been set to: ";
        FOR1(idir) { pout() << center[idir] << " "; }
        pout() << endl;

        if (nonperiodic_boundaries_exist)
        {
            // write out boundary conditions where non periodic - useful for
            // debug
            BoundaryConditions::write_boundary_conditions(boundary_params);
        }

        // Misc
        pp.load("ignore_checkpoint_name_mismatch",
                ignore_checkpoint_name_mismatch, false);
        pp.load("allow_user_restart", allow_user_restart, false);
        std::string default_restart_trigger_file = "GRCHOMBO_RESTART";
        pp.load("restart_trigger_file", restart_trigger_file,
                default_restart_trigger_file);

        pp.load("max_level", max_level, 0);
        // the reference ratio is hard coded to 2 on all levels
        // in principle it can be set to other values, but this is
        // not recommended since we do not test GRChombo with other
        // refinement ratios - use other values at your own risk
        ref_ratios.resize(max_level + 1);
        ref_ratios.assign(2);
        // read in frequency of regrid on each levels, needs
        // max_level + 1 entries (although never regrids on max_level+1)
        pp.getarr("regrid_interval", regrid_interval, 0, max_level + 1);

        // time stepping outputs and regrid data
        pp.load("checkpoint_interval", checkpoint_interval, 1);
        pp.load("chk_prefix", checkpoint_prefix);
        pp.load("plot_interval", plot_interval, 0);
        pp.load("plot_prefix", plot_prefix);
        pp.load("stop_time", stop_time, 1.0);
        pp.load("max_steps", max_steps, 1000000);
        pp.load("write_plot_ghosts", write_plot_ghosts, false);

        // alias the weird chombo names to something more descriptive
        // for these box params, and default to some reasonable values
        if (pp.contains("max_grid_size"))
        {
            pp.load("max_grid_size", max_grid_size);
        }
        else
        {
            pp.load("max_box_size", max_grid_size, 64);
        }
        if (pp.contains("block_factor"))
        {
            pp.load("block_factor", block_factor);
        }
        else
        {
            pp.load("min_box_size", block_factor, 8);
        }
    }

    // General parameters
    int verbosity;
    std::vector<int> N_vect;
    double L;                               // Physical sidelength of the grid
    std::array<double, CH_SPACEDIM> center; // grid center
    IntVect ivN;                 // The number of grid cells in each dimension
    double coarsest_dx;          // The coarsest resolution
    int max_level;               // the max number of regriddings to do
    int num_ghosts;              // must be at least 3 for KO dissipation
    int tag_buffer_size;         // Amount the tagged region is grown by
    Vector<int> ref_ratios;      // ref ratios between levels
    Vector<int> regrid_interval; // steps between regrid at each level
    int max_steps;
    bool ignore_checkpoint_name_mismatch;   // ignore mismatch of variable names
                                            // between restart file and program
    double dt_multiplier, stop_time;        // The Courant factor and stop time
    int checkpoint_interval, plot_interval; // Steps between outputs
    int max_grid_size, block_factor;        // max and min box sizes
    double fill_ratio; // determines how fussy the regridding is about tags
    std::string checkpoint_prefix, plot_prefix; // naming of files
    bool write_plot_ghosts;
    bool allow_user_restart; // allows the user to restart the simulation.
    std::string restart_trigger_file; // name of file to trigger user restart

    // Boundary conditions
    std::array<bool, CH_SPACEDIM> isPeriodic;     // periodicity
    BoundaryConditions::params_t boundary_params; // set boundaries in each dir
    bool nonperiodic_boundaries_exist;
    bool symmetric_boundaries_exist;

    // For tagging
    double regrid_threshold;
};

#endif /* CHOMBOPARAMETERS_HPP_ */

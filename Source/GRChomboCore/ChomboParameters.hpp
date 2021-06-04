/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHOMBOPARAMETERS_HPP_
#define CHOMBOPARAMETERS_HPP_

// Chombo includes
#include "MayDay.H"
#include "Misc.H"

// General includes
#include "ArrayTools.hpp"
#include "BoundaryConditions.hpp"
#include "FilesystemTools.hpp"
#include "GRParmParse.hpp"
#include "UserVariables.hpp"
#include "VariableType.hpp"
#include "unistd.h" // gives 'access'
#include <algorithm>
#include <string>

// Chombo namespace
#include "UsingNamespace.H"

class ChomboParameters
{
  public:
    ChomboParameters(GRParmParse &pp)
    {
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // must be before any pout() in the code to setPoutBaseName
        read_filesystem_params(pp);

        pp.load("verbosity", verbosity, 0);
        // Grid setup
        pp.load("max_spatial_derivative_order", max_spatial_derivative_order,
                4);
        pp.load("num_ghosts", num_ghosts,
                (max_spatial_derivative_order == 6) ? 4 : 3);
        pp.load("tag_buffer_size", tag_buffer_size, 3);
        pp.load("grid_buffer_size", grid_buffer_size, 8);
        pp.load("dt_multiplier", dt_multiplier, 0.25);
        pp.load("fill_ratio", fill_ratio, 0.7);

        // Periodicity and boundaries
        boundary_params.read_params(pp);

        // L's, N's and center
        read_grid_params(pp);

#ifdef CH_USE_HDF5
        pp.load("ignore_checkpoint_name_mismatch",
                ignore_checkpoint_name_mismatch, false);
#endif

        pp.load("max_level", max_level, 0);
        // the reference ratio is hard coded to 2
        // in principle it can be set to other values, but this is
        // not recommended since we do not test GRChombo with other
        // refinement ratios - use other values at your own risk
        ref_ratios.resize(max_level + 1);
        ref_ratios.assign(2);
        pp.getarr("regrid_interval", regrid_interval, 0, max_level);
        // Regridding on max_level does nothing but Chombo's AMR class
        // expects this Vector to be of length max_level + 1
        // so just set the final value to 0.
        regrid_interval.resize(max_level + 1);
        regrid_interval[max_level] = 0;

        if (pp.contains("regrid_thresholds"))
        {
            pout() << "Using multiple regrid thresholds." << std::endl;
            // As for regrid_interval, the last element is irrelevant
            pp.getarr("regrid_thresholds", regrid_thresholds, 0, max_level);
            regrid_thresholds.resize(max_level + 1);
            regrid_thresholds[max_level] = regrid_thresholds[max_level - 1];
        }
        else
        {
            pout() << "Using single regrid threshold." << std::endl;
            double regrid_threshold;
            pp.load("regrid_threshold", regrid_threshold, 0.5);
            regrid_thresholds = Vector<double>(max_level + 1, regrid_threshold);
        }

        // time stepping outputs and regrid data
        pp.load("checkpoint_interval", checkpoint_interval, 1);
        pp.load("plot_interval", plot_interval, 0);
        pp.load("stop_time", stop_time, 1000000.0);
        pp.load("max_steps", max_steps, 1000000);
#ifdef CH_USE_HDF5
        pp.load("write_plot_ghosts", write_plot_ghosts, false);

        // load vars to write to plot files
        UserVariables::load_vars_to_vector(pp, "plot_vars", "num_plot_vars",
                                           plot_vars, num_plot_vars);
#endif
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

        if (pp.contains("check_params"))
            just_check_params = true;

        pp.load("print_progress_only_to_rank_0", print_progress_only_to_rank_0,
                false);
    }

    void read_filesystem_params(GRParmParse &pp)
    {
        // In this function, cannot use default value - it may print a 'default
        // message' to pout and a 'setPoutBaseName' must happen before
        restart_from_checkpoint = pp.contains("restart_file");
#ifdef CH_USE_HDF5
        if (restart_from_checkpoint)
        {
            pp.load("restart_file", restart_file);
        }
        pp.load("chk_prefix", checkpoint_prefix);
        pp.load("plot_prefix", plot_prefix);
#endif

#ifdef CH_MPI
        // Again, cannot use default value
        if (pp.contains("pout_prefix"))
            pp.load("pout_prefix", pout_prefix);
        else
            pout_prefix = "pout";
#endif

        std::string default_path = "";
        if (pp.contains("output_path"))
            pp.load("output_path", output_path);
        else
            output_path = default_path;

#ifdef CH_MPI
        // user sets the 'subpath', we prepend 'output_path'
        if (pp.contains("pout_subpath"))
            pp.load("pout_subpath", pout_path);
        else
            pout_path = default_path;
#endif

#ifdef CH_USE_HDF5
        // user sets the 'subpath', we prepend 'output_path'
        if (pp.contains("hdf5_subpath"))
            pp.load("hdf5_subpath", hdf5_path);
        else
            hdf5_path = default_path;
#endif

        // add backslash to paths
        if (!output_path.empty() && output_path.back() != '/')
            output_path += "/";
#ifdef CH_MPI
        if (!pout_path.empty() && pout_path.back() != '/')
            pout_path += "/";
#endif
#ifdef CH_USE_HDF5
        if (!hdf5_path.empty() && hdf5_path.back() != '/')
            hdf5_path += "/";
#endif

        if (output_path != "./" && !output_path.empty())
        {
#ifdef CH_MPI
            pout_path = output_path + pout_path;
#endif
#ifdef CH_USE_HDF5
            hdf5_path = output_path + hdf5_path;
#endif
        }

#ifdef CH_MPI
        // change pout base name!
        if (!FilesystemTools::directory_exists(pout_path))
            FilesystemTools::mkdir_recursive(pout_path);
        setPoutBaseName(pout_path + pout_prefix);
#endif

        // only create hdf5 directory in setupAMRObject (when it becomes needed)
    }

    void read_grid_params(GRParmParse &pp)
    {
        // Grid N
        std::array<int, CH_SPACEDIM> Ni_full;
        std::array<int, CH_SPACEDIM> Ni;
        ivN = IntVect::Unit;

        // cannot contain both
        if ((pp.contains("N_full") && pp.contains("N")))
            MayDay::Error("Please only provide 'N' or 'N_full', not both");

        int N_full = -1;
        int N = -1;
        if (pp.contains("N_full"))
            pp.load("N_full", N_full);
        else if (pp.contains("N"))
            pp.load("N", N);

        // read all options (N, N_full, Ni_full and Ni) and then choose
        // accordingly
        FOR1(dir)
        {
            std::string name = ("N" + std::to_string(dir + 1));
            std::string name_full = ("N" + std::to_string(dir + 1) + "_full");
            Ni_full[dir] = -1;
            Ni[dir] = -1;

            // only one of them exists - this passes if none of the 4 exist, but
            // that is asserted below
            if (!((N_full > 0 || N > 0) && !pp.contains(name.c_str()) &&
                  !pp.contains(name_full.c_str())) &&
                !((N_full < 0 && N < 0) && !(pp.contains(name.c_str()) &&
                                             pp.contains(name_full.c_str()))))
                error("Please provide 'N' or 'N_full' or a set of "
                      "'N1/N1_full', 'N2/N2_full', 'N3/N3_full'");

            if (N_full < 0 && N < 0)
            {
                if (pp.contains(name_full.c_str()))
                    pp.load(name_full.c_str(), Ni_full[dir]);
                else
                    pp.load(name.c_str(), Ni[dir]);
            }
            if (N < 0 && N_full < 0 && Ni[dir] < 0 &&
                Ni_full[dir] < 0) // sanity check
                error("Please provide 'N' or 'N_full' or a set of "
                      "'N1/N1_full', 'N2/N2_full', 'N3/N3_full'");

            if (N_full > 0)
                Ni_full[dir] = N_full;
            else if (N > 0)
                Ni[dir] = N;

            if (Ni[dir] > 0)
            {
                if (boundary_params.lo_boundary[dir] ==
                        BoundaryConditions::REFLECTIVE_BC ||
                    boundary_params.hi_boundary[dir] ==
                        BoundaryConditions::REFLECTIVE_BC)

                    Ni_full[dir] = Ni[dir] * 2;
                else
                    Ni_full[dir] = Ni[dir];
            }
            else
            {
                if (boundary_params.lo_boundary[dir] ==
                        BoundaryConditions::REFLECTIVE_BC ||
                    boundary_params.hi_boundary[dir] ==
                        BoundaryConditions::REFLECTIVE_BC)
                {
                    check_parameter("N" + std::to_string(dir) + "_full",
                                    Ni_full[dir], Ni_full[dir] % 2 == 0,
                                    "must be a multiple of 2");

                    Ni[dir] = Ni_full[dir] / 2;
                }
                else
                    Ni[dir] = Ni_full[dir];
            }
            ivN[dir] = Ni[dir] - 1;
        }
        int max_N_full = *std::max_element(Ni_full.begin(), Ni_full.end());
        int max_N = ivN.max() + 1;

        // Grid L
        // cannot contain both
        if ((pp.contains("L_full") && pp.contains("L")))
            error("Please only provide 'L' or 'L_full', not both");

        double L_full = -1.;
        if (pp.contains("L_full"))
            pp.load("L_full", L_full);
        else
            pp.load("L", L, 1.0);

        if (L_full > 0.)
            // necessary for some reflective BC cases, as 'L' is the
            // length of the longest side of the box
            L = (L_full * max_N) / max_N_full;

        coarsest_dx = L / max_N;
        coarsest_dt = coarsest_dx * dt_multiplier;

        // grid spacing params
        dx.fill(coarsest_dx);
        origin.fill(coarsest_dx / 2.0);

        // These aren't parameters but used in parameter checks
        FOR1(idir)
        {
            reflective_domain_lo[idir] = ((boundary_params.lo_boundary[idir] ==
                                           BoundaryConditions::REFLECTIVE_BC)
                                              ? -1.0
                                              : 0.0) *
                                         (ivN[idir] + 1) * coarsest_dx;
            reflective_domain_hi[idir] = ((boundary_params.hi_boundary[idir] ==
                                           BoundaryConditions::REFLECTIVE_BC)
                                              ? 2.0
                                              : 1.0) *
                                         (ivN[idir] + 1) * coarsest_dx;
        }

        // First work out the default center ignoring reflective BCs
        // but taking into account different grid lengths in each direction
        std::array<double, CH_SPACEDIM> default_center;
#if CH_SPACEDIM == 3
        default_center = {0.5 * Ni[0] * coarsest_dx, 0.5 * Ni[1] * coarsest_dx,
                          0.5 * Ni[2] * coarsest_dx};
#elif CH_SPACEDIM == 2
        default_center = {0.5 * Ni[0] * coarsest_dx, 0.5 * Ni[1] * coarsest_dx};
#endif
        // Now take into account reflective BCs
        FOR1(idir)
        {
            if ((boundary_params.lo_boundary[idir] ==
                 BoundaryConditions::REFLECTIVE_BC) &&
                (boundary_params.hi_boundary[idir] !=
                 BoundaryConditions::REFLECTIVE_BC))
                default_center[idir] = 0.;
            else if ((boundary_params.hi_boundary[idir] ==
                      BoundaryConditions::REFLECTIVE_BC) &&
                     (boundary_params.lo_boundary[idir] !=
                      BoundaryConditions::REFLECTIVE_BC))
                default_center[idir] = coarsest_dx * Ni[idir];
        }

        pp.load("center", center, default_center); // default to center
    }

    void check_params()
    {
        check_parameter("L", L, L > 0.0, "must be > 0.0");
        check_parameter("max_level", max_level, max_level >= 0, "must be >= 0");
        check_parameter("max_spatial_derivative_order",
                        max_spatial_derivative_order,
                        max_spatial_derivative_order == 4 ||
                            max_spatial_derivative_order == 6,
                        "only 4 and 6 are supported");
        // the following check assumes you will be taking one-sided derivatives
        // of the order given by max_spatial_derivative_order
        check_parameter(
            "num_ghosts", num_ghosts,
            (num_ghosts >= ((max_spatial_derivative_order == 6) ? 4 : 3)) &&
                (num_ghosts <= block_factor),
            "must be >= 3 (4th order derivatives) or 4 (6th order derivatives) "
            "and <= min_box_size (aka block_factor)");
        check_parameter("tag_buffer_size", tag_buffer_size,
                        tag_buffer_size >= 0, "must be >= 0");
        // assume ref_ratio is always 2
        check_parameter(
            "grid_buffer_size", grid_buffer_size,
            grid_buffer_size >= ceil(num_ghosts / 2.0),
            "must be >= ceil(num_ghosts/max_ref_ratio) for proper nesting");

        // check the restart_file exists and can be read if restarting from a
        // checkpoint
#ifdef CH_USE_HDF5
        if (restart_from_checkpoint)
        {
            bool restart_file_exists =
                (access((hdf5_path + restart_file).c_str(), R_OK) == 0);
            check_parameter("restart_file", restart_file, restart_file_exists,
                            "file cannot be opened for reading");
        }
#endif

        check_parameter("dt_multiplier", dt_multiplier, dt_multiplier > 0.0,
                        "must be > 0.0");
        check_parameter("max_grid_size/max_box_size", max_grid_size,
                        max_grid_size >= 0, "must be >= 0");
        check_parameter("block_factor/min_box_size", block_factor,
                        block_factor >= 1, "must be >= 1");
        check_parameter("block_factor/min_box_size", block_factor,
                        Misc::isPower2(block_factor), "must be a power of 2");
        // note that this also enforces block_factor <= max_grid_size
        // if max_grid_size > 0
        check_parameter("block_factor/min_box_size", block_factor,
                        max_grid_size % block_factor == 0,
                        "must divide max_grid_size/max_box_size = " +
                            std::to_string(max_grid_size));
        FOR1(idir)
        {
            std::string Ni_string = "N" + std::to_string(idir + 1);
            std::string invalid_message = "must divide " + Ni_string;
            if (boundary_params.reflective_boundaries_exist)
            {
                invalid_message += " (or " + Ni_string + "_full/2)";
            }
            invalid_message += " = " + std::to_string(ivN[idir] + 1);
            check_parameter("block_factor/min_box_size", block_factor,
                            (ivN[idir] + 1) % block_factor == 0,
                            invalid_message);
        }
        check_parameter("fill_ratio", fill_ratio,
                        (fill_ratio > 0.0) && (fill_ratio <= 1.0),
                        "must be > 0 and <= 1");
        // (MR); while this would technically work (any plot files would just
        // overwrite a checkpoint file), I think a user would only ever do
        // this unintentinally
#ifdef CH_USE_HDF5
        check_parameter("plot_prefix", plot_prefix,
                        plot_interval <= 0 || plot_prefix != checkpoint_prefix,
                        "should be different to checkpoint_prefix");
#endif

        check_parameter("output_path", output_path,
                        FilesystemTools::directory_exists(output_path),
                        "should be a valid directory");
        // pout directory exists - we create it in read_filesystem_params()
        // can't check hdf5 directory yet - only created after

        if (boundary_params.reflective_boundaries_exist)
        {
            for (int ivar = 0; ivar < NUM_VARS; ++ivar)
            {
                std::string name = "vars_parity[c_" +
                                   UserVariables::variable_names[ivar] + "]";
                int var_parity = boundary_params.vars_parity[ivar];
                check_parameter(name, var_parity,
                                var_parity >= BoundaryConditions::EVEN &&
                                    var_parity < BoundaryConditions::UNDEFINED,
                                "parity type undefined");
            }
            for (int ivar = 0; ivar < NUM_DIAGNOSTIC_VARS; ++ivar)
            {
                std::string name = "vars_parity_diagnostic[c_" +
                                   DiagnosticVariables::variable_names[ivar] +
                                   "]";
                int var_parity = boundary_params.vars_parity_diagnostic[ivar];
                check_parameter(name, var_parity,
                                var_parity >= BoundaryConditions::EVEN &&
                                    var_parity <= BoundaryConditions::UNDEFINED,
                                "parity type undefined");
            }
        }
    }

    // General parameters
    int verbosity;
    double L;                               // Physical sidelength of the grid
    std::array<double, CH_SPACEDIM> center; // grid center
    IntVect ivN; // The number of grid cells in each dimension
    double coarsest_dx,
        coarsest_dt; // The coarsest resolution in space and time
    int max_level;   // the max number of regriddings to do
    int max_spatial_derivative_order; // The maximum order of the spatial
                                      // derivatives - does nothing
                                      // in Chombo but can be used in examples
    int num_ghosts;         // min dependent on max_spatial_derivative_order
    int tag_buffer_size;    // Amount the tagged region is grown by
    int grid_buffer_size;   // Number of cells between level
    Vector<int> ref_ratios; // ref ratios between levels
    // boundaries.
    Vector<int> regrid_interval; // steps between regrid at each level
    int max_steps;
    bool restart_from_checkpoint; // whether or not to restart or start afresh
#ifdef CH_USE_HDF5
    std::string restart_file;             // The path to the restart_file
    bool ignore_checkpoint_name_mismatch; // ignore mismatch of variable names
                                          // between restart file and program
#endif
    double dt_multiplier, stop_time;        // The Courant factor and stop time
    int checkpoint_interval, plot_interval; // Steps between outputs
    int max_grid_size, block_factor;        // max and min box sizes
    double fill_ratio; // determines how fussy the regridding is about tags
#ifdef CH_USE_HDF5
    std::string checkpoint_prefix, plot_prefix; // naming of files
#endif
    std::string output_path; // base path to use for all files
#ifdef CH_MPI
    std::string pout_prefix; // pout file prefix
    std::string pout_path;   // base path for pout files
#endif
#ifdef CH_USE_HDF5
    std::string hdf5_path; // base path for pout files
    bool write_plot_ghosts;
    int num_plot_vars;
    std::vector<std::pair<int, VariableType>>
        plot_vars; // vars to write to plot file
#endif

    std::array<double, CH_SPACEDIM> origin,
        dx; // location of coarsest origin and dx

    // Boundary conditions
    BoundaryConditions::params_t boundary_params; // set boundaries in each dir

    // For tagging
    Vector<double> regrid_thresholds;

    // For checking parameters and then exiting rather before instantiating
    // GRAMR (or child) object
    bool just_check_params = false;
    bool print_progress_only_to_rank_0;

  protected:
    // the low and high corners of the domain taking into account reflective BCs
    // only used in parameter checks hence protected
    std::array<double, CH_SPACEDIM> reflective_domain_lo, reflective_domain_hi;

    // use this error function instead of MayDay::error as this will only
    // print from rank 0
    void error(const std::string &a_error_message)
    {
        if (procID() == 0)
        {
            MayDay::Error(a_error_message.c_str());
        }
    }

    template <typename T>
    void check_parameter(const std::string &a_name, T a_value,
                         const bool a_valid,
                         const std::string &a_invalid_explanation)
    {
        if (a_valid)
            return;
        else
        {
            std::ostringstream error_message_ss;
            error_message_ss << "Parameter: " << a_name << " = " << a_value
                             << " is invalid: " << a_invalid_explanation;
            error(error_message_ss.str());
        }
    }

    template <typename T>
    void warn_parameter(const std::string &a_name, T a_value,
                        const bool a_nowarn,
                        const std::string &a_warning_explanation)
    {
        if (a_nowarn)
            return;
        else
        {
            // only print the warning from rank 0
            if (procID() == 0)
            {
                std::ostringstream warning_message_ss;
                warning_message_ss << "Parameter: " << a_name << " = "
                                   << a_value
                                   << " warning: " << a_warning_explanation;
                MayDay::Warning(warning_message_ss.str().c_str());
            }
        }
    }

    template <typename T, size_t N>
    void check_array_parameter(const std::string &a_name,
                               const std::array<T, N> &a_value,
                               const bool a_valid,
                               const std::string &a_invalid_explanation)
    {
        std::string value_str = ArrayTools::to_string(a_value);
        check_parameter(a_name, value_str, a_valid, a_invalid_explanation);
    }

    template <typename T, size_t N>
    void warn_array_parameter(const std::string &a_name,
                              const std::array<T, N> &a_value,
                              const bool a_nowarn,
                              const std::string &a_warning_explanation)
    {
        std::string value_str = ArrayTools::to_string(a_value);
        check_parameter(a_name, value_str, a_nowarn, a_warning_explanation);
    }
};

#endif /* CHOMBOPARAMETERS_HPP_ */

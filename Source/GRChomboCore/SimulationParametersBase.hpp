/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERSBASE_HPP_
#define SIMULATIONPARAMETERSBASE_HPP_

// General includes
#include "BoundaryConditions.hpp"
#include "CCZ4RHS.hpp"
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"
#include "SphericalExtraction.hpp"
#include <limits>

// add this type alias here for backwards compatibility
using extraction_params_t = SphericalExtraction::params_t;

class SimulationParametersBase : public ChomboParameters
{
  public:
    SimulationParametersBase(GRParmParse &pp) : ChomboParameters(pp)
    {
        read_params(pp);
        check_params();
    }

  private:
    void read_params(GRParmParse &pp)
    {
        // Lapse evolution
        pp.load("lapse_advec_coeff", ccz4_params.lapse_advec_coeff, 1.0);
        pp.load("lapse_coeff", ccz4_params.lapse_coeff, 2.0);
        pp.load("lapse_power", ccz4_params.lapse_power, 1.0);

        // Shift Evolution
        pp.load("shift_advec_coeff", ccz4_params.shift_advec_coeff, 0.0);
        pp.load("shift_Gamma_coeff", ccz4_params.shift_Gamma_coeff, 0.75);
        pp.load("eta", ccz4_params.eta, 1.0);

        // CCZ4 parameters
        pp.load("formulation", formulation, 0);
        pp.load("kappa1", ccz4_base_params.kappa1, 0.1);
        pp.load("kappa2", ccz4_base_params.kappa2, 0.0);
        pp.load("kappa3", ccz4_base_params.kappa3, 1.0);
        pp.load("covariantZ4", ccz4_base_params.covariantZ4, true);
        ccz4_params.kappa1 = ccz4_base_params.kappa1;
        ccz4_params.kappa2 = ccz4_base_params.kappa2;
        ccz4_params.kappa3 = ccz4_base_params.kappa3;
        ccz4_params.covariantZ4 = ccz4_base_params.covariantZ4;

        // Dissipation
        pp.load("sigma", sigma, 0.1);

        // Nan Check and min chi and lapse values
        pp.load("nan_check", nan_check, true);
        pp.load("min_chi", min_chi, 1e-4);
        pp.load("min_lapse", min_lapse, 1e-4);

        // directory to store data (extraction files, puncture data, constraint
        // norms)
        pp.load("data_subpath", data_path, std::string(""));
        if (!data_path.empty() && data_path.back() != '/')
            data_path += "/";
        if (output_path != "./" && !output_path.empty())
            data_path = output_path + data_path;

        // Extraction params
        pp.load("activate_extraction", activate_extraction, false);

        if (activate_extraction)
        {
            pp.load("num_extraction_radii",
                    extraction_params.num_extraction_radii, 1);

            // Check for multiple extraction radii, otherwise load single
            // radius/level (for backwards compatibility).
            if (pp.contains("extraction_levels"))
            {
                pp.load("extraction_levels",
                        extraction_params.extraction_levels,
                        extraction_params.num_extraction_radii);
            }
            else
            {
                pp.load("extraction_level", extraction_params.extraction_levels,
                        1, 0);
            }
            if (pp.contains("extraction_radii"))
            {
                pp.load("extraction_radii", extraction_params.extraction_radii,
                        extraction_params.num_extraction_radii);
            }
            else
            {
                pp.load("extraction_radius", extraction_params.extraction_radii,
                        1, 0.1);
            }

            pp.load("num_points_phi", extraction_params.num_points_phi, 2);
            pp.load("num_points_theta", extraction_params.num_points_theta, 5);
            if (extraction_params.num_points_theta % 2 == 0)
            {
                extraction_params.num_points_theta += 1;
                pout() << "Parameter: num_points_theta incompatible with "
                          "Simpson's "
                       << "rule so increased by 1.\n";
            }
            pp.load("extraction_center", extraction_params.center, center);

            if (pp.contains("modes"))
            {
                pp.load("num_modes", extraction_params.num_modes);
                std::vector<int> extraction_modes_vect(
                    2 * extraction_params.num_modes);
                pp.load("modes", extraction_modes_vect,
                        2 * extraction_params.num_modes);
                extraction_params.modes.resize(extraction_params.num_modes);
                for (int i = 0; i < extraction_params.num_modes; ++i)
                {
                    extraction_params.modes[i].first =
                        extraction_modes_vect[2 * i];
                    extraction_params.modes[i].second =
                        extraction_modes_vect[2 * i + 1];
                }
            }
            else
            {
                // by default extraction (l,m) = (2,0), (2,1) and (2,2)
                extraction_params.num_modes = 3;
                extraction_params.modes.resize(3);
                for (int i = 0; i < 3; ++i)
                {
                    extraction_params.modes[i].first = 2;
                    extraction_params.modes[i].second = i;
                }
            }

            pp.load("write_extraction", extraction_params.write_extraction,
                    false);

            std::string extraction_path;
            pp.load("extraction_subpath", extraction_path, data_path);
            if (!extraction_path.empty() && extraction_path.back() != '/')
                extraction_path += "/";
            if (output_path != "./" && !output_path.empty())
                extraction_path = output_path + extraction_path;

            extraction_params.data_path = data_path;
            extraction_params.extraction_path = extraction_path;

            // default names to Weyl extraction
            pp.load("extraction_file_prefix",
                    extraction_params.extraction_file_prefix,
                    std::string("Weyl4_extraction_"));
            pp.load("integral_file_prefix",
                    extraction_params.integral_file_prefix,
                    std::string("Weyl4_mode_"));

            if (pp.contains("extraction_extrapolation_order"))
            {
                int extrapolation_order;
                pp.load("extraction_extrapolation_order", extrapolation_order);
                pp.load("extraction_extrapolation_radii",
                        extraction_params.radii_idxs_for_extrapolation,
                        extrapolation_order);
                extraction_params.validate_extrapolation_radii();
            }
        }
    }

    void check_params()
    {
        check_parameter("dt_multiplier", dt_multiplier, dt_multiplier < 1.0,
                        "must be < 1.0 for stability");
        warn_parameter("dt_multiplier", dt_multiplier, dt_multiplier <= 0.5,
                       "is unlikely to be stable for > 0.5");

        check_parameter("sigma", sigma,
                        (sigma >= 0.0) && (sigma <= 2.0 / dt_multiplier),
                        "must be >= 0.0 and <= 2 / dt_multiplier for stability "
                        "(see Alcubierre p344)");
        warn_parameter("nan_check", nan_check, nan_check,
                       "should not normally be disabled");
        // not sure these are necessary hence commented out
        // check_parameter("min_chi", min_chi, (min_chi >= 0.0), "must be >=
        // 0.0"); check_parameter("min_lapse", min_lapse, (min_lapse >= 0.0)
        // "must be >= 0.0");
        check_parameter("formulation", formulation,
                        (formulation == CCZ4RHS<>::USE_CCZ4) ||
                            (formulation == CCZ4RHS<>::USE_BSSN),
                        "must be 0 or 1");
        if (formulation == CCZ4RHS<>::USE_CCZ4)
        {
            warn_parameter(
                "kappa1", ccz4_params.kappa1, ccz4_params.kappa1 > 0.0,
                "should be greater than 0.0 to damp constraints (see "
                "arXiv:1106.2254).");
            warn_parameter("kappa2", ccz4_params.kappa2,
                           ccz4_params.kappa2 > -1.0,
                           "should be greater than -1.0 to damp constraints "
                           "(see arXiv:1106.2254)");
        }
        else if (formulation == CCZ4RHS<>::USE_BSSN)
        {
            // maybe we should just set these to zero and print a warning
            // in the BSSN case
            warn_parameter("kappa1", ccz4_params.kappa1,
                           ccz4_params.kappa1 == 0.0,
                           "setting to 0.0 as required for BSSN");
            warn_parameter("kappa2", ccz4_params.kappa2,
                           ccz4_params.kappa2 == 0.0,
                           "setting to 0.0 as required for BSSN");
            warn_parameter("kappa3", ccz4_params.kappa3,
                           ccz4_params.kappa3 == 0.0,
                           "setting to 0.0 as required for BSSN");
            // no warning necessary for ccz4_params.covariantZ4
            ccz4_params.kappa1 = 0.0;
            ccz4_params.kappa2 = 0.0;
            ccz4_params.kappa3 = 0.0;
        }

        // only warn for gauge parameters as there are legitimate cases you may
        // want to deviate from the norm
        warn_parameter("lapse_advec_coeff", ccz4_params.lapse_advec_coeff,
                       min(abs(ccz4_params.lapse_advec_coeff),
                           abs(ccz4_params.lapse_advec_coeff - 1.0)) <
                           std::numeric_limits<double>::epsilon(),
                       "usually set to 0.0 or 1.0");
        warn_parameter("lapse_power", ccz4_params.lapse_power,
                       abs(ccz4_params.lapse_power - 1.0) <
                           std::numeric_limits<double>::epsilon(),
                       "set to 1.0 for 1+log slicing");
        warn_parameter("lapse_coeff", ccz4_params.lapse_coeff,
                       abs(ccz4_params.lapse_coeff - 2.0) <
                           std::numeric_limits<double>::epsilon(),
                       "set to 2.0 for 1+log slicing");
        warn_parameter("shift_Gamma_coeff", ccz4_params.shift_Gamma_coeff,
                       abs(ccz4_params.shift_Gamma_coeff - 0.75) <
                           std::numeric_limits<double>::epsilon(),
                       "usually set to 0.75");
        warn_parameter("shift_advec_coeff", ccz4_params.shift_advec_coeff,
                       min(abs(ccz4_params.shift_advec_coeff),
                           abs(ccz4_params.shift_advec_coeff - 1.0)) <
                           std::numeric_limits<double>::epsilon(),
                       "usually set to 0.0 or 1.0");
        warn_parameter("eta", ccz4_params.eta,
                       ccz4_params.eta > 0.1 && ccz4_params.eta < 10,
                       "usually O(1/M_ADM) so typically O(1) in code units");

        // Now extraction parameters
        if (activate_extraction)
        {
            check_parameter(
                "num_extraction_radii", extraction_params.num_extraction_radii,
                extraction_params.num_extraction_radii > 0,
                "must be bigger than 0 when activate_extraction = 1");

            FOR1(idir)
            {
                std::string center_name =
                    "extraction_center[" + std::to_string(idir) + "]";
                double center_in_dir = extraction_params.center[idir];
                check_parameter(
                    center_name, center_in_dir,
                    (center_in_dir >= reflective_domain_lo[idir]) &&
                        (center_in_dir <= reflective_domain_hi[idir]),
                    "must be in the computational domain after "
                    "applying reflective symmetry");
                for (int iradius = 0;
                     iradius < extraction_params.num_extraction_radii;
                     ++iradius)
                {
                    std::string radius_name =
                        "extraction_radii[" + std::to_string(iradius) + "]";
                    double radius = extraction_params.extraction_radii[iradius];
                    if (idir == 0)
                        check_parameter(radius_name, radius, radius >= 0.0,
                                        "must be >= 0.0");
                    check_parameter(
                        radius_name, radius,
                        (center_in_dir - radius >=
                         reflective_domain_lo[idir]) &&
                            (center_in_dir + radius <=
                             reflective_domain_hi[idir]),
                        "extraction sphere must lie within the computational "
                        "domain after applying reflective symmetry");
                }
            }
            for (int imode = 0; imode < extraction_params.num_modes; ++imode)
            {
                auto &mode = extraction_params.modes[imode];
                int l = mode.first;
                int m = mode.second;
                std::string mode_name = "modes[" + std::to_string(imode) + "]";
                std::string value_str = "(" + std::to_string(mode.first) +
                                        ", " + std::to_string(mode.second) +
                                        ")";
                check_parameter(
                    mode_name, value_str, (l >= 2) && (abs(m) <= l),
                    "l must be >= 2 and m must satisfy -l <= m <= l");
            }
            int extrapolation_order =
                extraction_params.radii_idxs_for_extrapolation.size();
            if (extrapolation_order > 1)
                pout() << "Extraction with extrapolation order "
                       << extrapolation_order << std::endl;
            else if (extrapolation_order ==
                     1) // 1 radius only because other were invalid
                pout() << "Extraction radii not valid." << std::endl;
        }
    }

  protected:
    // This is just the CCZ4 damping parameters in case you want to use
    // a different gauge (with different parameters)
    CCZ4_base_params_t ccz4_base_params;

  public:
    double sigma; // Kreiss-Oliger dissipation parameter

    bool nan_check;

    double min_chi, min_lapse;

    int formulation; // Whether to use BSSN or CCZ4

    // Collection of parameters necessary for the CCZ4 RHS
    // Note the gauge parameters are specific to MovingPunctureGauge
    // If you are using a different gauge, you need to load your parameters
    // in your own SimulationParameters class.
    CCZ4_params_t<> ccz4_params;

    bool activate_extraction;
    SphericalExtraction::params_t extraction_params;

    std::string data_path;
};

#endif /* SIMULATIONPARAMETERSBASE_HPP_ */

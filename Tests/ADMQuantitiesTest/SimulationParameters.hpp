/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"

// Problem specific includes
#include "KerrBH.hpp"
#include "SphericalExtraction.hpp"

class SimulationParameters : public ChomboParameters
{
  public:
    SimulationParameters(GRParmParse &pp) : ChomboParameters(pp)
    {
        read_params(pp);
    }

  private:
    void read_params(GRParmParse &pp)
    {
        // Initial Kerr data
        pp.load("kerr_mass", kerr_params.mass);
        pp.load("kerr_spin", kerr_params.spin);
        pp.load("kerr_center", kerr_params.center, center);

        // Extraction params
        pp.load("num_extraction_radii", extraction_params.num_extraction_radii,
                1);

        if (pp.contains("extraction_radii"))
        {
            pp.load("extraction_radii", extraction_params.extraction_radii,
                    extraction_params.num_extraction_radii);
        }
        else
        {
            pp.load("extraction_radius", extraction_params.extraction_radii, 1,
                    0.1);
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
        pp.load("write_extraction", extraction_params.write_extraction, false);

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

  public:
    KerrBH::params_t kerr_params;
    SphericalExtraction::params_t extraction_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CYLINDRICALEXTRACTION_HPP_
#define CYLINDRICALEXTRACTION_HPP_

#include "CylindricalGeometry.hpp"
#include "SurfaceExtraction.hpp"

//! A child class of SurfaceExtraction for extraction on cylindrical shells
class CylindricalExtraction : public SurfaceExtraction<CylindricalGeometry>
{
  public:
    struct params_t : surface_extraction_params_t
    {
        int &num_extraction_radii = num_surfaces;
        std::vector<double> &extraction_radii = surface_param_values;
        int &num_points_z = num_points_u;
        int &num_points_phi = num_points_v;
        std::array<double, CH_SPACEDIM> center; //!< the center of the
                                                //!< cylindrical shells
        double z_length;
    };
    const std::array<double, CH_SPACEDIM> m_center;
    double m_z_length;

    CylindricalExtraction(const params_t &a_params, double a_dt, double a_time,
                          bool a_first_step, double a_restart_time = 0.0)
        : SurfaceExtraction({a_params.center, a_params.z_length}, a_params,
                            a_dt, a_time, a_first_step, a_restart_time),
          m_center(a_params.center), m_z_length(a_params.z_length)
    {
    }

    CylindricalExtraction(const params_t &a_params,
                          const std::vector<vars_t> &a_vars, double a_dt,
                          double a_time, bool a_first_step,
                          double a_restart_time = 0.0)
        : CylindricalExtraction(a_params, a_dt, a_time, a_first_step,
                                a_restart_time)
    {
        add_vars(a_vars);
    }

    CylindricalExtraction(const params_t &a_params,
                          const std::vector<int> &a_evolution_vars, double a_dt,
                          double a_time, bool a_first_step,
                          double a_restart_time = 0.0)
        : CylindricalExtraction(a_params, a_dt, a_time, a_first_step,
                                a_restart_time)
    {
        add_evolution_vars(a_evolution_vars);
    }

    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the fields on the cylinders
        extract(a_interpolator);

        if (m_params.write_extraction)
        {
            write_extraction("CylinderExtractionOut_");
        }
    }
};

#endif /* CYLINDRICALEXTRACTION_HPP_ */

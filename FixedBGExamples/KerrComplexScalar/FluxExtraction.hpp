/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FLUXEXTRACTION_HPP_
#define FLUXEXTRACTION_HPP_

#include "SphericalExtraction.hpp"
//!  The class allows extraction of the values of the Flux components on
//!  spherical shells at specified radii, and integration over those shells
/*!
   The class allows the user to extract data from the grid for the Flux
   components over spherical shells at specified radii. The values may then be
   written to an output file, or integrated across the surfaces.
*/
class FluxExtraction : public SphericalExtraction
{
  public:
    //! The constructor
    FluxExtraction(SphericalExtraction::params_t &a_params, double a_dt,
                    double a_time, bool a_first_step,
                    double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                               a_restart_time)
    {
        add_var(c_Edot, VariableType::diagnostic);
        add_var(c_Jdot, VariableType::diagnostic);
    }

    //! The old constructor which assumes it is called in specificPostTimeStep
    //! so the first time step is when m_time == m_dt
    FluxExtraction(SphericalExtraction::params_t a_params, double a_dt,
                    double a_time, double a_restart_time = 0.0)
        : FluxExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                          a_restart_time)
    {
    }

    // the references of the vars as used in the integrator
    enum M_VARS
    {
        m_Edot,
        m_Jdot
    };

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the Flux scalars on the spheres
        extract(a_interpolator);

        // this would write out the values at every point on the sphere
        if (m_params.write_extraction)
        {
            write_extraction("Flux4ExtractionOut_");
        }

        // Setup to integrate Edot and Jdot
        std::vector<std::vector<double>> Flux_integrals(2);
        add_var_integrand(m_Edot, Flux_integrals[m_Edot],
                          IntegrationMethod::simpson);
        add_var_integrand(m_Jdot, Flux_integrals[m_Jdot],
                          IntegrationMethod::simpson);

        // do the integration over the surface
        integrate();

        // write the integrals
        std::vector<std::string> labels(2);
        labels[m_Edot] = "Edot";
        labels[m_Jdot] = "Jdot";
        std::string filename = "Flux_integrals";
        write_integrals(filename, Flux_integrals, labels);
    }
};

#endif /* FLUXEXTRACTION_HPP_ */

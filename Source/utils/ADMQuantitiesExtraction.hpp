/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMQUANTITIESEXTRACTION_HPP_
#define ADMQUANTITIESEXTRACTION_HPP_

#include "SphericalExtraction.hpp"
//!  The class allows extraction of the values of the ADM mass and angular
//!  momentum over spherical shells at specified radii, and integration over
//!  those shells
/*!
   The class allows the user to extract data from the grid for the ADM mass and
   angular momentum over spherical shells at specified radii. The values may
   then be written to an output file, or integrated across the surfaces.
*/
class ADMQuantitiesExtraction : public SphericalExtraction
{
  public:
    //! The constructor
    ADMQuantitiesExtraction(SphericalExtraction::params_t &a_params,
                            double a_dt, double a_time, bool a_first_step,
                            double a_restart_time = 0.0, int a_c_Madm = -1,
                            int a_c_Jadm = -1)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time),
          m_c_Madm(a_c_Madm), m_c_Jadm(a_c_Jadm)
    {
        if (m_c_Madm >= 0)
            add_var(m_c_Madm, VariableType::diagnostic);
        if (m_c_Jadm >= 0)
            add_var(m_c_Jadm, VariableType::diagnostic);
    }

    //! The old constructor which assumes it is called in specificPostTimeStep
    //! so the first time step is when m_time == m_dt
    ADMQuantitiesExtraction(SphericalExtraction::params_t a_params, double a_dt,
                            double a_time, double a_restart_time = 0.0,
                            int a_c_Madm = -1, int a_c_Jadm = -1)
        : ADMQuantitiesExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                                  a_restart_time, a_c_Madm, a_c_Jadm)
    {
    }

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the ADM mass and spin on the spheres
        extract(a_interpolator);

        if (m_params.write_extraction)
            write_extraction("ADMQuantitiesExtractionOut_");

        int num_integrals = (m_c_Madm >= 0) + (m_c_Jadm >= 0);
        if (!num_integrals)
            return;

        int J_index = m_c_Madm >= 0 ? 1 : 0;

        std::vector<std::vector<double>> out_integrals(num_integrals);

        if (m_c_Madm >= 0)
            add_var_integrand(0, out_integrals[0], IntegrationMethod::simpson);
        if (m_c_Jadm >= 0)
            add_var_integrand(J_index, out_integrals[J_index],
                              IntegrationMethod::simpson);

        // do the integration over the surface
        integrate();
        std::vector<std::string> labels(num_integrals);
        if (m_c_Madm >= 0)
            labels[0] = "M_adm";
        if (m_c_Jadm >= 0)
            labels[J_index] = "J_adm";
        write_integrals("IntegralADMQuantities", out_integrals, labels);
    }

  private:
    const int m_c_Madm, m_c_Jadm;
};

#endif /* ADMQUANTITIESEXTRACTION_HPP_ */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WEYLEXTRACTION_HPP_
#define WEYLEXTRACTION_HPP_

#include "SphericalExtraction.hpp"

//!  The class allows extraction of the values of the Weyl scalar components on
//!  spherical shells at specified radii, and integration over those shells
/*!
   The class allows the user to extract data from the grid for the Weyl
   components over spherical shells at specified radii. The values may then be
   written to an output file, or integrated across the surfaces.
*/
class WeylExtraction : public SphericalExtraction
{
  public:
    //! The constructor
    WeylExtraction(SphericalExtraction::params_t &a_params, double a_dt,
                   double a_time, bool a_first_step,
                   double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        add_var(c_Weyl4_Re, VariableType::diagnostic);
        add_var(c_Weyl4_Im, VariableType::diagnostic);
    }

    //! The old constructor which assumes it is called in specificPostTimeStep
    //! so the first time step is when m_time == m_dt
    WeylExtraction(SphericalExtraction::params_t a_params, double a_dt,
                   double a_time, double a_restart_time = 0.0)
        : WeylExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                         a_restart_time)
    {
    }

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the Weyl scalars on the spheres
        extract(a_interpolator);

        if (m_params.write_extraction)
            write_extraction(m_params.extraction_file_prefix);

        // now calculate and write the requested spherical harmonic modes
        std::vector<std::pair<std::vector<double>, std::vector<double>>>
            mode_integrals(m_num_modes);

        // note that this is normalised by multiplying by radius
        auto normalised_Weyl4_complex = [](std::vector<double> Weyl4_reim_parts,
                                           double r, double, double) {
            // here the std::vector<double> passed will just have
            // the real and imaginary parts of the Weyl4 scalar as its
            // only components
            return std::make_pair(r * Weyl4_reim_parts[0],
                                  r * Weyl4_reim_parts[1]);
        };

        // add the modes that will be integrated
        for (int imode = 0; imode < m_num_modes; ++imode)
        {
            const auto &mode = m_modes[imode];
            constexpr int es = -2;
            add_mode_integrand(es, mode.first, mode.second,
                               normalised_Weyl4_complex, mode_integrals[imode]);
        }

        // do the integration over the surface
        integrate();

        // write the integrals
        for (int imode = 0; imode < m_num_modes; ++imode)
        {
            const auto &mode = m_modes[imode];
            std::string integrals_filename = m_params.integral_file_prefix +
                                             std::to_string(mode.first) +
                                             std::to_string(mode.second);
            std::vector<std::vector<double>> integrals_for_writing = {
                std::move(mode_integrals[imode].first),
                std::move(mode_integrals[imode].second)};
            std::vector<std::string> labels = {"integral Re", "integral Im"};
            write_integrals(integrals_filename, integrals_for_writing, labels);
            write_integral(integrals_filename, integrals_for_writing[0],
                           labels[0]);
        }
    }
};

#endif /* WEYLEXTRACTION_HPP_ */

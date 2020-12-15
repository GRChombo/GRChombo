/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXSCALARPOTENTIAL_HPP_
#define COMPLEXSCALARPOTENTIAL_HPP_

#include "simd.hpp"

class ComplexScalarPotential
{
  public:
    struct params_t
    {
        double scalar_mass;
    };

  private:
    const params_t m_params;

  public:
    //! The constructor
    ComplexScalarPotential(const params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi_Re,
                           data_t &dVdphi_Im, const vars_t<data_t> &vars) const
    {
        const double m = m_params.scalar_mass;

        // The potential value at phi
        V_of_phi = 0.5 * m * m * vars.phi_Re * vars.phi_Re +
                   0.5 * m * m * vars.phi_Im * vars.phi_Im;

        // The potential gradient at phi wrt real and im fields
        dVdphi_Re = m * m * vars.phi_Re;
        dVdphi_Im = m * m * vars.phi_Im;
    }
};

#endif /* COMPLEXSCALARPOTENTIAL_HPP_ */

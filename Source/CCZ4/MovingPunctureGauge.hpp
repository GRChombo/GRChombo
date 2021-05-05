/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MOVINGPUNCTUREGAUGE_HPP_
#define MOVINGPUNCTUREGAUGE_HPP_


#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

/// This is an example of a gauge class that can be used in the CCZ4RHS compute
/// class
/**
 * This class implements a slightly more generic version of the moving puncture
 * gauge. In particular it uses a Bona-Masso slicing condition of the form
 * f(lapse) = -c*lapse^(p-2)
 * and a Gamma-driver shift condition
 **/
class MovingPunctureGauge
{
  public:
    struct params_t
    {
        // lapse params:
        double lapse_advec_coeff = 0.; //!< Switches advection terms in
                                       //! the lapse condition on/off
        double lapse_power = 1.; //!< The power p in \f$\partial_t \alpha = - c
                                 //!\alpha^p(K-2\Theta)\f$
        double lapse_coeff = 2.; //!< The coefficient c in \f$\partial_t \alpha
                                 //!= -c \alpha^p(K-2\Theta)\f$
        // shift params:
        double shift_Gamma_coeff = 0.75; //!< Gives the F in \f$\partial_t
                                         //!  \beta^i =  F B^i\f$
        double shift_advec_coeff = 0.;   //!< Switches advection terms in the
                                         //! shift condition on/off
        double eta = 1.; //!< The eta in \f$\partial_t B^i = \partial_t \tilde
                         //!\Gamma - \eta B^i\f$
    };

  protected:
    params_t m_params;

  public:
    MovingPunctureGauge(const params_t &a_params) : m_params(a_params) {}

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    inline void rhs_gauge(vars_t<data_t> &rhs, const vars_t<data_t> &vars,
                          const vars_t<Tensor<1, data_t>> &d1,
                          const diff2_vars_t<Tensor<2, data_t>> &d2,
                          const vars_t<data_t> &advec,
                          const Cell<data_t> &current_cell) const
    {
        rhs.lapse = m_params.lapse_advec_coeff * advec.lapse -
                    m_params.lapse_coeff *
                        pow(vars.lapse, m_params.lapse_power) *
                        (vars.K - 2 * vars.Theta);
        FOR1(i)
        {
            rhs.shift[i] = m_params.shift_advec_coeff * advec.shift[i] +
                           m_params.shift_Gamma_coeff * vars.B[i];
            rhs.B[i] = m_params.shift_advec_coeff * advec.B[i] -
                       m_params.shift_advec_coeff * advec.Gamma[i] +
                       rhs.Gamma[i] - m_params.eta * vars.B[i];
        }
    }
};

#endif /* MOVINGPUNCTUREGAUGE_HPP_ */

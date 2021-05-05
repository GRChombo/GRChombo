/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTEGRATEDMOVINGPUNCTUREGAUGE_HPP_
#define INTEGRATEDMOVINGPUNCTUREGAUGE_HPP_

#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "MovingPunctureGauge.hpp"
#include "Tensor.hpp"

/// This is an example of a gauge class that can be used in the CCZ4RHS compute
/// class
/**
 * This class implements a slightly more generic version of the moving puncture
 * gauge. In particular it uses a Bona-Masso slicing condition of the form
 * f(lapse) = -c*lapse^(p-2)
 * and an Integrated version of the Gamma-driver shift condition
 * (see details in arXiv:gr-qc/0605030)
 **/
class IntegratedMovingPunctureGauge
{
  public:
    using params_t = MovingPunctureGauge::params_t;

  protected:
    params_t m_params;

    /// Vars needed internally in 'compute'
    template <class data_t> struct Vars
    {
        Tensor<1, data_t> shift;
        Tensor<1, data_t> Gamma; //!< Conformal connection functions

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(
                mapping_function, GRInterval<c_shift1, c_shift3>(), shift);
            VarsTools::define_enum_mapping(
                mapping_function, GRInterval<c_Gamma1, c_Gamma3>(), Gamma);
        }
    };

  public:
    IntegratedMovingPunctureGauge(const params_t &a_params) : m_params(a_params)
    {
    }

    // set the initial B^i to the initial condition equivalent to:
    // \partial_t shift - advec_coeff * advec.shift = 0
    // Include in your Example in GRAMRLevel::initial_data as:
    // fillAllGhosts();
    // BoxLoops::loop(IntegratedMovingPunctureGauge(m_p.ccz4_params),
    // m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();

        Tensor<1, data_t> B;
        FOR1(i)
        {
            B[i] = m_params.shift_Gamma_coeff * vars.Gamma[i] -
                   m_params.eta * vars.shift[i];
        }

        current_cell.store_vars(B, GRInterval<c_B1, c_B3>());
    }

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
                           m_params.shift_Gamma_coeff * vars.Gamma[i] -
                           m_params.eta * vars.shift[i] - vars.B[i];
            rhs.B[i] = 0.; // static, stays the same to save initial condition
        }
    }
};

#endif /* INTEGRATEDMOVINGPUNCTUREGAUGE_HPP_ */

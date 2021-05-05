/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CCZ4RHS_HPP_)
#error "This file should only be included through CCZ4RHS.hpp"
#endif

#ifndef CCZ4RHS_IMPL_HPP_
#define CCZ4RHS_IMPL_HPP_

#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"

template <class gauge_t, class deriv_t>
inline CCZ4RHS<gauge_t, deriv_t>::CCZ4RHS(
    CCZ4_params_t<typename gauge_t::params_t> a_params,
    const gauge_t &a_gauge,
    double a_dx, double a_sigma, int a_formulation,
    double a_cosmological_constant)
    : m_params(a_params), m_gauge(a_gauge), m_sigma(a_sigma),
      m_formulation(a_formulation),
      m_cosmological_constant(a_cosmological_constant), m_deriv(a_dx)
{
    // A user who wants to use BSSN should also have damping paramters = 0
    if (m_formulation == USE_BSSN)
    {
        if ((m_params.kappa1 != 0.) || (m_params.kappa2 != 0.) ||
            (m_params.kappa3 != 0.))
        {
            MayDay::Error("BSSN formulation is selected - CCZ4 kappa values "
                          "should be set to zero in params");
        }
    }
    if (m_formulation > USE_BSSN)
        MayDay::Error("The requested formulation is not supported");
}

template <class gauge_t, class deriv_t>
template <typename dummy_t>
inline CCZ4RHS<gauge_t, deriv_t>::CCZ4RHS(
    CCZ4_params_t<typename gauge_t::params_t> a_params, double a_dx,
    double a_sigma, int a_formulation, double a_cosmological_constant,
    typename std::enable_if_t<std::is_constructible<dummy_t, params_t>::value>
        *)
    : CCZ4RHS(a_params, gauge_t(a_params), a_dx, a_sigma, a_formulation,
              a_cosmological_constant)
{
}

template <class gauge_t, class deriv_t>
template <class data_t>
void CCZ4RHS<gauge_t, deriv_t>::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
        m_deriv.template advection<Vars>(current_cell, vars.shift);

    Vars<data_t> rhs;
    rhs_equation(rhs, vars, d1, d2, advec, current_cell);

    m_deriv.add_dissipation(rhs, current_cell, m_sigma);

    current_cell.store_vars(rhs); // Write the rhs into the output FArrayBox
}

template <class gauge_t, class deriv_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
void CCZ4RHS<gauge_t, deriv_t>::rhs_equation(
    vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec,
    const Cell<data_t> &current_cell) const
{
    using namespace TensorAlgebra;

    auto h_UU = compute_inverse_sym(vars.h);
    auto chris = compute_christoffel(d1.h, h_UU);

    Tensor<1, data_t> Z_over_chi;
    Tensor<1, data_t> Z;

    if (m_formulation == USE_BSSN)
    {
        FOR1(i) Z_over_chi[i] = 0.0;
    }
    else
    {
        FOR1(i) Z_over_chi[i] = 0.5 * (vars.Gamma[i] - chris.contracted[i]);
    }
    FOR1(i) Z[i] = vars.chi * Z_over_chi[i];

    auto ricci =
        CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, Z_over_chi);

    data_t divshift = compute_trace(d1.shift);
    data_t Z_dot_d1lapse = compute_dot_product(Z, d1.lapse);
    data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);

    Tensor<2, data_t> covdtilde2lapse;
    Tensor<2, data_t> covd2lapse;
    FOR2(k, l)
    {
        covdtilde2lapse[k][l] = d2.lapse[k][l];
        FOR1(m) { covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m]; }
        covd2lapse[k][l] =
            vars.chi * covdtilde2lapse[k][l] +
            0.5 * (d1.lapse[k] * d1.chi[l] + d1.chi[k] * d1.lapse[l] -
                   vars.h[k][l] * dlapse_dot_dchi);
    }

    data_t tr_covd2lapse = -(GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR1(i)
    {
        tr_covd2lapse -= vars.chi * chris.contracted[i] * d1.lapse[i];
        FOR1(j)
        {
            tr_covd2lapse += h_UU[i][j] * (vars.chi * d2.lapse[i][j] +
                                           d1.lapse[i] * d1.chi[j]);
        }
    }

    Tensor<2, data_t> A_UU = raise_all(vars.A, h_UU);

    // A^{ij} A_{ij}. - Note the abuse of the compute trace function.
    data_t tr_A2 = compute_trace(vars.A, A_UU);
    rhs.chi = advec.chi +
              (2.0 / GR_SPACEDIM) * vars.chi * (vars.lapse * vars.K - divshift);
    FOR2(i, j)
    {
        rhs.h[i][j] = advec.h[i][j] - 2.0 * vars.lapse * vars.A[i][j] -
                      (2.0 / GR_SPACEDIM) * vars.h[i][j] * divshift;
        FOR1(k)
        {
            rhs.h[i][j] +=
                vars.h[k][i] * d1.shift[k][j] + vars.h[k][j] * d1.shift[k][i];
        }
    }

    Tensor<2, data_t> Adot_TF;
    FOR2(i, j)
    {
        Adot_TF[i][j] =
            -covd2lapse[i][j] + vars.chi * vars.lapse * ricci.LL[i][j];
    }
    make_trace_free(Adot_TF, vars.h, h_UU);

    FOR2(i, j)
    {
        rhs.A[i][j] = advec.A[i][j] + Adot_TF[i][j] +
                      vars.A[i][j] * (vars.lapse * (vars.K - 2 * vars.Theta) -
                                      (2.0 / GR_SPACEDIM) * divshift);
        FOR1(k)
        {
            rhs.A[i][j] +=
                vars.A[k][i] * d1.shift[k][j] + vars.A[k][j] * d1.shift[k][i];
            FOR1(l)
            {
                rhs.A[i][j] -=
                    2 * vars.lapse * h_UU[k][l] * vars.A[i][k] * vars.A[l][j];
            }
        }
    }

    data_t kappa1_times_lapse;
    if (m_params.covariantZ4)
        kappa1_times_lapse = m_params.kappa1;
    else
        kappa1_times_lapse = m_params.kappa1 * vars.lapse;

    if (m_formulation == USE_BSSN)
    {
        rhs.Theta = 0; // ensure the Theta of CCZ4 remains at zero
        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        rhs.K = advec.K + vars.lapse * (tr_A2 + vars.K * vars.K / GR_SPACEDIM) -
                tr_covd2lapse;
        rhs.K += -2 * vars.lapse * m_cosmological_constant / (GR_SPACEDIM - 1.);
    }
    else
    {
        rhs.Theta =
            advec.Theta +
            0.5 * vars.lapse *
                (ricci.scalar - tr_A2 +
                 ((GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) * vars.K * vars.K -
                 2 * vars.Theta * vars.K) -
            0.5 * vars.Theta * kappa1_times_lapse *
                ((GR_SPACEDIM + 1) + m_params.kappa2 * (GR_SPACEDIM - 1)) -
            Z_dot_d1lapse;

        rhs.Theta += -vars.lapse * m_cosmological_constant;
        rhs.K =
            advec.K +
            vars.lapse * (ricci.scalar + vars.K * (vars.K - 2 * vars.Theta)) -
            kappa1_times_lapse * GR_SPACEDIM * (1 + m_params.kappa2) *
                vars.Theta -
            tr_covd2lapse;
        rhs.K += -2 * vars.lapse * GR_SPACEDIM / (GR_SPACEDIM - 1.) *
                 m_cosmological_constant;
    }

    Tensor<1, data_t> Gammadot;
    FOR1(i)
    {
        Gammadot[i] = (2.0 / GR_SPACEDIM) *
                          (divshift * (chris.contracted[i] +
                                       2 * m_params.kappa3 * Z_over_chi[i]) -
                           2 * vars.lapse * vars.K * Z_over_chi[i]) -
                      2 * kappa1_times_lapse * Z_over_chi[i];
        FOR1(j)
        {
            Gammadot[i] +=
                2 * h_UU[i][j] *
                    (vars.lapse * d1.Theta[j] - vars.Theta * d1.lapse[j]) -
                2 * A_UU[i][j] * d1.lapse[j] -
                vars.lapse * ((2 * (GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                                  h_UU[i][j] * d1.K[j] +
                              GR_SPACEDIM * A_UU[i][j] * d1.chi[j] / vars.chi) -
                (chris.contracted[j] + 2 * m_params.kappa3 * Z_over_chi[j]) *
                    d1.shift[i][j];

            FOR1(k)
            {
                Gammadot[i] +=
                    2 * vars.lapse * chris.ULL[i][j][k] * A_UU[j][k] +
                    h_UU[j][k] * d2.shift[i][j][k] +
                    ((GR_SPACEDIM - 2.0) / (double)GR_SPACEDIM) * h_UU[i][j] *
                        d2.shift[k][j][k];
            }
        }
    }

    FOR1(i) { rhs.Gamma[i] = advec.Gamma[i] + Gammadot[i]; }

    m_gauge.rhs_gauge(rhs, vars, d1, d2, advec, current_cell);
}

#endif /* CCZ4RHS_IMPL_HPP_ */

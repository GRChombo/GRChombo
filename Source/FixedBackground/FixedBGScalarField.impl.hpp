/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(FIXEDBGSCALARFIELD_HPP_)
#error "This file should only be included through FixedBGScalarField.hpp"
#endif

#ifndef FIXEDBGSCALARFIELD_IMPL_HPP_
#define FIXEDBGSCALARFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> FixedBGScalarField<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars, const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &gamma_UU,
    const Tensor<3, data_t> &chris_phys_ULL) const
{
    emtensor_t<data_t> out;

    // call the function which computes the em tensor excluding the potential
    emtensor_excl_potential(out, vars, metric_vars, d1, gamma_UU,
                            chris_phys_ULL);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi = 0.0;
    my_potential.compute_potential(V_of_phi, dVdphi, vars);

    out.rho += V_of_phi;
    out.S += -3.0 * V_of_phi;
    FOR2(i, j) { out.Sij[i][j] += -metric_vars.gamma[i][j] * V_of_phi; }

    return out;
}

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void FixedBGScalarField<potential_t>::emtensor_excl_potential(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &gamma_UU, const Tensor<3, data_t> &chris_phys_ULL)
{
    // Useful quantity Vt
    data_t Vt = -vars.Pi * vars.Pi;
    FOR2(i, j) { Vt += gamma_UU[i][j] * d1.phi[i] * d1.phi[j]; }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] =
            -0.5 * metric_vars.gamma[i][j] * Vt + d1.phi[i] * d1.phi[j];
    }

    // S = Tr_S_ij
    out.S = TensorAlgebra::compute_trace(out.Sij, gamma_UU);

    // S_i (note lower index) = - n^a T_a0
    FOR1(i) { out.Si[i] = -d1.phi[i] * vars.Pi; }

    // rho = n^a n^b T_ab
    out.rho = vars.Pi * vars.Pi + 0.5 * Vt;
}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FixedBGScalarField<potential_t>::matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // call the function for the rhs excluding the potential
    matter_rhs_excl_potential(total_rhs, vars, metric_vars, d1, d2, advec);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi = 0.0;

    // compute potential
    my_potential.compute_potential(V_of_phi, dVdphi, vars);

    // adjust RHS for the potential term
    total_rhs.Pi += -metric_vars.lapse * dVdphi;
}

// the RHS excluding the potential terms
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FixedBGScalarField<potential_t>::matter_rhs_excl_potential(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec)
{
    using namespace TensorAlgebra;

    const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
    const auto chris_phys = compute_christoffel(metric_vars.d1_gamma, gamma_UU);

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs.phi = metric_vars.lapse * vars.Pi + advec.phi;
    rhs.Pi = metric_vars.lapse * metric_vars.K * vars.Pi + advec.Pi;

    FOR2(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs.Pi += gamma_UU[i][j] * (metric_vars.lapse * d2.phi[i][j] +
                                    metric_vars.d1_lapse[i] * d1.phi[j]);
        FOR1(k)
        {
            rhs.Pi += -metric_vars.lapse * gamma_UU[i][j] *
                      chris_phys.ULL[k][i][j] * d1.phi[k];
        }
    }
}

#endif /* FIXEDBGSCALARFIELD_IMPL_HPP_ */

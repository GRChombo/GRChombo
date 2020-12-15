/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(FIXEDBGCOMPLEXSCALARFIELD_HPP_)
#error "This file should only be included through FixedBGComplexScalarField.hpp"
#endif

#ifndef FIXEDBGCOMPLEXSCALARFIELD_IMPL_HPP_
#define FIXEDBGCOMPLEXSCALARFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> FixedBGComplexScalarField<potential_t>::compute_emtensor(
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
    data_t dVdphi_Re = 0.0;
    data_t dVdphi_Im = 0.0;
    my_potential.compute_potential(V_of_phi, dVdphi_Re, dVdphi_Im, vars);

    out.rho += V_of_phi;
    out.S += -3.0 * V_of_phi;
    FOR2(i, j) { out.Sij[i][j] += -metric_vars.gamma[i][j] * V_of_phi; }

    return out;
}

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void FixedBGComplexScalarField<potential_t>::emtensor_excl_potential(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &gamma_UU, const Tensor<3, data_t> &chris_phys_ULL)
{
    // Useful quantity Vt
    data_t Vt = -vars.Pi_Re * vars.Pi_Re - vars.Pi_Im * vars.Pi_Im;
    FOR2(i, j)
    {
        Vt += gamma_UU[i][j] * d1.phi_Re[i] * d1.phi_Re[j] +
              gamma_UU[i][j] * d1.phi_Im[i] * d1.phi_Im[j];
    }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] = -0.5 * metric_vars.gamma[i][j] * Vt +
                        d1.phi_Re[i] * d1.phi_Re[j] +
                        d1.phi_Im[i] * d1.phi_Im[j];
    }

    // S = Tr_S_ij
    out.S = TensorAlgebra::compute_trace(out.Sij, gamma_UU);

    // S_i (note lower index) = - n^a T_a0
    FOR1(i)
    {
        out.Si[i] = -d1.phi_Re[i] * vars.Pi_Re - d1.phi_Im[i] * vars.Pi_Im;
    }

    // rho = n^a n^b T_ab
    out.rho = vars.Pi_Re * vars.Pi_Re + vars.Pi_Im * vars.Pi_Im + 0.5 * Vt;
}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FixedBGComplexScalarField<potential_t>::matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // call the function for the rhs excluding the potential
    matter_rhs_excl_potential(total_rhs, vars, metric_vars, d1, d2, advec);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi_Re = 0.0;
    data_t dVdphi_Im = 0.0;

    // compute potential
    my_potential.compute_potential(V_of_phi, dVdphi_Re, dVdphi_Im, vars);

    // adjust RHS for the potential term
    total_rhs.Pi_Re += -metric_vars.lapse * dVdphi_Re;
    total_rhs.Pi_Im += -metric_vars.lapse * dVdphi_Im;
}

// the RHS excluding the potential terms
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FixedBGComplexScalarField<potential_t>::matter_rhs_excl_potential(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec)
{
    using namespace TensorAlgebra;

    const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
    const auto chris_phys = compute_christoffel(metric_vars.d1_gamma, gamma_UU);

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs.phi_Re = metric_vars.lapse * vars.Pi_Re + advec.phi_Re;
    rhs.phi_Im = metric_vars.lapse * vars.Pi_Im + advec.phi_Im;
    rhs.Pi_Re = metric_vars.lapse * metric_vars.K * vars.Pi_Re + advec.Pi_Re;
    rhs.Pi_Im = metric_vars.lapse * metric_vars.K * vars.Pi_Im + advec.Pi_Im;

    FOR2(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs.Pi_Re += gamma_UU[i][j] * (metric_vars.lapse * d2.phi_Re[i][j] +
                                       metric_vars.d1_lapse[i] * d1.phi_Re[j]);
        rhs.Pi_Im += gamma_UU[i][j] * (metric_vars.lapse * d2.phi_Im[i][j] +
                                       metric_vars.d1_lapse[i] * d1.phi_Im[j]);

        FOR1(k)
        {
            rhs.Pi_Re += -metric_vars.lapse * gamma_UU[i][j] *
                         (chris_phys.ULL[k][i][j] * d1.phi_Re[k]);
            rhs.Pi_Im += -metric_vars.lapse * gamma_UU[i][j] *
                         (chris_phys.ULL[k][i][j] * d1.phi_Im[k]);
        }
    }
}

#endif /* FIXEDBGCOMPLEXSCALARFIELD_IMPL_HPP_ */

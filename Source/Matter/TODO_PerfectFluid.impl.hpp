/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(PERFECTFLUID_HPP_)
#error "This file should only be included through PerfectFluid.hpp"
#endif

#ifndef PERFECTFLUID_IMPL_HPP_
#define PERFECTFLUID_IMPL_HPP_
#include "DimensionDefinitions.hpp"

// Calculate the stress energy tensor elements
template <class eos_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> PerfectFluid<eos_t>::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;

    // Copy the field vars into FluidObject
    FluidObject<data_t> vars_fl;
    vars_fl.density = vars.density;
    vars_fl.energy = vars.energy;                                                     //FIXME: check that I should be using vars_fl, and not vars

    // set the eos values
    data_t pressure = 0.0;   // P = P ( density, energy)
    data_t enthalpy = 0.0;   // h = 1 + energy + pressure/density

    // compute potential and add constributions to EM Tensor
    my_eos.compute_eos(pressure, enthalpy, vars);

    { //emtensor_excl_potential  TODO:Is it okay here, or create a function?

        FOR2(i, j) { Vt += vars.chi * h_UU[i][j] * d1_phi[i] * d1_phi[j]; }           //FIXME: check that I should be using vars, and not vars_fl

        // Calculate components of EM Tensor
        // S_ij = T_ij
        FOR2(i, j)
        {
            out.Sij[i][j] =
              vars.density * vars.u[i] * vars.u[j] +
              pressure * vars.h[i][j] / vars.chi;
        }

        // S_i (note lower index) = - n^a T_ai
        FOR1(i) { out.Si[i] = vars.lapse * vars.density * enthalpy *
                    vars.u[i] * vars.u[0]; }

        // S = Tr_S_ij
        out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);


        // rho = n^a n^b T_ab
        out.rho = vars.lapse * vars.lapse *
                 (vars.density * enthalpy * vars.u[0] * vars.u[0] -
                  pressure) ;
    }

    return out;
}

// // Calculate the stress energy tensor elements
// template <class potential_t>
// template <class data_t, template <typename> class vars_t>
// void ScalarField<potential_t>::emtensor_excl_potential(
//     emtensor_t<data_t> &out, const vars_t<data_t> &vars,
//     const SFObject<data_t> &vars_sf, const Tensor<1, data_t> &d1_phi,
//     const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL)
// {
//     // Useful quantity Vt
//     data_t Vt = -vars_sf.Pi * vars_sf.Pi;
//     FOR2(i, j) { Vt += vars.chi * h_UU[i][j] * d1_phi[i] * d1_phi[j]; }
//
//     // Calculate components of EM Tensor
//     // S_ij = T_ij
//     FOR2(i, j)
//     {
//         out.Sij[i][j] =
//             -0.5 * vars.h[i][j] * Vt / vars.chi + d1_phi[i] * d1_phi[j];
//     }
//
//     // S = Tr_S_ij
//     out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);
//
//     // S_i (note lower index) = - n^a T_ai
//     FOR1(i) { out.Si[i] = -d1_phi[i] * vars_sf.Pi; }
//
//     // rho = n^a n^b T_ab
//     out.rho = vars_sf.Pi * vars_sf.Pi + 0.5 * Vt;
// }


// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void PerfectFluid<eos_t>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // first get the non potential part of the rhs
    // this may seem a bit long winded, but it makes the function
    // work for more multiple fields

    // the rhs vars
    FluidObject<data_t> rhs_fl;
    rhs_fl.W = 0;
    rhs_fl.D = 0;
    rhs_fl.E = 0;
    rhs_fl.Z[0] = 0;

    // advection terms
    FluidObject<data_t> advec_fl;

    // the vars
    FluidObject<data_t> vars_fl;
    vars_fl.W = vars.W;
    vars_fl.D = vars.D;
    vars_fl.E = vars.E;
    vars_fl.Z[0] = vars.Z[0];

    FOR1(i) {
      advec_fl.Z[i] = advec.Z[i];

      vars_fl.V[i] = vars.V[i];
      vars_fl.Z[i] = vars.Z[i];

      rhs_fl.V[i] = 0;
      rhs_fl.Z[i] = 0;
    }

    // set the eos values
    data_t pressure = 0.0;   // P = P ( density, energy)
    data_t enthalpy = 0.0;   // h = 1 + energy + pressure/density

    // compute potential and add constributions to EM Tensor
    my_eos.compute_eos(pressure, enthalpy, vars);

    // useful variable
    data_t dt_W = 0.0;                                                                        // FIXME: need to be defined

    {  //matter_rhs_excl_potential  //TODO: create indp function?
      using namespace TensorAlgebra;

      const auto h_UU = compute_inverse_sym(vars.h);
      const auto chris = compute_christoffel(d1.h, h_UU);

      // evolution equations for scalar field and (minus) its conjugate momentum
      FOR1(i)
      {
          rhs_sf.D += - d1.D[i] * vars_fl.V[i] - vars_fl.D * d1.V[i][i];
          rhs_sf.E += - d1.E[i] * vars_fl.V[i] - vars_fl.D * d1.V[i][i] +
                   ( - d1.W[i] * vars_fl.V[i] - vars_fl.W * d1.V[i][i] - dt_W)
                   * pressure;

          rhs_sf.Z[i] +=  vars_fl.Z[0] * vars.lapse * d1.lapse[i] +
                      - vars.lapse * pow(vars.chi, 1.5) * d1.preassure[i]                 // FIXME:   first deriv. of pressure!!  var not exist yet
                      +  advec_sf.Z[i];
      }

      FOR2(i, j)
      {
          // includes non conformal parts of chris not included in chris_ULL
          rhs_sf.Z[i] += - d1.V[j][j] * vars_fl.Z[i] - d1.Z[j][i] * vars_fl.V[j];


          FOR1(k)
          {
              rhs_sf.Z[i] +=  0.5 * vars_fl.Z[j] * vars_fl.Z[k] *
                              d1.h_UU[i][j][k] / vars_fl.Z[0];                            // FIXME:  d1.h_UU  var not exist yet
          }
      }
  }
}

//                                                                                        //FIXME: Stopped  coding here!

// // the RHS excluding the potential terms
// template <class potential_t>
// template <class data_t, template <typename> class vars_t>
// void ScalarField<potential_t>::matter_rhs_excl_potential(
//     SFObject<data_t> &rhs_sf, const vars_t<data_t> &vars,
//     const SFObject<data_t> &vars_sf, const vars_t<Tensor<1, data_t>> &d1,
//     const Tensor<1, data_t> &d1_phi, const Tensor<2, data_t> &d2_phi,
//     const SFObject<data_t> &advec_sf)
// {
//     using namespace TensorAlgebra;
//
//     const auto h_UU = compute_inverse_sym(vars.h);
//     const auto chris = compute_christoffel(d1.h, h_UU);
//
//     // evolution equations for scalar field and (minus) its conjugate momentum
//     rhs_sf.phi = vars.lapse * vars_sf.Pi + advec_sf.phi;
//     rhs_sf.Pi = vars.lapse * vars.K * vars_sf.Pi + advec_sf.Pi;
//
//     FOR2(i, j)
//     {
//         // includes non conformal parts of chris not included in chris_ULL
//         rhs_sf.Pi += h_UU[i][j] * (-0.5 * d1.chi[j] * vars.lapse * d1_phi[i] +
//                                    vars.chi * vars.lapse * d2_phi[i][j] +
//                                    vars.chi * d1.lapse[i] * d1_phi[j]);
//         FOR1(k)
//         {
//             rhs_sf.Pi += -vars.chi * vars.lapse * h_UU[i][j] *
//                          chris.ULL[k][i][j] * d1_phi[k];
//         }
//     }
// }

#endif /* PERFECTFLUID_IMPL_HPP_ */

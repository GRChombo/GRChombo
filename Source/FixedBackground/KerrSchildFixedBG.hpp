/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef KERRSCHILDFIXEDBG_HPP_
#define KERRSCHILDFIXEDBG_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the initial conditions for a Kerr Schild BH
//! https://arxiv.org/pdf/gr-qc/9805023.pdf
class KerrSchildFixedBG
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double mass = 1.0;                      //!<< The mass of the BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the BH
        double spin = 0.0;                      //!< The spin param a = J / M
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

    const params_t m_params;
    const double m_dx;

    KerrSchildFixedBG(params_t a_params, double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
        // check this spin param is sensible
        if ((m_params.spin > m_params.mass) || (m_params.spin < -m_params.mass))
        {
            MayDay::Error(
                "The dimensionless spin parameter must be in the range "
                "-1.0 < spin < 1.0");
        }
    }

    /// This just calculates chi which helps with regridding, debug etc
    /// it is only done once on setup as the BG is fixed
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // get position and set vars
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        Vars<data_t> metric_vars;
        compute_metric_background(metric_vars, current_cell);

        // calculate and save chi
        data_t chi = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        chi = pow(chi, -1.0 / 3.0);
        current_cell.store_vars(chi, c_chi);
    }

    // Kerr Schild solution
    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Cell<data_t> &current_cell) const
    {
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

        // black hole params - mass M and spin a
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid including effect of spin
        // on x direction (length contraction)
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t rho = coords.get_radius();
        const data_t rho2 = rho * rho;

        // the Kerr Schild radius r
        const data_t r2 = 0.5 * (rho2 - a2) +
                          sqrt(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z);
        const data_t r = sqrt(r2);
        const data_t cos_theta = z / r;

        // find the H and el quantities (el decomposed into space and time)
        data_t H = M * r / (r2 + a2 * cos_theta * cos_theta);
        const Tensor<1, data_t> el = {(r * x + a * y) / (r2 + a2),
                                      (r * y - a * x) / (r2 + a2), z / r};
        const data_t el_t = 1.0;

        // Calculate the gradients in el and H
        Tensor<1, data_t> dHdx;
        Tensor<1, data_t> dltdx;
        Tensor<2, data_t> dldx;
        get_KS_derivs(dHdx, dldx, dltdx, H, coords);

        // populate ADM vars
        vars.lapse = pow(1.0 + 2.0 * H * el_t * el_t, -0.5);
        FOR2(i, j)
        {
            vars.gamma[i][j] =
                TensorAlgebra::delta(i, j) + 2.0 * H * el[i] * el[j];
        }
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(vars.gamma);
        FOR1(i)
        {
            vars.shift[i] = 0;
            FOR1(j)
            {
                vars.shift[i] += gamma_UU[i][j] * 2.0 * H * el[j] * el_t;
            }
        }

        // Calculate partial derivative of spatial metric
        FOR3(i, j, k)
        {
            vars.d1_gamma[i][j][k] =
                2.0 * (el[i] * el[j] * dHdx[k] + H * el[i] * dldx[j][k] +
                       H * el[j] * dldx[i][k]);
        }

        // calculate derivs of lapse and shift
        FOR1(i)
        {
            vars.d1_lapse[i] = -pow(vars.lapse, 3.0) * el_t *
                               (el_t * dHdx[i] + 2.0 * H * dltdx[i]);
        }

        // use the fact that shift^i = lapse^2 * shift_i
        FOR2(i, j)
        {
            vars.d1_shift[i][j] =
                2.0 * el_t * dHdx[j] * pow(vars.lapse, 2.0) * el[i] +
                4.0 * el_t * H * vars.lapse * vars.d1_lapse[j] * el[i] +
                2.0 * el_t * H * pow(vars.lapse, 2.0) * dldx[i][j] +
                2.0 * dltdx[j] * H * pow(vars.lapse, 2.0) * el[i];
        }

        // calculate the extrinsic curvature, using the fact that
        // 2 * lapse * K_ij = D_i \beta_j + D_j \beta_i - dgamma_ij dt
        // and dgamma_ij dt = 0 in chosen fixed gauge
        const auto chris_phys = compute_christoffel(vars.d1_gamma, gamma_UU);
        FOR2(i, j)
        {
            vars.K_tensor[i][j] = 0.0;
            FOR1(k)
            {
                vars.K_tensor[i][j] +=
                    vars.gamma[k][j] * vars.d1_shift[k][i] +
                    vars.gamma[k][i] * vars.d1_shift[k][j] +
                    (vars.d1_gamma[k][i][j] + vars.d1_gamma[k][j][i]) *
                        vars.shift[k];
                FOR1(m)
                {
                    vars.K_tensor[i][j] += -2.0 * chris_phys.ULL[k][i][j] *
                                           vars.gamma[k][m] * vars.shift[m];
                }
            }
            vars.K_tensor[i][j] *= 0.5 / vars.lapse;
        }
        vars.K = compute_trace(gamma_UU, vars.K_tensor);
    }

  protected:
    /// Work out the gradients of the quantities H and el appearing in the Kerr
    /// Schild solution
    template <class data_t>
    void get_KS_derivs(Tensor<1, data_t> &dHdx, Tensor<2, data_t> &dldx,
                       Tensor<1, data_t> &dltdx, const data_t &H,
                       const Coordinates<data_t> &coords) const
    {
        // black hole params - mass M and boost v
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid, and useful quantities
        Tensor<1, data_t> x;
        x[0] = coords.x;
        x[1] = coords.y;
        x[2] = coords.z;
        const double z = coords.z;
        const data_t rho = coords.get_radius();
        const data_t rho2 = rho * rho;

        // the Kerr Schild radius r
        const data_t r2 = 0.5 * (rho2 - a2) +
                          sqrt(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z);
        const data_t r = sqrt(r2);
        const data_t cos_theta = z / r;
        const data_t cos_theta2 = cos_theta * cos_theta;

        using namespace TensorAlgebra;
        // derivatives of r wrt actual grid coords
        Tensor<1, data_t> drhodx;
        FOR1(i) { drhodx[i] = x[i] / rho; }

        Tensor<1, data_t> drdx;
        FOR1(i)
        {
            drdx[i] =
                0.5 / r *
                (rho * drhodx[i] +
                 0.5 / sqrt(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z) *
                     (drhodx[i] * rho * (rho2 - a2) +
                      delta(i, 2) * 2.0 * a2 * z));
        }

        Tensor<1, data_t> dcosthetadx;
        FOR1(i) { dcosthetadx[i] = -z / r2 * drdx[i] + delta(i, 2) / r; }

        FOR1(i)
        {
            dHdx[i] = H * (drdx[i] / r -
                           2.0 / (r2 + a2 * cos_theta2) *
                               (r * drdx[i] + a2 * cos_theta * dcosthetadx[i]));
        }

        // note to use convention as in rest of tensors the last index is the
        // derivative index so these are d_i l_j
        FOR1(i)
        {
            // first the el_x comp
            dldx[0][i] =
                (x[0] * drdx[i] + r * delta(i, 0) + a * delta(i, 1) -
                 2.0 * r * drdx[i] * (r * x[0] + a * x[1]) / (r2 + a2)) /
                (r2 + a2);
            // now the el_y comp
            dldx[1][i] =
                (x[1] * drdx[i] + r * delta(i, 1) - a * delta(i, 0) -
                 2.0 * r * drdx[i] * (r * x[1] - a * x[0]) / (r2 + a2)) /
                (r2 + a2);
            // now the el_z comp
            dldx[2][i] = -x[2] * drdx[i] / r2 + delta(i, 2) / r;
        }

        // then dltdi
        FOR1(i) { dltdx[i] = 0.0; }
    }

  public:
    // used to decide when to excise - ie when within the horizon of the BH
    // note that this is not templated over data_t
    double excise(const Cell<double> &current_cell) const
    {
        // black hole params - mass M and spin a
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid
        const Coordinates<double> coords(current_cell, m_dx, m_params.center);
        const double x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const double r_plus = M + sqrt(M * M - a2);
        const double r_minus = M - sqrt(M * M - a2);

        const double outer_horizon =
            (x * x + y * y) / (2.0 * M * r_plus) + z * z / r_plus / r_plus;

        const double inner_horizon =
            (x * x + y * y) / (2.0 * M * r_minus) + z * z / r_minus / r_minus;

        // value less than 1 indicates we are within the horizon
        return sqrt(outer_horizon);
    }
};

#endif /* KERRSCHILDFIXEDBG_HPP_ */

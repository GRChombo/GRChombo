/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(TWOPUNCTURESINITIALDATA_HPP_)
#error "This file should only be included through TwoPuncturesInitialData.hpp"
#endif

#ifndef TWOPUNCTURESINITIALDATA_IMPL_HPP_
#define TWOPUNCTURESINITIALDATA_IMPL_HPP_

void TwoPuncturesInitialData::compute(Cell<double> current_cell) const
{
    Vars<double> vars;
    // Set only the non-zero components explicitly below
    VarsTools::assign(vars, 0.);

    Coordinates<double> coords(current_cell, m_dx, m_center);
    Tensor<2, double> h_phys, extrinsic_K;
    Tensor<1, double> shift, Z3;
    double lapse, Theta;

    interpolate_tp_vars(coords, h_phys, extrinsic_K, lapse, shift, Theta, Z3);

    // analytically set Bowen-York properties below (e.g. conformal flatness,
    // tracefree K)

    // metric variables
    vars.chi = 1.0 / h_phys[0][0];
    FOR1(i) { vars.h[i][i] = 1.0; }

    // extrinsic curvature
    FOR2(i, j) { vars.A[i][j] = vars.chi * extrinsic_K[i][j]; }
    // conformal flatness means h_UU = h
    TensorAlgebra::make_trace_free(vars.A, vars.h, vars.h);

    // gauge
    vars.lapse = lapse;

    current_cell.store_vars(vars);
}

void TwoPuncturesInitialData::interpolate_tp_vars(
    const Coordinates<double> &coords, Tensor<2, double> &out_h_phys,
    Tensor<2, double> &out_extrinsic_K, double &out_lapse,
    Tensor<1, double> &out_shift, double &out_Theta,
    Tensor<1, double> &out_Z3) const
{
    double coords_array[CH_SPACEDIM];
    coords_array[0] = coords.x;
    coords_array[1] = coords.y;
    coords_array[2] = coords.z;

    using namespace TP::Z4VectorShortcuts;
    double TP_state[Qlen];
    m_two_punctures.Interpolate(coords_array, TP_state);

    // metric
    out_h_phys[0][0] = TP_state[g11];
    out_h_phys[0][1] = out_h_phys[1][0] = TP_state[g12];
    out_h_phys[0][2] = out_h_phys[2][0] = TP_state[g13];
    out_h_phys[1][1] = TP_state[g22];
    out_h_phys[1][2] = out_h_phys[2][1] = TP_state[g23];
    out_h_phys[2][2] = TP_state[g33];

    // extrinsic curvature
    out_extrinsic_K[0][0] = TP_state[K11];
    out_extrinsic_K[0][1] = out_extrinsic_K[1][0] = TP_state[K12];
    out_extrinsic_K[0][2] = out_extrinsic_K[2][0] = TP_state[K13];
    out_extrinsic_K[1][1] = TP_state[K22];
    out_extrinsic_K[1][2] = out_extrinsic_K[2][1] = TP_state[K23];
    out_extrinsic_K[2][2] = TP_state[K33];

    // Z4 vector
    out_Z3[0] = TP_state[Z1];
    out_Z3[1] = TP_state[Z2];
    out_Z3[2] = TP_state[Z3];
    out_Theta = TP_state[Theta];

    // gauge
    out_lapse = TP_state[lapse];
    out_shift[0] = TP_state[shift1];
    out_shift[1] = TP_state[shift2];
    out_shift[2] = TP_state[shift3];
}

#endif /* TWOPUNCTURESINITIALDATA_IMPL_HPP_ */

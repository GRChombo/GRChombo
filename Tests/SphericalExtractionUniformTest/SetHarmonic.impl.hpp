/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SETHARMONIC_HPP_)
#error "This file should only be included through SetHarmonic.hpp"
#endif

#ifndef SETHARMONIC_IMPL_HPP_
#define SETHARMONIC_IMPL_HPP_

#include "DebuggingTools.hpp"
#include "SetHarmonic.hpp"
#include "SphericalHarmonics.hpp"
#include "simd.hpp"

template <class data_t>
void SetHarmonic::compute(Cell<data_t> current_cell) const
{

    Coordinates<data_t> coords(current_cell, m_dx, m_center);

    using namespace SphericalHarmonics;
    auto Y_lm = spin_Y_lm(coords.x, coords.y, coords.z, m_es, m_el, m_em);
    data_t out_Re = Y_lm.Real;
    data_t out_Im = Y_lm.Im;

    current_cell.store_vars(out_Re, m_var_Re);
    current_cell.store_vars(out_Im, m_var_Im);
}

#endif /* SETHARMONIC_IMPL_HPP_ */

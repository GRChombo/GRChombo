/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETHARMONIC_HPP_
#define SETHARMONIC_HPP_

#include "BoxLoops.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

// This compute class sets two vars to the real and imaginary parts of the
// s = a_es, l = a_el, m = a_em spin-weighted spherical harmonic
class SetHarmonic
{
  public:
    SetHarmonic(int a_var_Re, int a_var_Im, int a_es, int a_el, int a_em,
                const std::array<double, CH_SPACEDIM> &a_center, double a_dx)
        : m_var_Re(a_var_Re), m_var_Im(a_var_Im), m_es(a_es), m_el(a_el),
          m_em(a_em), m_dx(a_dx), m_center(a_center)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const int m_var_Re;
    const int m_var_Im;
    const int m_es;
    const int m_el;
    const int m_em;
    const double m_dx;
    const std::array<double, CH_SPACEDIM> &m_center;
};

#include "SetHarmonic.impl.hpp"

#endif /* SETHARMONIC_HPP_ */

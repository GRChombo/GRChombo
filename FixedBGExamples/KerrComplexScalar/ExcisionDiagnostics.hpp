/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONDIAGNOSTICS_HPP_
#define EXCISIONDIAGNOSTICS_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Does excision for fixed BG BH solutions
//! Note that it is does not using simd so one must set disable_simd()
template <class matter_t, class background_t> class ExcisionDiagnostics
{
  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center
    const FourthOrderDerivatives m_deriv;
    const background_t m_background;
    const double m_inner_r;
    const double m_outer_r;

  public:
    ExcisionDiagnostics(const double a_dx,
                        const std::array<double, CH_SPACEDIM> a_center,
                        background_t a_background, const double a_inner_r,
                        const double a_outer_r)
        : m_dx(a_dx), m_deriv(m_dx), m_center(a_center),
          m_background(a_background), m_inner_r(a_inner_r), m_outer_r(a_outer_r)
    {
    }

    void compute(const Cell<double> current_cell) const
    {
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        double horizon_distance = m_background.excise(current_cell);
        if (coords.get_radius() < m_inner_r || coords.get_radius() > m_outer_r)
        {
            current_cell.store_vars(0.0, c_rho);
            current_cell.store_vars(0.0, c_rhoJ);
        } // else do nothing
    }
};

#endif /* EXCISIONDIAGNOSTICS_HPP_ */

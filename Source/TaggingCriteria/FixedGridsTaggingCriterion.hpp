/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDGRIDSTAGGINGCRITERION_HPP_
#define FIXEDGRIDSTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class FixedGridsTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const FourthOrderDerivatives m_deriv;
    const int m_level;
    const std::array<double, CH_SPACEDIM> m_center;

  public:
    FixedGridsTaggingCriterion(const double dx, const int a_level,
                               const double a_L,
                               const std::array<double, CH_SPACEDIM> a_center)
        : m_dx(dx), m_deriv(dx), m_level(a_level), m_L(a_L),
          m_center(a_center){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        data_t criterion = 0.0;
        // make sure the inner part is regridded around the horizon
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4
        double ratio = pow(2.0, -(m_level + 2.0));
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const data_t max_abs_xy = simd_max(abs(coords.x), abs(coords.y));
        const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));
        auto regrid = simd_compare_lt(max_abs_xyz, m_L * ratio);
        criterion = simd_conditional(regrid, 100.0, criterion);

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* FIXEDGRIDSTAGGINGCRITERION_HPP_ */

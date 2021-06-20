/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTEGRATIONMETHOD_HPP_
#define INTEGRATIONMETHOD_HPP_

// Chombo includes
#include "CH_assert.H"

// Other includes
#include <utility>
#include <vector>

// Chombo namespace
#include "UsingNamespace.H"

//! A class to store and return the weights associated to a Newton-Cotes formula
//! for numerical integration/quadrature which can be closed (i.e. includes the
//! endpoints) or open (does not include the end points). This is used by
//! SurfaceExtraction for integration over extraction surfaces.
class IntegrationMethod
{
  private:
    std::vector<double> m_weights;
    int m_num_weights;
    bool m_is_closed;

  public:
    //! Constructor
    IntegrationMethod(const std::vector<double> &a_weights,
                      bool a_is_closed = true)
        : m_weights(a_weights), m_num_weights(a_weights.size()),
          m_is_closed(a_is_closed)
    {
        CH_assert(a_weights.size() > 0);
    }

    //! Checks that this integration method is suitable given the number of
    //! points and periodicity
    inline bool is_valid(int a_num_points, bool a_is_periodic) const
    {
        if (m_is_closed && !a_is_periodic)
        {
            return (a_num_points % m_num_weights == 1 || m_num_weights == 1);
        }
        else
        {
            return (a_num_points % m_num_weights == 0);
        }
    }

    //! Returns whether this IntegrationMethod is closed or not
    inline bool is_closed() const { return m_is_closed; }

    //! Returns the weight for a point with given index
    inline double weight(int a_index, int a_num_points,
                         bool a_is_periodic) const
    {
        const int weight_index = a_index % m_num_weights;
        const bool endpoint =
            (a_index == 0 || a_index == a_num_points - 1) && !a_is_periodic;
        // if this is a closed formula, not a geometry endpoint but at the edge
        // of the formula, need to double the weight as this is how Newton-Cotes
        // formulae are combined.
        if (m_is_closed && !endpoint && weight_index == 0)
            return 2.0 * m_weights[weight_index];
        else // otherwise we just use the weight from the formula
            return m_weights[weight_index];
    }

    // Closed default methods
    static const IntegrationMethod trapezium;
    static const IntegrationMethod simpson;
    static const IntegrationMethod simpson38;
    static const IntegrationMethod boole;

    // Open default methods
    static const IntegrationMethod midpoint;
    static const IntegrationMethod open_2nd_order;
    static const IntegrationMethod open_3rd_order;
    static const IntegrationMethod open_4th_order;
};

#endif /* INTEGRATIONMETHOD_HPP_ */

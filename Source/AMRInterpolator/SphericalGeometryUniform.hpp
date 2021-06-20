/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHERICALGEOMETRYUNIFORM_HPP_
#define SPHERICALGEOMETRYUNIFORM_HPP_

// Chombo includes
#include "AlwaysInline.hpp"
#include "MayDay.H"

// Other includes
#include "AlwaysInline.hpp"
#include "IntegrationMethod.hpp"
#include "MayDay.H"
#include "SphericalGeometry.hpp"
#include <array>
#include <cmath>
#include <string>

// Chombo namespace
#include "UsingNamespace.H"

//! This SurfaceGeometry template class provides spherical shell geometry
//! implementation for the SurfaceExtraction class
//! u = u, v = phi
class SphericalGeometryUniform
{
  public:
    SphericalGeometryUniform(const std::array<double, CH_SPACEDIM> &a_center)
        : m_center(a_center)
    {
#if CH_SPACEDIM == 3 // for now force 'Z', could be an argument later on
        m_up_dir = SphericalGeometry::Z;
#elif CH_SPACEDIM == 2 // in Cartoon method, assume 'X' as 'z' axis
        m_up_dir = SphericalGeometry::X;
#endif
    }

    ALWAYS_INLINE SphericalGeometry::UP_DIR get_up_dir() const
    {
        return m_up_dir;
    }

    static ALWAYS_INLINE double get_domain_u_min() { return -1.; }
    static ALWAYS_INLINE double get_domain_u_max() { return 1.; }
    static ALWAYS_INLINE double get_domain_v_min() { return 0.; }
    static ALWAYS_INLINE double get_domain_v_max() { return 2. * M_PI; }

    //! returns the grid spacing in u
    static ALWAYS_INLINE double du(int a_num_points_u)
    {
        return (get_domain_u_max() - get_domain_u_min()) /
               (double)(a_num_points_u);
    }

    //! returns the grid spacing in phi
    static ALWAYS_INLINE double dv(int a_num_points_phi)
    {
        return (get_domain_v_max() - get_domain_v_min()) /
               ((double)a_num_points_phi);
    }

    //! returns the u coordinate associated to the u/u index
    static ALWAYS_INLINE double u(int a_iu, int a_num_points_u)
    {
        return (a_iu + 0.5) * du(a_num_points_u) + get_domain_u_min();
    }

    //! returns the phi coordinate associated to the phi/v index
    static ALWAYS_INLINE double v(int a_iphi, int a_num_points_phi)
    {
        return a_iphi * dv(a_num_points_phi) + get_domain_v_min();
    }

    static ALWAYS_INLINE bool is_u_periodic() { return false; }
    static ALWAYS_INLINE bool is_v_periodic() { return true; }

    static ALWAYS_INLINE std::string param_name() { return "r"; }

    static ALWAYS_INLINE std::string u_name() { return "z/r"; }

    static ALWAYS_INLINE std::string v_name() { return "phi"; }

    //! returns the area element on a sphere with radius a_radius at the point
    //! (a_u, a_phi)
    static ALWAYS_INLINE double area_element(double a_radius, double a_u,
                                             double a_phi)
    {
        return a_radius * a_radius;
    }

    //! returns the Cartesian coordinate in direction a_dir with specified
    //! radius, u and phi.
    ALWAYS_INLINE double get_grid_coord(int a_dir, double a_radius, double a_u,
                                        double a_phi = 0.) const
    {
        double center = (a_dir < CH_SPACEDIM ? m_center[a_dir] : 0.);
        switch ((a_dir + 2 - m_up_dir) % 3)
        {
        case (0):
            return center + a_radius * sqrt(1. - a_u * a_u) * cos(a_phi);
        case (1):
            return center + a_radius * sqrt(1. - a_u * a_u) * sin(a_phi);
        case (2):
            return center + a_radius * a_u;
        default:
            MayDay::Error("SphericalGeometryUniform: Direction not supported");
        }
    }

    static const IntegrationMethod &
    get_recommended_integration_method_u(int a_num_points_u)
    {
        static const IntegrationMethod &open_2nd_order =
            IntegrationMethod::open_2nd_order;
        static const IntegrationMethod &midpoint = IntegrationMethod::midpoint;
        if (open_2nd_order.is_valid(a_num_points_u, is_u_periodic()))
            return open_2nd_order;
        MayDay::Warning(
            "Use a multipel of 3 in 'u' points to use milne's rule. "
            "Defaulting to midpoint.");
        return midpoint;
    }

    static const IntegrationMethod &
    get_recommended_integration_method_v(int a_num_points_phi)
    {
        static const IntegrationMethod &simpson = IntegrationMethod::simpson;
        static const IntegrationMethod &trapezium =
            IntegrationMethod::trapezium;
        if (simpson.is_valid(a_num_points_phi, is_v_periodic()))
            return simpson;
        MayDay::Warning("Use an even number of 'v' points to use simpson rule. "
                        "Defaulting to trapezium.");
        return trapezium;
    }

  protected:
    std::array<double, CH_SPACEDIM> m_center;
    SphericalGeometry::UP_DIR m_up_dir;
};

#endif /* SPHERICALGEOMETRYUNIFORM_HPP_ */

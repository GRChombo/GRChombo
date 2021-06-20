/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CYLINDRICALGEOMETRY_HPP_
#define CYLINDRICALGEOMETRY_HPP_

// Chombo includes
#include "MayDay.H"

// Other includes
#include <array>
#include <cmath>
#include <string>

// Chombo namespace
#include "UsingNamespace.H"

//! This class, derived from the SurfaceGeometry class, provides cylindrical
//! shell geometry implementation for the SurfaceExtraction class v = phi
class CylindricalGeometry
{
  private:
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_z_length;

  public:
    CylindricalGeometry(const std::array<double, CH_SPACEDIM> &a_center,
                        double a_z_length)
        : m_center(a_center), m_z_length(a_z_length)
    {
    }

    ALWAYS_INLINE double get_domain_u_min() const { return 0.; }
    ALWAYS_INLINE double get_domain_u_max() const { return m_z_length; }
    ALWAYS_INLINE double get_domain_v_min() const { return 0.; }
    ALWAYS_INLINE double get_domain_v_max() const { return 2. * M_PI; }

    //! returns the grid spacing in z
    ALWAYS_INLINE double du(int a_num_points_z) const
    {
        return (get_domain_u_max() - get_domain_u_min()) / (a_num_points_z - 1);
    }

    //! returns the grid spacing in phi
    ALWAYS_INLINE double dv(int a_num_points_phi) const
    {
        return (get_domain_v_max() - get_domain_v_min()) /
               ((double)a_num_points_phi);
    }

    //! returns the z coordinate associated to the theta/u index
    ALWAYS_INLINE double u(int a_iz, int a_num_points_z) const
    {
        return a_iz * du(a_num_points_z) -
               ((get_domain_v_max() - get_domain_v_min()) / 2.0);
    }

    //! returns the phi coordinate associated to the phi/v index
    ALWAYS_INLINE double v(int a_iphi, int a_num_points_phi) const
    {
        return a_iphi * dv(a_num_points_phi);
    }

    static ALWAYS_INLINE bool is_u_periodic() { return false; }
    static ALWAYS_INLINE bool is_v_periodic() { return true; }
    //! returns the Cartesian coordinate in direction a_dir with specified
    //! radius, z and phi.
    ALWAYS_INLINE double get_grid_coord(int a_dir, double a_radius, double a_z,
                                        double a_phi) const
    {
        switch (a_dir)
        {
        case (0):
            return m_center[0] + a_radius * cos(a_phi);
        case (1):
            return m_center[1] + a_radius * sin(a_phi);
        case (2):
            return m_center[2] + a_z;
        default:
            MayDay::Error("CylindricalGeometry: Direction not supported");
        }
    }

    //! returns the area element on a cylinder with radius a_radius at the point
    //! (a_z, a_phi)
    static ALWAYS_INLINE double area_element(double a_radius, double a_z,
                                             double a_phi)
    {
        return a_radius;
    }

    static ALWAYS_INLINE std::string param_name() { return "r"; }
    static ALWAYS_INLINE std::string u_name() { return "z"; }
    static ALWAYS_INLINE std::string v_name() { return "phi"; }
};

#endif /* CYLINDRICALGEOMETRY_HPP_ */

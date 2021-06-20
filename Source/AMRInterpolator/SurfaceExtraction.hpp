/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SURFACEEXTRACTION_HPP_
#define SURFACEEXTRACTION_HPP_

// Chombo includes
#include "CH_assert.H"
#include "SPMD.H"

// Other inclues
#include "AMRInterpolator.hpp"
#include "DimensionDefinitions.hpp"
#include "FilesystemTools.hpp"
#include "IntegrationMethod.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "SmallDataIO.hpp" // for writing data
#include "UserVariables.hpp"
#include <algorithm>
#include <array>
#include <functional>
#include <tuple>
#include <utility>
#include <vector>

// Chombo namespace
#include "UsingNamespace.H"

struct surface_extraction_params_t
{
    int num_surfaces; //!< number of surfaces over which to extraction
    std::vector<double>
        surface_param_values; //!< the values of the
                              //!< parameter that gives the required
                              //!< surfaces with SurfaceGeom geometry (e.g.
                              //!< radii for spherical shells)
    int num_points_u;         //!< the number of points for the first parameter
                              //!< that parameterises each surface
    int num_points_v;         //!< the number of points for the second parameter
                              //!< that parameterises each surfaces
    std::vector<int> extraction_levels; //!< the level on which to do the
                                        //!< extraction for each surface
    bool write_extraction; //!< whether or not to write the extracted data

    std::string data_path, integral_file_prefix;
    std::string extraction_path, extraction_file_prefix;

    int min_extraction_level() const
    {
        return *(std::min_element(extraction_levels.begin(),
                                  extraction_levels.end()));
    }
};

//! This class extracts grid variables on 2 dimensional surfaces each
//! parameterised by u and v with different surfaces given by level sets of
//! another parameter
template <class SurfaceGeometry> class SurfaceExtraction
{
  public:
    using vars_t = std::tuple<int, VariableType, Derivative>;

  protected:
    const SurfaceGeometry m_geom; //!< the geometry class which knows about
                                  //!< the particular surface
    const surface_extraction_params_t m_params;
    std::vector<std::tuple<int, VariableType, Derivative>>
        m_vars; //!< the vector of pairs of
    //!< variables and derivatives to extract
    const double m_dt;
    const double m_time;
    const bool m_first_step;
    const double m_restart_time;
    const int m_num_interp_points; //!< the total number of points this rank
                                   //!< will extract (0 on ranks > 0)
    const double m_du; //!< the grid spacing in u (used in integrate)
    const double m_dv; //!< the grid spacing in v (used in integrate)

    std::vector<std::vector<double>> m_interp_data;
    std::array<std::vector<double>, CH_SPACEDIM> m_interp_coords;
    // this is the really long type used for integrands
    // the vector<double> is a vector of all the extracted variables at that
    // point in the order they were added
    using integrand_t =
        std::function<double(std::vector<double> &, double, double, double)>;
    std::vector<integrand_t> m_integrands;
    std::vector<std::array<IntegrationMethod, 2>> m_integration_methods;
    std::vector<std::reference_wrapper<std::vector<double>>> m_integrals;
    std::vector<bool> m_broadcast_integrals;

    bool m_done_extraction; //!< whether or not the extract function has been
                            //!< called for this object

    //! returns the flattened index for m_interp_data and m_interp_coords
    //! associated to given surface, u and v indices
#if CH_SPACEDIM == 3
    int index(int a_isurface, int a_iu, int a_iv) const
    {
        return a_isurface * m_params.num_points_u * m_params.num_points_v +
               a_iu * m_params.num_points_v + a_iv;
    }
#elif CH_SPACEDIM == 2
    int index(int a_isurface, int a_iu) const
    {
        return a_isurface * m_params.num_points_u + a_iu;
    }
#endif

  public:
    //! Normal constructor which requires vars to be added after construction
    //! using add_var or add_vars
    SurfaceExtraction(const SurfaceGeometry &a_geom,
                      const surface_extraction_params_t &a_params, double a_dt,
                      double a_time, bool a_first_step,
                      double a_restart_time = 0.0);

    //! add a single variable or derivative of variable
    void add_var(int a_var,
                 const VariableType var_type = VariableType::evolution,
                 const Derivative &a_deriv = Derivative::LOCAL);

    //! add a vector of variables/derivatives of variables
    void add_vars(const std::vector<vars_t> &a_vars);

    //! add a vector of evolution variables (no derivatives)
    void add_evolution_vars(const std::vector<int> &a_vars);

    //! add a vector of diagnostic variables (no derivatives)
    void add_diagnostic_vars(const std::vector<int> &a_vars);

    //! Alternative constructor with a predefined vector of variables and
    //! derivatives
    SurfaceExtraction(
        const SurfaceGeometry &a_geom,
        const surface_extraction_params_t &a_params,
        const std::vector<std::tuple<int, VariableType, Derivative>> &a_vars,
        double a_dt, double a_time, bool a_first_step,
        double a_restart_time = 0.0);

    //! Another alternative constructor with a predefined vector of variables
    //! no derivatives
    SurfaceExtraction(const SurfaceGeometry &a_geom,
                      const surface_extraction_params_t &a_params,
                      const std::vector<int> &a_vars, double a_dt,
                      double a_time, bool a_first_step,
                      double a_restart_time = 0.0);

    //! Do the extraction
    template <typename InterpAlgo>
    void extract(AMRInterpolator<InterpAlgo> *a_interpolator);

    //! Add an integrand dependent on the interpolated data over the surface
    //! for integrate() to integrate over.
    //! Note the area_element is already included from the SurfaceGeometry
    //! template class
    //! The last argument is whether to broadcast the result to all MPI ranks
    //! or just keep on rank 0. Most use cases won't need this set to true.
    void add_integrand(
        const integrand_t &a_integrand, std::vector<double> &out_integrals,
        const IntegrationMethod &a_method_u = IntegrationMethod::trapezium,
        const IntegrationMethod &a_method_v = IntegrationMethod::trapezium,
        const bool a_broadcast_integral = false);

    //! Add an integrand which is just a single var. The a_var argument should
    //! correspond to the order in which the desired var was added to this
    //! object with add_var
    //! The last argument is whether to broadcast the result to all MPI ranks
    //! or just keep on rank 0. Most use cases won't need this set to true.
    void add_var_integrand(
        int a_var, std::vector<double> &out_integrals,
        const IntegrationMethod &a_method_u = IntegrationMethod::trapezium,
        const IntegrationMethod &a_method_v = IntegrationMethod::trapezium,
        const bool a_broadcast_integral = false);

    //! Integrate the integrands added using add_integrand
    void integrate();

    //! This integrate function can be used if you only want to integrate one
    //! integrand. It calls add_integrand() and integrate()
    //! The last argument is whether to broadcast the result to all MPI ranks
    //! or just keep on rank 0. Most use cases won't need this set to true.
    std::vector<double> integrate(
        integrand_t a_integrand,
        const IntegrationMethod &a_method_u = IntegrationMethod::trapezium,
        const IntegrationMethod &a_method_v = IntegrationMethod::trapezium,
        const bool a_broadcast_integral = false);

    //! Write the interpolated data to a file with a block for each surface
    void write_extraction(std::string a_file_prefix) const;

    //! write some integrals to a file at this timestep
    void write_integrals(const std::string &a_filename,
                         const std::vector<std::vector<double>> &a_integrals,
                         const std::vector<std::string> &a_labels = {}) const;

    //! convenience caller for write_integrals in the case of just integral per
    //! surface
    void write_integral(const std::string &a_filename,
                        const std::vector<double> a_integrals,
                        const std::string a_label = "") const;
};

#include "SurfaceExtraction.impl.hpp"

#endif /* SURFACEEXTRACTION_HPP_ */

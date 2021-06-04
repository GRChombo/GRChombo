/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SURFACEEXTRACTION_HPP_)
#error "This file should only be included through SurfaceExtraction.hpp"
#endif

#ifndef SURFACEEXTRACTION_IMPL_HPP_
#define SURFACEEXTRACTION_IMPL_HPP_

//! Normal constructor which requires vars to be added after construction
//! using add_var or add_vars
template <class SurfaceGeometry>
SurfaceExtraction<SurfaceGeometry>::SurfaceExtraction(
    const SurfaceGeometry &a_geom, const params_t &a_params, double a_dt,
    double a_time, bool a_first_step, double a_restart_time)
    : m_geom(a_geom), m_params(a_params), m_dt(a_dt), m_time(a_time),
      m_first_step(a_first_step), m_restart_time(a_restart_time),
      m_num_interp_points((procID() == 0)
                              ? m_params.num_surfaces * m_params.num_points_u *
                                    m_params.num_points_v
                              : 0),
      m_du(m_geom.du(m_params.num_points_u)),
      m_dv(m_geom.dv(m_params.num_points_v)), m_done_extraction(false)
{
    // check folders only in first two timesteps
    // (or at m_first_step if this is not the first two timesteps)
    if (m_time < m_restart_time + 1.5 * m_dt || m_first_step)
    {
        if (!FilesystemTools::directory_exists(m_params.data_path))
            FilesystemTools::mkdir_recursive(m_params.data_path);

        if (m_params.write_extraction &&
            !FilesystemTools::directory_exists(m_params.extraction_path))
            FilesystemTools::mkdir_recursive(m_params.extraction_path);
    }

    // only interp points on rank 0
    if (procID() == 0)
    {
        FOR1(idir) { m_interp_coords[idir].resize(m_num_interp_points); }

        for (int isurface = 0; isurface < m_params.num_surfaces; ++isurface)
        {
            double surface_param_value =
                m_params.surface_param_values[isurface];
            for (int iu = 0; iu < m_params.num_points_u; ++iu)
            {
                double u = m_geom.u(iu, m_params.num_points_u);
                for (int iv = 0; iv < m_params.num_points_v; ++iv)
                {
                    double v = m_geom.v(iv, m_params.num_points_v);
                    FOR1(idir)
                    {
                        int idx = index(isurface, iu, iv);
                        m_interp_coords[idir][idx] = m_geom.get_grid_coord(
                            idir, surface_param_value, u, v);
                    }
                }
            }
        }
    }

    m_params.validate_extrapolation_radii();
}

//! add a single variable or derivative of variable
template <class SurfaceGeometry>
void SurfaceExtraction<SurfaceGeometry>::add_var(int a_var,
                                                 const VariableType a_var_type,
                                                 const Derivative &a_deriv)
{
    CH_assert(!m_done_extraction);
    m_vars.push_back(std::make_tuple(a_var, a_var_type, a_deriv));
    // m_num_interp_points is 0 on ranks > 0
    m_interp_data.emplace_back(m_num_interp_points);
}

//! add a vector of variables/derivatives of variables
template <class SurfaceGeometry>
void SurfaceExtraction<SurfaceGeometry>::add_vars(
    const std::vector<std::tuple<int, VariableType, Derivative>> &a_vars)
{
    for (auto var : a_vars)
    {
        add_var(std::get<0>(var), std::get<1>(var), std::get<2>(var));
    }
}

//! add a vector of evolutionvariables (no derivatives)
template <class SurfaceGeometry>
void SurfaceExtraction<SurfaceGeometry>::add_evolution_vars(
    const std::vector<int> &a_vars)
{
    for (auto var : a_vars)
    {
        add_var(var, VariableType::evolution);
    }
}

//! add a vector of evolutionvariables (no derivatives)
template <class SurfaceGeometry>
void SurfaceExtraction<SurfaceGeometry>::add_diagnostic_vars(
    const std::vector<int> &a_vars)
{
    for (auto var : a_vars)
    {
        add_var(var, VariableType::diagnostic);
    }
}

//! Alternative constructor with a predefined vector of variables and
//! derivatives
template <class SurfaceGeometry>
SurfaceExtraction<SurfaceGeometry>::SurfaceExtraction(
    const SurfaceGeometry &a_geom, const params_t &a_params,
    const std::vector<std::tuple<int, VariableType, Derivative>> &a_vars,
    double a_dt, double a_time, bool a_first_step, double a_restart_time)
    : SurfaceExtraction<SurfaceGeometry>(a_geom, a_params, a_dt, a_time,
                                         a_first_step, a_restart_time)
{
    add_vars(a_vars);
}

//! Another alternative constructor with a predefined vector of variables
//! no derivatives
template <class SurfaceGeometry>
SurfaceExtraction<SurfaceGeometry>::SurfaceExtraction(
    const SurfaceGeometry &a_geom, const params_t &a_params,
    const std::vector<int> &a_vars, double a_dt, double a_time,
    bool a_first_step, double a_restart_time)
    : SurfaceExtraction<SurfaceExtraction>(a_geom, a_params, a_dt, a_time,
                                           a_first_step, a_restart_time)
{
    add_evolution_vars(a_vars);
}

//! Do the extraction
template <class SurfaceGeometry>
template <typename InterpAlgo>
void SurfaceExtraction<SurfaceGeometry>::extract(
    AMRInterpolator<InterpAlgo> *a_interpolator)
{
    if (a_interpolator == nullptr)
    {
        MayDay::Error("SurfaceExtraction: invalid AMRInterpolator pointer");
    }
    // m_num_interp_points is 0 on ranks > 0
    InterpolationQuery query(m_num_interp_points);
    FOR1(idir) { query.setCoords(idir, m_interp_coords[idir].data()); }
    for (int ivar = 0; ivar < m_vars.size(); ++ivar)
    {
        // note the difference in order between the m_vars tuple in this class
        // and the InterpolationQuery::out_t type
        query.addComp(std::get<0>(m_vars[ivar]), m_interp_data[ivar].data(),
                      std::get<2>(m_vars[ivar]), std::get<1>(m_vars[ivar]));
    }

    // submit the query
    a_interpolator->interp(query);
    m_done_extraction = true;
}

//! Add an integrand (which must of type integrand_t) for integrate() to
//! integrate over. Note the area_element is already included from the
//! SurfaceGeometry template class
template <class SurfaceGeometry>
void SurfaceExtraction<SurfaceGeometry>::add_integrand(
    const integrand_t &a_integrand, std::vector<double> &out_integrals,
    const IntegrationMethod &a_method_u, const IntegrationMethod &a_method_v,
    const bool a_broadcast_integral)
{
    // if broadcasting integral to all ranks, all ranks need to know about it
    m_broadcast_integrals.push_back(a_broadcast_integral);
    if (a_broadcast_integral || procID() == 0)
    {
        // resize the out_integrals and store a reference to it
        out_integrals.resize(m_params.num_surfaces);
        std::fill(out_integrals.begin(), out_integrals.end(), 0.0);
    }
    m_integrals.push_back(std::ref(out_integrals));

    // only rank 0 actually does the integration and needs to know about the
    // integrand and integration method
    if (procID() == 0)
    {
        // store the integrand
        m_integrands.push_back(a_integrand);

        // check if integration methods are valid given periodicity and number
        // of points
        bool valid_u =
            a_method_u.is_valid(m_params.num_points_u, m_geom.is_u_periodic());
        bool valid_v =
            a_method_v.is_valid(m_params.num_points_v, m_geom.is_v_periodic());

        // default to using the trapezium rule if provided methods are not valid
        IntegrationMethod method_u = IntegrationMethod::trapezium;
        IntegrationMethod method_v = IntegrationMethod::trapezium;
        if (!valid_u)
        {
            MayDay::Warning(
                "SurfaceExtraction<SurfaceGeometry>::integrate: Provided "
                "IntegrationMethod for u is not valid with\nthis num_points_u; "
                "reverting to trapezium rule.");
        }
        else
        {
            method_u = a_method_u;
        }
        if (!valid_v)
        {
            MayDay::Warning(
                "SurfaceExtraction<SurfaceGeometry>::integrate: Provided "
                "IntegrationMethod for v is not valid with\nthis num_points_v; "
                "reverting to trapezium rule.");
        }
        else
        {
            method_v = a_method_v;
        }
        m_integration_methods.push_back({method_u, method_v});
    }
}

//! Add an integrand which is just a single var. The a_var argument should
//! correspond to the order in which the desired var was added to this object
//! with add_var
template <class SurfaceGeometry>
void SurfaceExtraction<SurfaceGeometry>::add_var_integrand(
    int a_var, std::vector<double> &out_integrals,
    const IntegrationMethod &a_method_u, const IntegrationMethod &a_method_v,
    const bool a_broadcast_integral)
{
    CH_assert(a_var >= 0 && a_var < m_vars.size());
    integrand_t var_integrand = [var = a_var](std::vector<double> &data, double,
                                              double,
                                              double) { return data[var]; };
    add_integrand(var_integrand, out_integrals, a_method_u, a_method_v,
                  a_broadcast_integral);
}

template <class SurfaceGeometry>
void SurfaceExtraction<SurfaceGeometry>::integrate()
{
    CH_assert(m_done_extraction);
    if (procID() == 0)
    {
        // note this condition won't be true on other ranks
        CH_assert(m_integrands.size() == m_integration_methods.size() &&
                  m_integrals.size() > 0);
        int num_integrals = m_integrals.size();

        for (int isurface = 0; isurface < m_params.num_surfaces; ++isurface)
        {
            double surface_param_value =
                m_params.surface_param_values[isurface];
            for (int iu = 0; iu < m_params.num_points_u; ++iu)
            {
                double u = m_geom.u(iu, m_params.num_points_u);
                std::vector<double> inner_integral(num_integrals, 0.0);
                for (int iv = 0; iv < m_params.num_points_v; ++iv)
                {
                    double v = m_geom.v(iv, m_params.num_points_v);
                    std::vector<double> data_here(m_vars.size());
                    for (int ivar = 0; ivar < m_vars.size(); ++ivar)
                    {
                        data_here[ivar] =
                            m_interp_data[ivar][index(isurface, iu, iv)];
                    }
                    for (int iintegral = 0; iintegral < num_integrals;
                         ++iintegral)
                    {
                        auto integrand = m_integrands[iintegral];
                        double integrand_with_area_element =
                            integrand(data_here, surface_param_value, u, v) *
                            m_geom.area_element(surface_param_value, u, v);
                        double weight =
                            m_integration_methods[iintegral][1].weight(
                                iv, m_params.num_points_v,
                                m_geom.is_v_periodic());
                        inner_integral[iintegral] +=
                            weight * m_dv * integrand_with_area_element;
                    }
                }
                for (int iintegral = 0; iintegral < num_integrals; ++iintegral)
                {
                    double weight = m_integration_methods[iintegral][0].weight(
                        iu, m_params.num_points_u, m_geom.is_u_periodic());
                    (m_integrals[iintegral].get())[isurface] +=
                        weight * m_du * inner_integral[iintegral];
                }
            }
        }
    }

    // now broadcast result to non-zero ranks if requested
    for (int iintegral = 0; iintegral < m_integrals.size(); ++iintegral)
    {
        if (m_broadcast_integrals[iintegral])
        {
            Vector<double> broadcast_Vector;
            if (procID() == 0)
                broadcast_Vector = m_integrals[iintegral].get();
            broadcast(broadcast_Vector, 0);
            if (procID() != 0)
            {
                m_integrals[iintegral].get() = broadcast_Vector.stdVector();
            }
        }
    }
}

//! Integrate some integrand dependent on the interpolated data over the
//! surface. The integrand function should be of the signature
//! double integrand(std::vector<double> data_here,
//!     double a_surface_param_value, double a_u, double a_v)
//! where data_here is a vector of all the interpolated variables at the
//! point specified by the other arguments.
template <class SurfaceGeometry>
std::vector<double> SurfaceExtraction<SurfaceGeometry>::integrate(
    integrand_t a_integrand, const IntegrationMethod &a_method_u,
    const IntegrationMethod &a_method_v, const bool a_broadcast_integral)
{
    m_integrands.clear();
    m_integration_methods.clear();
    m_integrals.clear();

    std::vector<double> out_integrals(m_params.num_surfaces, 0.0);
    add_integrand(a_integrand, out_integrals, a_method_u, a_method_v,
                  a_broadcast_integral);
    integrate();

    return out_integrals;
}

//! Write the interpolated data to a file with a block for each surface
template <class SurfaceGeometry>
void SurfaceExtraction<SurfaceGeometry>::write_extraction(
    std::string a_file_prefix) const
{
    CH_assert(m_done_extraction);
    if (procID() == 0)
    {
        SmallDataIO extraction_file(m_params.extraction_path + a_file_prefix,
                                    m_dt, m_time, m_restart_time,
                                    SmallDataIO::NEW, m_first_step);

        for (int isurface = 0; isurface < m_params.num_surfaces; ++isurface)
        {
            // Write headers
            std::vector<std::string> header1_strings = {
                "time = " + std::to_string(m_time) + ",",
                m_geom.param_name() + " = " +
                    std::to_string(m_params.surface_param_values[isurface])};
            extraction_file.write_header_line(header1_strings, "");
            std::vector<std::string> components(m_vars.size());
            for (int ivar = 0; ivar < m_vars.size(); ++ivar)
            {
                if (std::get<2>(m_vars[ivar]) != Derivative::LOCAL)
                {
                    components[ivar] =
                        Derivative::name(std::get<2>(m_vars[ivar])) + "_";
                }
                else
                {
                    components[ivar] = "";
                }
                if (std::get<1>(m_vars[ivar]) == VariableType::evolution)
                {
                    components[ivar] +=
                        UserVariables::variable_names[std::get<0>(
                            m_vars[ivar])];
                }
                else
                {
                    components[ivar] +=
                        DiagnosticVariables::variable_names[std::get<0>(
                            m_vars[ivar])];
                }
            }
            std::vector<std::string> coords = {m_geom.u_name(),
                                               m_geom.v_name()};
            extraction_file.write_header_line(components, coords);

            // Now the data
            for (int iu = 0; iu < m_params.num_points_u; ++iu)
            {
                double u = m_geom.u(iu, m_params.num_points_u);
                for (int iv = 0; iv < m_params.num_points_v; ++iv)
                {
                    double v = m_geom.v(iv, m_params.num_points_v);
                    int idx = index(isurface, iu, iv);
                    std::vector<double> data(m_vars.size());
                    for (int ivar = 0; ivar < m_vars.size(); ++ivar)
                    {
                        data[ivar] = m_interp_data[ivar][idx];
                    }

                    extraction_file.write_data_line(data, {u, v});
                }
            }
            extraction_file.line_break();
        }
    }
}

//! write some integrals to a file at this timestep
template <class SurfaceGeometry>
void SurfaceExtraction<SurfaceGeometry>::write_integrals(
    const std::string &a_filename,
    const std::vector<std::vector<double>> &a_integrals,
    const std::vector<std::string> &a_labels) const
{
    if (procID() == 0)
    {
        std::vector<double> extrapolations =
            richardson_extrapolation(a_integrals);
        bool extrapolation_done = extrapolations.size() > 0;

        const int num_integrals_per_surface = a_integrals.size();
        const int total_num_surfaces =
            m_params.num_surfaces + (extrapolation_done > 0 ? 1 : 0);

        // if labels are provided there must be the same number of labels as
        // there are integrals
        if (!a_labels.empty())
        {
            CH_assert(num_integrals_per_surface == a_labels.size());
        }
        // each inner vector element of a_integrals must have the same number of
        // elements as there are surfaces (i.e. one integral per surface)
        for (auto vect : a_integrals)
        {
            CH_assert(vect.size() == m_params.num_surfaces);
        }
        // open file for writing
        SmallDataIO integral_file(m_params.data_path + a_filename, m_dt, m_time,
                                  m_restart_time, SmallDataIO::APPEND,
                                  m_first_step);

        // remove any duplicate data if this is a restart
        integral_file.remove_duplicate_time_data();

        if (m_first_step)
        {
            // make header strings
            std::vector<std::string> header1_strings(num_integrals_per_surface *
                                                     total_num_surfaces);
            std::vector<std::string> header2_strings(num_integrals_per_surface *
                                                     total_num_surfaces);
            for (int iintegral = 0; iintegral < num_integrals_per_surface;
                 ++iintegral)
            {
                for (int isurface = 0; isurface < m_params.num_surfaces;
                     ++isurface)
                {
                    int idx = isurface * num_integrals_per_surface + iintegral;
                    if (a_labels.empty())
                        header1_strings[idx] = "";
                    else
                        header1_strings[idx] = a_labels[iintegral];
                    header2_strings[idx] =
                        std::to_string(m_params.surface_param_values[isurface]);
                }
                if (extrapolation_done)
                {
                    int idx =
                        m_params.num_surfaces * num_integrals_per_surface +
                        iintegral;
                    header1_strings[idx] = a_labels[iintegral];
                    header2_strings[idx] = "Infinity";
                }
            }
            std::string pre_header2_string = m_geom.param_name() + " = ";

            // write headers
            integral_file.write_header_line(header1_strings);
            integral_file.write_header_line(header2_strings,
                                            pre_header2_string);
        }

        // make vector of data for writing
        std::vector<double> data_for_writing(num_integrals_per_surface *
                                             total_num_surfaces);
        for (int iintegral = 0; iintegral < num_integrals_per_surface;
             ++iintegral)
        {
            for (int isurface = 0; isurface < m_params.num_surfaces; ++isurface)
            {
                int idx = isurface * num_integrals_per_surface + iintegral;
                data_for_writing[idx] = a_integrals[iintegral][isurface];
            }
            if (extrapolation_done)
            {
                int idx = m_params.num_surfaces * num_integrals_per_surface +
                          iintegral;
                data_for_writing[idx] = extrapolations[iintegral];
            }
        }

        // write data
        integral_file.write_time_data_line(data_for_writing);
    }
}

//! convenience caller for write_integrals in the case of just one integral per
//! surface
template <class SurfaceGeometry>
void SurfaceExtraction<SurfaceGeometry>::write_integral(
    const std::string &a_filename, const std::vector<double> &a_integrals,
    const std::string &a_label) const
{
    write_integrals(a_filename, {a_integrals}, {a_label});
}

template <class SurfaceGeometry>
void SurfaceExtraction<
    SurfaceGeometry>::params_t::validate_extrapolation_radii()
{
    std::vector<int> valid_radii;
    for (int i = 0; i < radii_idxs_for_extrapolation.size(); ++i)
    {
        // if valid
        if (radii_idxs_for_extrapolation[i] < num_surfaces &&
            radii_idxs_for_extrapolation[i] > -num_surfaces)
        {
            // allow negative indices, such that '-1' is 'last', '-2' is 'one
            // before last'
            if (radii_idxs_for_extrapolation[i] < 0)
                radii_idxs_for_extrapolation[i] += num_surfaces;

            // if not repeated already
            if (i == 0 || std::find(radii_idxs_for_extrapolation.begin(),
                                    radii_idxs_for_extrapolation.begin() + i,
                                    radii_idxs_for_extrapolation[i]) ==
                              radii_idxs_for_extrapolation.begin() + i)
            {
                valid_radii.push_back(radii_idxs_for_extrapolation[i]);
            }
        }
    }
    radii_idxs_for_extrapolation = valid_radii;
}

// function to apply richardson extrapolation of 2nd or 3rd order to calculated
// integrals
template <class SurfaceGeometry>
std::vector<double>
SurfaceExtraction<SurfaceGeometry>::richardson_extrapolation(
    const std::vector<std::vector<double>> &integrals) const
{
    CH_TIME("SurfaceExtraction::richardson_extrapolation");

    // already validated the radii
    int extrapolation_order = m_params.radii_idxs_for_extrapolation.size();
    int num_comps = integrals.size();

    if (extrapolation_order < 2 || num_comps == 0 ||
        integrals[0].size() < extrapolation_order)
    {
        return std::vector<double>();
    }

    std::vector<double> extrapolations(num_comps);

    if (extrapolation_order == 3)
    {
        int i1 = m_params.radii_idxs_for_extrapolation[0];
        int i2 = m_params.radii_idxs_for_extrapolation[1];
        int i3 = m_params.radii_idxs_for_extrapolation[2];
        double r1 = m_params.surface_param_values[i1];
        double r2 = m_params.surface_param_values[i2];
        double r3 = m_params.surface_param_values[i3];
        double c2 = r2 / r1 * (1. - r3 / r1) / (1. - r3 / r2);
        double c3 = -r3 / r1 * (1. - r2 / r1) / (r2 / r3 - 1.);
        for (int icomp = 0; icomp < num_comps; ++icomp)
        {
            double comp_inf =
                (integrals[icomp][i1] - c2 * integrals[icomp][i2] -
                 c3 * integrals[icomp][i3]) /
                (1. - c2 - c3);
            extrapolations[icomp] = comp_inf;
        }
    }
    else if (extrapolation_order == 2)
    {
        int i1 = m_params.radii_idxs_for_extrapolation[0];
        int i2 = m_params.radii_idxs_for_extrapolation[1];
        double r1 = m_params.surface_param_values[i1];
        double r2 = m_params.surface_param_values[i2];
        for (int icomp = 0; icomp < num_comps; ++icomp)
        {
            double comp_inf =
                (integrals[icomp][i1] - r2 / r1 * integrals[icomp][i2]) /
                (1. - r2 / r1);
            extrapolations[icomp] = comp_inf;
        }
    }

    return extrapolations;
}

#endif /* SURFACEEXTRACTION_IMPL_HPP_ */

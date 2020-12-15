/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "BoundaryConditions.hpp"
#include "FArrayBox.H"
#include "ProblemDomain.H"
#include "RealVect.H"
#include <algorithm>
#include <array>
#include <map>
#include <string>

BoundaryConditions::params_t::params_t()
{
    // set defaults
    hi_boundary.fill(STATIC_BC);
    lo_boundary.fill(STATIC_BC);
    is_periodic.fill(true);
    nonperiodic_boundaries_exist = false;
    boundary_solution_enforced = false;
    boundary_rhs_enforced = false;
    reflective_boundaries_exist = false;
    sommerfeld_boundaries_exist = false;
    vars_parity.fill(BoundaryConditions::UNDEFINED);
    vars_parity_diagnostic.fill(BoundaryConditions::EXTRAPOLATING_BC);
    vars_asymptotic_values.fill(0.0);
    extrapolation_order = 1;
    mixed_bc_extrapolating_vars.clear();
    mixed_bc_sommerfeld_vars.clear();
}

void BoundaryConditions::params_t::set_is_periodic(
    const std::array<bool, CH_SPACEDIM> &a_is_periodic)
{
    is_periodic = a_is_periodic;
    FOR1(idir)
    {
        if (!is_periodic[idir])
            nonperiodic_boundaries_exist = true;
    }
}
void BoundaryConditions::params_t::set_hi_boundary(
    const std::array<int, CH_SPACEDIM> &a_hi_boundary)
{
    FOR1(idir)
    {
        if (!is_periodic[idir])
        {
            hi_boundary[idir] = a_hi_boundary[idir];
            if (hi_boundary[idir] == REFLECTIVE_BC)
            {
                boundary_rhs_enforced = true;
                boundary_solution_enforced = true;
                reflective_boundaries_exist = true;
            }
            else if (hi_boundary[idir] == SOMMERFELD_BC)
            {
                boundary_rhs_enforced = true;
                sommerfeld_boundaries_exist = true;
            }
            else if (hi_boundary[idir] == EXTRAPOLATING_BC)
            {
                boundary_rhs_enforced = true;
                boundary_solution_enforced = true;
                extrapolating_boundaries_exist = true;
            }
            else if (hi_boundary[idir] == MIXED_BC)
            {
                boundary_rhs_enforced = true;
                boundary_solution_enforced = true;
                sommerfeld_boundaries_exist = true;
                extrapolating_boundaries_exist = true;
                mixed_boundaries_exist = true;
            }
        }
    }
}
void BoundaryConditions::params_t::set_lo_boundary(
    const std::array<int, CH_SPACEDIM> &a_lo_boundary)
{
    FOR1(idir)
    {
        if (!is_periodic[idir])
        {
            lo_boundary[idir] = a_lo_boundary[idir];
            if (lo_boundary[idir] == REFLECTIVE_BC)
            {
                boundary_solution_enforced = true;
                reflective_boundaries_exist = true;
            }
            else if (lo_boundary[idir] == SOMMERFELD_BC)
            {
                boundary_rhs_enforced = true;
                sommerfeld_boundaries_exist = true;
            }
            else if (lo_boundary[idir] == EXTRAPOLATING_BC)
            {
                boundary_rhs_enforced = true;
                boundary_solution_enforced = true;
                extrapolating_boundaries_exist = true;
            }
            else if (lo_boundary[idir] == MIXED_BC)
            {
                boundary_rhs_enforced = true;
                boundary_solution_enforced = true;
                sommerfeld_boundaries_exist = true;
                extrapolating_boundaries_exist = true;
                mixed_boundaries_exist = true;
            }
        }
    }
}

void BoundaryConditions::params_t::read_params(GRParmParse &pp)
{
    using namespace UserVariables; // for loading the arrays

    // still load even if not contained, to ensure printout saying parameters
    // were set to their default values
    std::array<bool, CH_SPACEDIM> isPeriodic;
    pp.load("isPeriodic", isPeriodic, is_periodic);
    if (pp.contains("isPeriodic"))
        set_is_periodic(isPeriodic);

    std::array<int, CH_SPACEDIM> hiBoundary;
    pp.load("hi_boundary", hiBoundary, hi_boundary);
    if (pp.contains("hi_boundary"))
        set_hi_boundary(hiBoundary);

    std::array<int, CH_SPACEDIM> loBoundary;
    pp.load("lo_boundary", loBoundary, lo_boundary);
    if (pp.contains("lo_boundary"))
        set_lo_boundary(loBoundary);

    if (reflective_boundaries_exist)
    {
        pp.load("vars_parity", vars_parity);
        if (pp.contains("vars_parity_diagnostic"))
            pp.load("vars_parity_diagnostic", vars_parity_diagnostic);
    }
    if (sommerfeld_boundaries_exist)
    {
        if (pp.contains("num_nonzero_asymptotic_vars"))
        {
            int num_values = 0;
            std::vector<std::pair<int, VariableType>> nonzero_asymptotic_vars;
            load_vars_to_vector(pp, "nonzero_asymptotic_vars",
                                "num_nonzero_asymptotic_vars",
                                nonzero_asymptotic_vars, num_values);
            const double default_value = 0.0;
            load_values_to_array(pp, "nonzero_asymptotic_values",
                                 nonzero_asymptotic_vars,
                                 vars_asymptotic_values, default_value);
        }
        // for backwards compatibility, but above method should
        // be preferred in future, as is less error prone
        else
        {
            pp.load("vars_asymptotic_values", vars_asymptotic_values);
        }
    }
    if (extrapolating_boundaries_exist)
    {
        pp.load("extrapolation_order", extrapolation_order, 1);
    }
    if (mixed_boundaries_exist)
    {
        int num_extrapolating_vars = 0;
        std::vector<std::pair<int, VariableType>> extrapolating_vars;
        load_vars_to_vector(pp, "extrapolating_vars", "num_extrapolating_vars",
                            extrapolating_vars, num_extrapolating_vars);
        for (int icomp = 0; icomp < NUM_VARS; icomp++)
        {
            bool is_extrapolating = false;
            // if the variable is not in extrapolating vars, it
            // is assumed to be sommerfeld by default
            for (int icomp2 = 0; icomp2 < extrapolating_vars.size(); icomp2++)
            {
                if (icomp == extrapolating_vars[icomp2].first)
                {
                    // should be an evolution variable
                    CH_assert(extrapolating_vars[icomp2].second ==
                              VariableType::evolution);
                    mixed_bc_extrapolating_vars.push_back(icomp);
                    is_extrapolating = true;
                }
            }
            if (!is_extrapolating)
            {
                mixed_bc_sommerfeld_vars.push_back(icomp);
            }
        }
    }
    if (nonperiodic_boundaries_exist)
    {
        // write out boundary conditions where non periodic - useful for
        // debug
        write_boundary_conditions(*this);
    }
}

/// define function sets members and is_defined set to true
void BoundaryConditions::define(double a_dx,
                                std::array<double, CH_SPACEDIM> a_center,
                                const params_t &a_params,
                                ProblemDomain a_domain, int a_num_ghosts)
{
    m_dx = a_dx;
    m_params = a_params;
    m_domain = a_domain;
    m_domain_box = a_domain.domainBox();
    m_num_ghosts = a_num_ghosts;
    FOR1(i) { m_center[i] = a_center[i]; }
    // make a vector of the vars indices, [0,1,2...NUM_VARS]
    m_evolution_comps.resize(NUM_VARS);
    std::iota(m_evolution_comps.begin(), m_evolution_comps.end(), 0);
    m_diagnostic_comps.resize(NUM_DIAGNOSTIC_VARS);
    std::iota(m_diagnostic_comps.begin(), m_diagnostic_comps.end(), 0);
    is_defined = true;
}

/// change the asymptotic values of the variables for the Sommerfeld BCs
/// this will allow them to evolve during a simulation if necessary
void BoundaryConditions::set_vars_asymptotic_values(
    std::array<double, NUM_VARS> &vars_asymptotic_values)
{
    m_params.vars_asymptotic_values = vars_asymptotic_values;
}

void BoundaryConditions::write_reflective_conditions(int idir,
                                                     const params_t &a_params)
{
    pout() << "The variables that are parity odd in this direction are : "
           << endl;
    for (int icomp = 0; icomp < NUM_VARS; icomp++)
    {
        int parity = get_var_parity(icomp, idir, a_params);
        if (parity == -1)
        {
            pout() << UserVariables::variable_names[icomp] << "    ";
        }
    }
    for (int icomp = 0; icomp < NUM_DIAGNOSTIC_VARS; icomp++)
    {
        int parity =
            get_var_parity(icomp, idir, a_params, VariableType::diagnostic);
        if (parity == -1)
        {
            pout() << DiagnosticVariables::variable_names[icomp] << "    ";
        }
    }
}

void BoundaryConditions::write_sommerfeld_conditions(int idir,
                                                     const params_t &a_params)
{
    pout() << "The non zero asymptotic values of the variables "
              "in this direction are : "
           << endl;
    for (int icomp = 0; icomp < NUM_VARS; icomp++)
    {
        if (a_params.vars_asymptotic_values[icomp] != 0)
        {
            pout() << UserVariables::variable_names[icomp] << " = "
                   << a_params.vars_asymptotic_values[icomp] << "    ";
        }
    }
    // not done for diagnostics
}

void BoundaryConditions::write_mixed_conditions(int idir,
                                                const params_t &a_params)
{
    // check all the vars have been assigned a BC - this should always be the
    // case because of how the params are assigned
    CH_assert(a_params.mixed_bc_sommerfeld_vars.size() +
                  a_params.mixed_bc_extrapolating_vars.size() ==
              NUM_VARS);

    // now do the write out
    pout()
        << "The variables that use extrapolating bcs in this direction are : "
        << endl;
    for (int icomp = 0; icomp < NUM_VARS; icomp++)
    {
        std::vector<int> v = a_params.mixed_bc_extrapolating_vars;
        if (std::binary_search(v.begin(), v.end(), icomp))
        {
            pout() << UserVariables::variable_names[icomp] << "    ";
        }
    }
    pout() << endl;
    pout() << "The other variables all use Sommerfeld boundary conditions."
           << endl;
    write_sommerfeld_conditions(idir, a_params);
}

/// write out boundary params (used during setup for debugging)
void BoundaryConditions::write_boundary_conditions(const params_t &a_params)
{
    pout() << "You are using non periodic boundary conditions." << endl;
    pout() << "The boundary params chosen are:  " << endl;
    pout() << "---------------------------------" << endl;

    std::map<int, std::string> bc_names = {{STATIC_BC, "Static"},
                                           {SOMMERFELD_BC, "Sommerfeld"},
                                           {REFLECTIVE_BC, "Reflective"},
                                           {EXTRAPOLATING_BC, "Extrapolating"},
                                           {MIXED_BC, "Mixed"}};
    FOR1(idir)
    {
        if (!a_params.is_periodic[idir])
        {
            pout() << "- " << bc_names[a_params.hi_boundary[idir]]
                   << " boundaries in direction high " << idir << endl;
            // high directions
            if (a_params.hi_boundary[idir] == REFLECTIVE_BC)
            {
                write_reflective_conditions(idir, a_params);
            }
            else if (a_params.hi_boundary[idir] == SOMMERFELD_BC)
            {
                write_sommerfeld_conditions(idir, a_params);
            }
            else if (a_params.hi_boundary[idir] == MIXED_BC)
            {
                write_mixed_conditions(idir, a_params);
            }
            pout() << endl;

            // low directions
            pout() << "- " << bc_names[a_params.lo_boundary[idir]]
                   << " boundaries in direction low " << idir << endl;
            if (a_params.lo_boundary[idir] == REFLECTIVE_BC)
            {
                write_reflective_conditions(idir, a_params);
            }
            else if (a_params.lo_boundary[idir] == SOMMERFELD_BC)
            {
                write_sommerfeld_conditions(idir, a_params);
            }
            else if (a_params.lo_boundary[idir] == MIXED_BC)
            {
                write_mixed_conditions(idir, a_params);
            }
            pout() << endl;
        }
    }
    pout() << "---------------------------------" << endl;
}

/// The function which returns the parity of each of the vars in
/// UserVariables.hpp The parity should be defined in the params file, and
/// will be output to the pout files for checking at start/restart of
/// simulation (It is only required for reflective boundary conditions.)
int BoundaryConditions::get_var_parity(int a_comp, int a_dir,
                                       const VariableType var_type) const
{
    int var_parity = get_var_parity(a_comp, a_dir, m_params, var_type);

    return var_parity;
}

/// static version used for initial output of boundary values
int BoundaryConditions::get_var_parity(int a_comp, int a_dir,
                                       const params_t &a_params,
                                       const VariableType var_type)
{
    int comp_parity = (var_type == VariableType::evolution
                           ? a_params.vars_parity[a_comp]
                           : a_params.vars_parity_diagnostic[a_comp]);

    int vars_parity = 1;
    if ((a_dir == 0) && (comp_parity == ODD_X || comp_parity == ODD_XY ||
                         comp_parity == ODD_XZ || comp_parity == ODD_XYZ))
    {
        vars_parity = -1;
    }
    else if ((a_dir == 1) && (comp_parity == ODD_Y || comp_parity == ODD_XY ||
                              comp_parity == ODD_YZ || comp_parity == ODD_XYZ))
    {
        vars_parity = -1;
    }
    else if ((a_dir == 2) && (comp_parity == ODD_Z || comp_parity == ODD_XZ ||
                              comp_parity == ODD_YZ || comp_parity == ODD_XYZ))
    {
        vars_parity = -1;
    }
    return vars_parity;
}

/// Fill the rhs boundary values appropriately based on the params set
void BoundaryConditions::fill_rhs_boundaries(const Side::LoHiSide a_side,
                                             const GRLevelData &a_soln,
                                             GRLevelData &a_rhs)
{
    CH_assert(is_defined);
    CH_TIME("BoundaryConditions::fill_rhs_boundaries");

    // cycle through the directions, filling the rhs
    FOR1(idir)
    {
        // only do something if this direction is not periodic
        if (!m_params.is_periodic[idir])
        {
            int boundary_condition = get_boundary_condition(a_side, idir);
            fill_boundary_cells_dir(a_side, a_soln, a_rhs, idir,
                                    boundary_condition);
        }
    }
}

/// fill solution boundary conditions, e.g. after interpolation
void BoundaryConditions::fill_solution_boundaries(const Side::LoHiSide a_side,
                                                  GRLevelData &a_state)
{
    CH_assert(is_defined);
    CH_TIME("BoundaryConditions::fill_solution_boundaries");

    // cycle through the directions
    FOR1(idir)
    {
        // only do something if this direction is not periodic and solution
        // boundary enforced in this direction
        if (!m_params.is_periodic[idir])
        {
            int boundary_condition = get_boundary_condition(a_side, idir);

            // same copying of cells which we require for the rhs solution
            // but tell it we are not filling the rhs for mixed condition
            if ((boundary_condition == REFLECTIVE_BC) ||
                (boundary_condition == EXTRAPOLATING_BC) ||
                (boundary_condition == MIXED_BC))
            {
                const bool filling_rhs = false;
                fill_boundary_cells_dir(a_side, a_state, a_state, idir,
                                        boundary_condition,
                                        VariableType::evolution, filling_rhs);
            }
        }
    }
}

/// fill diagnostic boundaries
void BoundaryConditions::fill_diagnostic_boundaries(const Side::LoHiSide a_side,
                                                    GRLevelData &a_state)
{
    CH_assert(is_defined);
    CH_TIME("BoundaryConditions::fill_diagnostic_boundaries");

    // cycle through the directions
    FOR1(idir)
    {
        // only do something if this direction is not periodic
        if (!m_params.is_periodic[idir])
        {
            int boundary_condition = get_boundary_condition(a_side, idir);
            // for any non reflective BC, we just want to fill the ghosts with
            // something non nan so set the boundary condition to be
            // EXTRAPOLATING
            if (boundary_condition != REFLECTIVE_BC)
            {
                boundary_condition = EXTRAPOLATING_BC;
            }
            const bool filling_rhs = false;
            fill_boundary_cells_dir(a_side, a_state, a_state, idir,
                                    boundary_condition,
                                    VariableType::diagnostic, filling_rhs);
        }
    }
}

/// Fill the boundary values appropriately based on the params set
/// in the direction dir
void BoundaryConditions::fill_boundary_cells_dir(
    const Side::LoHiSide a_side, const GRLevelData &a_soln, GRLevelData &a_out,
    const int dir, const int boundary_condition, const VariableType var_type,
    const bool filling_rhs)
{
    const std::vector<int> &comps =
        (var_type == VariableType::evolution ? m_evolution_comps
                                             : m_diagnostic_comps);

    // iterate through the boxes, shared amongst threads
    DataIterator dit = a_out.dataIterator();
    int nbox = dit.size();
#pragma omp parallel for default(shared)
    for (int ibox = 0; ibox < nbox; ++ibox)
    {
        DataIndex dind = dit[ibox];
        FArrayBox &out_box = a_out[dind];
        const FArrayBox &soln_box = a_soln[dind];
        Box this_box = out_box.box();
        IntVect offset_lo = -this_box.smallEnd() + m_domain_box.smallEnd();
        IntVect offset_hi = +this_box.bigEnd() - m_domain_box.bigEnd();

        // reduce box to the intersection of the box and the
        // problem domain ie remove all outer ghost cells
        this_box &= m_domain_box;
        // get the boundary box (may be Empty)
        Box boundary_box =
            get_boundary_box(a_side, dir, offset_lo, offset_hi, this_box);

        // now we have the appropriate box, fill it!
        BoxIterator bit(boundary_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            switch (boundary_condition)
            {
            // simplest case - boundary values are set to zero
            case STATIC_BC:
            {
                for (int icomp : comps)
                {
                    out_box(iv, icomp) = 0.0;
                }
                break;
            }
            // Sommerfeld is outgoing radiation - only applies to rhs
            case SOMMERFELD_BC:
            {
                fill_sommerfeld_cell(out_box, soln_box, iv, comps);
                break;
            }
            // Enforce a reflective symmetry in some direction
            case REFLECTIVE_BC:
            {
                fill_reflective_cell(out_box, iv, a_side, dir, comps, var_type);
                break;
            }
            case EXTRAPOLATING_BC:
            {
                fill_extrapolating_cell(out_box, iv, a_side, dir, comps,
                                        m_params.extrapolation_order);
                break;
            }
            case MIXED_BC:
            {
                fill_extrapolating_cell(out_box, iv, a_side, dir,
                                        m_params.mixed_bc_extrapolating_vars,
                                        m_params.extrapolation_order);
                if (filling_rhs)
                {
                    fill_sommerfeld_cell(out_box, soln_box, iv,
                                         m_params.mixed_bc_sommerfeld_vars);
                }
                break;
            }
            default:
                MayDay::Error(
                    "BoundaryCondition::Supplied boundary not supported.");
            } // end switch
        }     // end iterate over box
    }         // end iterate over boxes
}

void BoundaryConditions::fill_sommerfeld_cell(
    FArrayBox &rhs_box, const FArrayBox &soln_box, const IntVect iv,
    const std::vector<int> &sommerfeld_comps) const
{
    // assumes an asymptotic value + radial waves and permits them
    // to exit grid with minimal reflections
    // get real position on the grid
    RealVect loc(iv + 0.5 * RealVect::Unit);
    loc *= m_dx;
    loc -= m_center;
    double radius_squared = 0.0;
    FOR1(i) { radius_squared += loc[i] * loc[i]; }
    double radius = sqrt(radius_squared);
    IntVect lo_local_offset = iv - soln_box.smallEnd();
    IntVect hi_local_offset = soln_box.bigEnd() - iv;

    // Apply Sommerfeld BCs to each variable in sommerfeld_comps
    for (int icomp : sommerfeld_comps)
    {
        rhs_box(iv, icomp) = 0.0;
        FOR1(idir2)
        {
            IntVect iv_offset1 = iv;
            IntVect iv_offset2 = iv;
            double d1;
            // bit of work to get the right stencils for near
            // the edges of the domain, only using second order
            // stencils for now
            if (lo_local_offset[idir2] < 1)
            {
                // near lo end
                iv_offset1[idir2] += +1;
                iv_offset2[idir2] += +2;
                d1 = 1.0 / m_dx *
                     (-1.5 * soln_box(iv, icomp) +
                      2.0 * soln_box(iv_offset1, icomp) -
                      0.5 * soln_box(iv_offset2, icomp));
            }
            else if (hi_local_offset[idir2] < 1)
            {
                // near hi end
                iv_offset1[idir2] += -1;
                iv_offset2[idir2] += -2;
                d1 = 1.0 / m_dx *
                     (+1.5 * soln_box(iv, icomp) -
                      2.0 * soln_box(iv_offset1, icomp) +
                      0.5 * soln_box(iv_offset2, icomp));
            }
            else
            {
                // normal case
                iv_offset1[idir2] += +1;
                iv_offset2[idir2] += -1;
                d1 =
                    0.5 / m_dx *
                    (soln_box(iv_offset1, icomp) - soln_box(iv_offset2, icomp));
            }

            // for each direction add dphidx * x^i / r
            rhs_box(iv, icomp) += -d1 * loc[idir2] / radius;
        }

        // asymptotic values - these need to have been set in
        // the params file
        rhs_box(iv, icomp) +=
            (m_params.vars_asymptotic_values[icomp] - soln_box(iv, icomp)) /
            radius;
    }
}

void BoundaryConditions::fill_reflective_cell(
    FArrayBox &out_box, const IntVect iv, const Side::LoHiSide a_side,
    const int dir, const std::vector<int> &reflective_comps,
    const VariableType var_type) const
{
    // assume boundary is a reflection of values within the grid
    // care must be taken with variable parity to maintain correct
    // values on reflection, e.g. x components of vectors are odd
    // parity in the x direction
    IntVect iv_copy = iv;
    /// where to copy the data from - mirror image in domain
    if (a_side == Side::Lo)
    {
        iv_copy[dir] = -iv[dir] - 1;
    }
    else
    {
        iv_copy[dir] = 2 * m_domain_box.bigEnd(dir) - iv[dir] + 1;
    }

    // replace value at iv with value at iv_copy
    for (int icomp : reflective_comps)
    {
        int parity = get_var_parity(icomp, dir, var_type);
        out_box(iv, icomp) = parity * out_box(iv_copy, icomp);
    }
}

void BoundaryConditions::fill_extrapolating_cell(
    FArrayBox &out_box, const IntVect iv, const Side::LoHiSide a_side,
    const int dir, const std::vector<int> &extrapolating_comps,
    const int order) const
{
    for (int icomp : extrapolating_comps)
    {
        // current radius
        double radius = Coordinates<double>::get_radius(
            iv, m_dx, {m_center[0], m_center[1], m_center[2]});

        // vector of 2 nearest values and radii within the grid
        std::array<double, 2> value_at_point;
        std::array<double, 2> r_at_point;
        // how many units are we from domain boundary?
        int units_from_edge = 0;
        if (a_side == Side::Hi)
        {
            // how many units are we from domain boundary?
            units_from_edge = iv[dir] - m_domain_box.bigEnd(dir);
            // vector of 2 nearest values and radii within the grid
            for (int i = 0; i < 2; i++)
            {
                IntVect iv_tmp = iv;
                iv_tmp[dir] += -units_from_edge - i;
                FOR1(idir)
                {
                    if (iv_tmp[idir] > m_domain_box.bigEnd(idir))
                    {
                        iv_tmp[idir] = m_domain_box.bigEnd(idir);
                    }
                    else if (iv_tmp[idir] < m_domain_box.smallEnd(idir))
                    {
                        iv_tmp[idir] = m_domain_box.smallEnd(idir);
                    }
                }
                value_at_point[i] = out_box(iv_tmp, icomp);
                r_at_point[i] = Coordinates<double>::get_radius(
                    iv_tmp, m_dx, {m_center[0], m_center[1], m_center[2]});
            }
        }
        else // Lo side
        {
            // how many units are we from domain boundary?
            units_from_edge = -iv[dir] + m_domain_box.smallEnd(dir);
            // vector of 2 nearest values within the grid
            for (int i = 0; i < 2; i++)
            {
                IntVect iv_tmp = iv;
                iv_tmp[dir] += units_from_edge + i;
                FOR1(idir)
                {
                    if (iv_tmp[idir] > m_domain_box.bigEnd(idir))
                    {
                        iv_tmp[idir] = m_domain_box.bigEnd(idir);
                    }
                    else if (iv_tmp[idir] < m_domain_box.smallEnd(idir))
                    {
                        iv_tmp[idir] = m_domain_box.smallEnd(idir);
                    }
                }
                value_at_point[i] = out_box(iv_tmp, icomp);
                r_at_point[i] = Coordinates<double>::get_radius(
                    iv_tmp, m_dx, {m_center[0], m_center[1], m_center[2]});
            }
        }

        // assume some radial dependence and fit it
        double analytic_change = 0.0;
        // comp = const
        if (order == 0)
        {
            analytic_change = 0.0;
        }
        // comp = B + A*r
        else if (order == 1)
        {
            double delta_r_in_domain = r_at_point[1] - r_at_point[0];
            double A =
                (value_at_point[1] - value_at_point[0]) / delta_r_in_domain;
            double delta_r_here = radius - r_at_point[0];
            analytic_change = A * delta_r_here;
        }
        // other orders not supported yet
        else
        {
            MayDay::Error("Order not supported for boundary extrapolation.");
        }

        // set the value here to the extrapolated value
        out_box(iv, icomp) = value_at_point[0] + analytic_change;
    }
}

/// Copy the boundary values from src to dest
/// NB only acts if same box layout of input and output data
void BoundaryConditions::copy_boundary_cells(const Side::LoHiSide a_side,
                                             const GRLevelData &a_src,
                                             GRLevelData &a_dest)
{
    CH_TIME("BoundaryConditions::copy_boundary_cells");

    CH_assert(is_defined);
    CH_assert(a_src.nComp() == NUM_VARS);
    if (a_src.boxLayout() == a_dest.boxLayout())
    {
        // cycle through the directions
        FOR1(idir)
        {
            // only do something if this direction is not periodic
            if (!m_params.is_periodic[idir])
            {
                // iterate through the boxes, shared amongst threads
                DataIterator dit = a_dest.dataIterator();
                int nbox = dit.size();
#pragma omp parallel for default(shared)
                for (int ibox = 0; ibox < nbox; ++ibox)
                {
                    DataIndex dind = dit[ibox];
                    FArrayBox &m_dest_box = a_dest[dind];
                    Box this_box = m_dest_box.box();
                    IntVect offset_lo =
                        -this_box.smallEnd() + m_domain_box.smallEnd();
                    IntVect offset_hi =
                        +this_box.bigEnd() - m_domain_box.bigEnd();

                    // reduce box to the intersection of the box and the
                    // problem domain ie remove all outer ghost cells
                    this_box &= m_domain_box;

                    // get the boundary box (may be Empty)
                    Box boundary_box = get_boundary_box(a_side, idir, offset_lo,
                                                        offset_hi, this_box);

                    BoxIterator bit(boundary_box);
                    for (bit.begin(); bit.ok(); ++bit)
                    {
                        IntVect iv = bit();
                        for (int icomp : m_evolution_comps)
                        {
                            m_dest_box(iv, icomp) = a_src[dind](iv, icomp);
                        }
                    } // end iterate over box
                }     // end iterate over boxes
            }         // end if(not periodic)
        }             // end iterate over spacedims
    }                 // end test for same box layout
}

/// Fill the fine boundary values in a_state
/// Required for interpolating onto finer levels at boundaries
void BoundaryConditions::interp_boundaries(GRLevelData &a_fine_state,
                                           GRLevelData &a_coarse_state,
                                           const Side::LoHiSide a_side)
{
    CH_assert(is_defined);
    CH_assert(a_fine_state.nComp() == NUM_VARS);
    CH_assert(a_coarse_state.nComp() == NUM_VARS);
    CH_TIME("BoundaryConditions::interp_boundaries");

    // cycle through the directions
    FOR1(idir)
    {
        // only do something if this direction is not periodic
        if (!m_params.is_periodic[idir])
        {
            // Ref ratio is always two
            int ref_ratio = 2;

            // create a coarsened fine layout and copy the coarse data onto
            // it
            DisjointBoxLayout coarsened_layout;
            coarsen(coarsened_layout, a_fine_state.disjointBoxLayout(),
                    ref_ratio * IntVect::Unit);
            GRLevelData coarsened_fine;
            coarsened_fine.define(coarsened_layout, NUM_VARS,
                                  m_num_ghosts * IntVect::Unit);
            Box coarse_domain_box = coarsen(m_domain_box, ref_ratio);

            // trick the copyTo into thinking the boundary cells are within
            // the domain by growing the domain
            Box grown_domain_box = coarse_domain_box;
            grown_domain_box.grow(m_num_ghosts * IntVect::Unit);
            Copier boundary_copier;
            boundary_copier.ghostDefine(
                a_coarse_state.disjointBoxLayout(),
                coarsened_fine.disjointBoxLayout(), grown_domain_box,
                m_num_ghosts * IntVect::Unit, m_num_ghosts * IntVect::Unit);
            a_coarse_state.copyTo(a_coarse_state.interval(), coarsened_fine,
                                  coarsened_fine.interval(), boundary_copier);

            // iterate through the coarse boxes, shared amongst threads
            DataIterator dit = coarsened_layout.dataIterator();
            int nbox = dit.size();
#pragma omp parallel for default(shared)
            for (int ibox = 0; ibox < nbox; ++ibox)
            {
                DataIndex dind = dit[ibox];
                FArrayBox &m_fine_box = a_fine_state[dind];
                FArrayBox &m_coarse_box = coarsened_fine[dind];
                Box this_box = m_coarse_box.box();
                Box fine_box = m_fine_box.box();
                IntVect offset_lo =
                    -this_box.smallEnd() + coarse_domain_box.smallEnd();
                IntVect offset_hi =
                    +this_box.bigEnd() - coarse_domain_box.bigEnd();

                // reduce box to the intersection of the box and the
                // problem domain ie remove all outer ghost cells
                this_box &= coarse_domain_box;

                // get the boundary box - remove one cell as we only want 2
                // coarse cells filled in each direction, to fill the 3 fine
                // cells on the level above
                Box boundary_box = get_boundary_box(a_side, idir, offset_lo,
                                                    offset_hi, this_box, 1);

                // define standard stencil for interp where not near
                // boundaries in other dirs
                IntVect default_offset =
                    IntVect::Zero + sign(a_side) * 2 * BASISV(idir);
                FourthOrderInterpStencil default_stencil(default_offset,
                                                         ref_ratio);

                // now interp the box from coarse to fine
                BoxIterator bit(boundary_box);
                for (bit.begin(); bit.ok(); ++bit)
                {
                    IntVect iv = bit();
                    IntVect lo_local_offset = iv - m_coarse_box.smallEnd();
                    IntVect hi_local_offset = m_coarse_box.bigEnd() - iv;

                    // bit of work to get the right stencils for near the
                    // edges of the box
                    bool near_boundary = false;
                    IntVect local_boundary_offset = IntVect::Zero;
                    FOR1(idir2)
                    {
                        if (idir2 == idir)
                        {
                            local_boundary_offset[idir2] =
                                default_offset[idir2];
                        }
                        else if ((idir2 != idir) &&
                                 (lo_local_offset[idir2] > 1) &&
                                 (hi_local_offset[idir2] > 1))
                        {
                            local_boundary_offset[idir2] = 0;
                        }
                        else if ((idir2 != idir) &&
                                 (lo_local_offset[idir2] == 1))
                        {
                            local_boundary_offset[idir2] = -2;
                            near_boundary = true;
                        }
                        else if ((idir2 != idir) &&
                                 (hi_local_offset[idir2] == 1))
                        {
                            local_boundary_offset[idir2] = +2;
                            near_boundary = true;
                        }
                        else
                        {
                            MayDay::Error(
                                "BoundaryConditions::define bad boxes");
                        }
                    }

                    // if not near the boundary use the default stencil,
                    // otherwise use the one calculated locally
                    if (!near_boundary)
                    {
                        default_stencil.fillFine(m_fine_box, m_coarse_box, iv);
                    }
                    else
                    {
                        FourthOrderInterpStencil local_stencil(
                            local_boundary_offset, ref_ratio);
                        local_stencil.fillFine(m_fine_box, m_coarse_box, iv);
                    }
                } // end loop box
            }     // end loop boxes
        }         // end if is_periodic
    }             // end loop idir
}

/// Get the boundary condition for a_dir and a_side
int BoundaryConditions::get_boundary_condition(const Side::LoHiSide a_side,
                                               const int a_dir)
{
    int boundary_condition = 0;
    if (a_side == Side::Lo)
    {
        boundary_condition = m_params.lo_boundary[a_dir];
    }
    else
    {
        boundary_condition = m_params.hi_boundary[a_dir];
    }
    return boundary_condition;
}

/// get the boundary box to fill if we are at a boundary
Box BoundaryConditions::get_boundary_box(
    const Side::LoHiSide a_side, const int a_dir, const IntVect &offset_lo,
    const IntVect &offset_hi, Box &this_ghostless_box, int shrink_for_coarse)
{
    // default constructor gives empty box
    Box boundary_box;

    // check if we are over the edges of the domain - are we a boundary box?
    // if so create the box of the cells we want to fill
    if (((a_side == Side::Hi) && (offset_hi[a_dir] > 0)) ||
        ((a_side == Side::Lo) && (offset_lo[a_dir] > 0)))
    {
        // Get just the boundary box to iterate over, m_num_ghosts ghost
        // cells unless we are filling the coarse cells in the interp case
        // where we want to fill only two coarse ghost cells (to cover 3
        // fine ones)
        if (a_side == Side::Lo)
        {
            boundary_box = adjCellLo(this_ghostless_box, a_dir,
                                     m_num_ghosts - shrink_for_coarse);
        }
        else
        {
            boundary_box = adjCellHi(this_ghostless_box, a_dir,
                                     m_num_ghosts - shrink_for_coarse);
        }

        // adjust for any offsets - catches the corners etc
        // but only want to fill them once, so y fills x, z fills y and x
        // etc. Required in periodic direction corners in cases where there
        // are mixed boundaries, (otherwise these corners are full of nans)
        FOR1(idir)
        {
            if (offset_lo[idir] > 0) // this direction is a low end boundary
            {
                if ((idir < a_dir) || (m_params.is_periodic[idir]))
                {
                    // grow it to fill the corners
                    boundary_box.growLo(idir, m_num_ghosts - shrink_for_coarse);
                }
            }
            else // cut off end ghost cell
            {
                if (idir != a_dir)
                {
                    boundary_box.growLo(idir, -shrink_for_coarse);
                }
            }

            if (offset_hi[idir] > 0) // this direction is a high end
                                     // boundary
            {
                if ((idir < a_dir) || (m_params.is_periodic[idir]))
                {
                    // grow it to fill the corners
                    boundary_box.growHi(idir, m_num_ghosts - shrink_for_coarse);
                }
            }
            else // cut off end ghost cell
            {
                if (idir != a_dir)
                {
                    boundary_box.growHi(idir, -shrink_for_coarse);
                }
            }
        }
    }
    return boundary_box;
}

/// Operator called by transform to grow the boxes where required
Box ExpandGridsToBoundaries::operator()(const Box &a_in_box)
{
    Box out_box = a_in_box;
    IntVect offset_lo =
        -out_box.smallEnd() + m_boundaries.m_domain_box.smallEnd();
    IntVect offset_hi = +out_box.bigEnd() - m_boundaries.m_domain_box.bigEnd();

    FOR1(idir)
    {
        if (!m_boundaries.m_params.is_periodic[idir])
        {
            if ((m_boundaries.get_boundary_condition(Side::Lo, idir) ==
                     BoundaryConditions::SOMMERFELD_BC ||
                 m_boundaries.get_boundary_condition(Side::Lo, idir) ==
                     BoundaryConditions::MIXED_BC) &&
                offset_lo[idir] == 0)
            {
                out_box.growLo(idir, m_boundaries.m_num_ghosts);
            }
            if ((m_boundaries.get_boundary_condition(Side::Hi, idir) ==
                     BoundaryConditions::SOMMERFELD_BC ||
                 m_boundaries.get_boundary_condition(Side::Hi, idir) ==
                     BoundaryConditions::MIXED_BC) &&
                offset_hi[idir] == 0)
            {
                out_box.growHi(idir, m_boundaries.m_num_ghosts);
            }
        }
    }
    return out_box;
}

/// This function takes a default constructed open DisjointBoxLayout and
/// grows the boxes lying along the boundary to include the boundaries if
/// necessary (i.e. in the Sommerfeld BC case). It is used to define the correct
/// DisjointBoxLayout for the exchange copier so that shared boundary ghosts are
/// exchanged correctly.
void BoundaryConditions::expand_grids_to_boundaries(
    DisjointBoxLayout &a_out_grids, const DisjointBoxLayout &a_in_grids)
{
    if (!a_in_grids.isClosed())
    {
        MayDay::Error("input to expand_grids_to_boundaries must be closed");
    }

    // Grow the problem domain to include the boundary ghosts
    ProblemDomain domain_with_boundaries = a_in_grids.physDomain();
    FOR1(idir)
    {
        if (!m_params.is_periodic[idir])
        {
            if ((get_boundary_condition(Side::Lo, idir) == SOMMERFELD_BC) ||
                (get_boundary_condition(Side::Lo, idir) == MIXED_BC))
            {
                domain_with_boundaries.growLo(idir, m_num_ghosts);
            }
            if ((get_boundary_condition(Side::Hi, idir) == SOMMERFELD_BC) ||
                (get_boundary_condition(Side::Hi, idir) == MIXED_BC))
            {
                domain_with_boundaries.growHi(idir, m_num_ghosts);
            }
        }
    }

    // Copy grids and apply transformation
    a_out_grids.deepCopy(a_in_grids, domain_with_boundaries);
    ExpandGridsToBoundaries expand_grids_to_boundaries(*this);
    a_out_grids.transform(expand_grids_to_boundaries);
    a_out_grids.close();
}

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(AMRINTERPOLATOR_HPP_)
#error "This file should only be included through AMRInterpolator.hpp"
#endif

#ifndef AMRINTERPOLATOR_IMPL_HPP_
#define AMRINTERPOLATOR_IMPL_HPP_

#include <sstream>

// A bit of Android-ism here, but it's really useful!
// Identifies the printout as originating from this class.
template <typename InterpAlgo>
const string AMRInterpolator<InterpAlgo>::TAG =
    "\x1b[32;1m[AMRInterpolator]\x1b[0m ";

template <typename InterpAlgo>
AMRInterpolator<InterpAlgo>::AMRInterpolator(
    const GRAMR &gr_amr, const std::array<double, CH_SPACEDIM> &coarsest_origin,
    const std::array<double, CH_SPACEDIM> &coarsest_dx, int verbosity)
    : AMRInterpolator(gr_amr, coarsest_origin, coarsest_dx,
                      BoundaryConditions::params_t(), verbosity)
{
}

template <typename InterpAlgo>
AMRInterpolator<InterpAlgo>::AMRInterpolator(
    const GRAMR &gr_amr, const std::array<double, CH_SPACEDIM> &coarsest_origin,
    const std::array<double, CH_SPACEDIM> &coarsest_dx,
    const BoundaryConditions::params_t &a_bc_params, int verbosity)
    : m_gr_amr(gr_amr), m_coarsest_origin(coarsest_origin),
      m_coarsest_dx(coarsest_dx),
      m_num_levels(const_cast<GRAMR &>(m_gr_amr).getAMRLevels().size()),
      m_verbosity(verbosity), m_bc_params(a_bc_params)
{
    set_reflective_BC();
}

template <typename InterpAlgo>
void AMRInterpolator<InterpAlgo>::refresh(const bool a_fill_ghosts)
{
    CH_TIME("AMRInterpolator::refresh");

    m_num_levels = const_cast<GRAMR &>(m_gr_amr).getAMRLevels().size();

    m_mem_level.clear();
    m_mem_box.clear();

    if (a_fill_ghosts)
    {
        fill_multilevel_ghosts(VariableType::evolution);
        if (NUM_DIAGNOSTIC_VARS > 0)
            fill_multilevel_ghosts(VariableType::diagnostic);
    }
}

template <typename InterpAlgo>
void AMRInterpolator<InterpAlgo>::fill_multilevel_ghosts(
    const VariableType a_var_type, const Interval &a_comps,
    const int a_min_level, const int a_max_level)
{
    m_gr_amr.fill_multilevel_ghosts(a_var_type, a_comps, a_min_level,
                                    a_max_level);
}

template <typename InterpAlgo>
const AMR &AMRInterpolator<InterpAlgo>::getAMR() const
{
    return m_gr_amr;
}

template <typename InterpAlgo>
const std::array<double, CH_SPACEDIM> &
AMRInterpolator<InterpAlgo>::get_coarsest_dx() const
{
    return m_coarsest_dx;
}
template <typename InterpAlgo>
const std::array<double, CH_SPACEDIM> &
AMRInterpolator<InterpAlgo>::get_coarsest_origin() const
{
    return m_coarsest_origin;
}

template <typename InterpAlgo>
bool AMRInterpolator<InterpAlgo>::get_boundary_reflective(Side::LoHiSide a_side,
                                                          int a_dir) const
{
    if (a_side == Side::Lo)
        return m_lo_boundary_reflective[a_dir];
    else
        return m_hi_boundary_reflective[a_dir];
}
template <typename InterpAlgo>
bool AMRInterpolator<InterpAlgo>::get_boundary_periodic(int a_dir) const
{
    return m_bc_params.is_periodic[a_dir];
}

template <typename InterpAlgo>
void AMRInterpolator<InterpAlgo>::limit_num_levels(unsigned int num_levels)
{
    CH_TIME("AMRInterpolator::limit_num_levels");

    int max_num_levels = const_cast<AMR &>(m_gr_amr).getAMRLevels().size();
    if (num_levels > max_num_levels || num_levels == 0)
    {
        m_num_levels = max_num_levels;
    }
    else
    {
        m_num_levels = num_levels;
    }
}

// The coordinates of the requested points are most likely stored as
// "struct-of-arrays", i.e. we get an array of x-coords, array of y-coords, etc.
// For maximum compatibility, this function accepts raw pointers for
// coordinates.
template <typename InterpAlgo>
void AMRInterpolator<InterpAlgo>::interp(InterpolationQuery &query)
{
    CH_TIME("AMRInterpolator::interp");

    if (m_verbosity)
    {
        pout() << TAG << "\x1b[32;1mInterpolating data\x1b[0m" << endl;

        for (typename InterpolationQuery::iterator it = query.compsBegin();
             it != query.compsEnd(); ++it)
        {
            const Derivative deriv = it->first;

            pout() << "    This rank is querying for interpolated D(";
            for (int i = 0; i < CH_SPACEDIM; ++i)
            {
                pout() << deriv[i];
                if (i < CH_SPACEDIM - 1)
                {
                    pout() << ",";
                }
            }
            pout() << ") data for " << it->second.size() << " components"
                   << endl;
        }

        pout() << "    Summary: " << query.numComps() << " datasets at "
               << query.m_num_points << " points" << endl;
    }

    // Compute the bounds and spacings for each level
    computeLevelLayouts();

    // Find boxes in which to interpolated each data point
    InterpolationLayout interp_layout = findBoxes(query);

    // Prepare MPI buffers
    prepareMPI(query, interp_layout);
    // Calculate interpolated values
    exchangeMPIQuery();
    calculateAnswers(query);
    exchangeMPIAnswer();

    int comp_idx = 0;

    for (typename InterpolationQuery::iterator deriv_it = query.compsBegin();
         deriv_it != query.compsEnd(); ++deriv_it)
    {
        typedef std::vector<typename InterpolationQuery::out_t> comps_t;
        comps_t &comps = deriv_it->second;

        for (typename comps_t::iterator it = comps.begin(); it != comps.end();
             ++it)
        {
            int comp = std::get<0>(*it);
            double *out = std::get<1>(*it);
            VariableType type = std::get<2>(*it);
            for (int point_idx = 0; point_idx < query.m_num_points; ++point_idx)
            {
                int parity = get_var_parity(comp, type, point_idx, query,
                                            deriv_it->first);
                out[point_idx] =
                    parity * m_query_data[comp_idx][m_mpi_mapping[point_idx]];
            }
            comp_idx++;
        }
    }
}

// Calculates each level's grid spacing and origin location
template <typename InterpAlgo>
void AMRInterpolator<InterpAlgo>::computeLevelLayouts()
{
    CH_TIME("AMRInterpolator::computeLevelLayouts");

    ostream &_pout = pout();

    if (m_verbosity)
    {
        _pout << TAG << "Entering computeLevelLayouts" << endl;
    }

    const Vector<AMRLevel *> &levels =
        const_cast<GRAMR &>(m_gr_amr).getAMRLevels();
    const int num_levels = m_num_levels; // levels.size();

    m_origin.resize(num_levels);
    m_dx.resize(num_levels);

    IntVect prev_small_end = IntVect::Zero;
    IntVect prev_big_end = IntVect::Zero;

    for (int level_idx = 0; level_idx < num_levels; ++level_idx)
    {
        AMRLevel &level = *levels[level_idx];
        const Box &domain_box = level.problemDomain().domainBox();
        const IntVect &small_end = domain_box.smallEnd();
        const IntVect &big_end = domain_box.bigEnd();

        if (level_idx == 0)
        {
            m_origin[level_idx] = m_coarsest_origin;
            m_dx[level_idx] = m_coarsest_dx;
        }
        else
        {
            const int ref_ratio = level.refRatio();

            for (int i = 0; i < CH_SPACEDIM; ++i)
            {
                // Sanity check
                CH_assert(ref_ratio ==
                          (big_end[i] - small_end[i] + 1) /
                              (prev_big_end[i] - prev_small_end[i] + 1));

                // This level's dx
                m_dx[level_idx][i] = m_dx[level_idx - 1][i] / ref_ratio;

                m_origin[level_idx][i] =
                    // coordinate of the previous level's origin
                    m_origin[level_idx - 1][i]
                    // ((edge of the boundary cell) + (distance away from
                    // boundary at the refinement level)) * (expressed in terms
                    // of previous level's grid coordinates)
                    + ((prev_small_end[i] - 0.5) + 1.0 / (2 * ref_ratio)) *
                          m_dx[level_idx - 1][i]
                    // offset from this level's grid origin to small end
                    - small_end[i] * m_dx[level_idx][i];
            }
        }

        if (m_verbosity >= 2)
        {
            _pout << "    Level " << level_idx << "\t"
                  << "dx=("
                  << D_TERM(m_dx[level_idx][0], << "," << m_dx[level_idx][1],
                                                << "," << m_dx[level_idx][2])
                  << ")\t"
                  << "grid_origin=("
                  << D_TERM(m_origin[level_idx][0],
                            << "," << m_origin[level_idx][1],
                            << "," << m_origin[level_idx][2])
                  << ")" << endl;
        }

        prev_small_end = small_end;
        prev_big_end = big_end;
    }

    if (m_verbosity)
    {
        _pout << TAG << "Leaving computeLevelLayouts" << endl;
    }
}

template <typename InterpAlgo>
typename AMRInterpolator<InterpAlgo>::InterpolationLayout
AMRInterpolator<InterpAlgo>::findBoxes(InterpolationQuery &query)
{
    CH_TIME("AMRInterpolator::findBoxes");

    ostream &_pout = pout();

    if (m_verbosity)
    {
        _pout << TAG << "Entering findBoxes" << endl;
    }

    const Vector<AMRLevel *> &levels =
        const_cast<GRAMR &>(m_gr_amr).getAMRLevels();
    const int num_levels = m_num_levels; // levels.size();

    std::array<double, CH_SPACEDIM> grid_coord;
    IntVect nearest;

    InterpolationLayout interp_layout(query.m_num_points);

    int points_found = 0;

    // Check memoised boxes first
    // FIXME: the memoisation gives bad results in the following case:
    // x1 is found on level n-1 in box A, A is memoised
    // x1 is not on level n
    // x2 is in box B on level n and box A on level n-1
    // The current algorithm fill find x2 in A rather than B as it should do
    // Switched off for now
    // for (int i = 0; i < m_mem_box.size(); ++i)
    //{
    //    const int level_idx = m_mem_level[i];
    //    const AMRLevel& level = *levels[level_idx];

    //    const LevelData<FArrayBox>& level_data = dynamic_cast<const
    //    InterpSource<>&>(level).getLevelData(); const DisjointBoxLayout&
    //    box_layout = level_data.disjointBoxLayout(); const Box& domain_box =
    //    level.problemDomain().domainBox();

    //    CH_assert(box_layout.isClosed() && box_layout.isSorted());

    //    const IntVect& small_end = domain_box.smallEnd();
    //    const IntVect& big_end = domain_box.bigEnd();

    //    const int box_idx = m_mem_box[i];
    //    const LayoutIterator& layout_it  = box_layout.layoutIterator();
    //    const LayoutIndex& layout_idx = layout_it[box_idx];
    //    const Box& box = box_layout[layout_idx];
    //    int rank = box_layout.procID(layout_idx);

    //    for (int point_idx = 0; point_idx < query.m_num_points; ++point_idx)
    //    {
    //        // Skip points that have already been found in another box
    //        if (interp_layout.level_idx[point_idx] > -1)
    //        {
    //            continue;
    //        }
    //        else
    //        {
    //            // Calculate "grid coordinates" for current point
    //            for (int i = 0; i < CH_SPACEDIM; ++i)
    //            {
    //                grid_coord[i] = (query.m_coords[i][point_idx] -
    //                m_origin[level_idx][i]) / m_dx[level_idx][i];
    //
    //                // point lies beyond the "small end" of the whole domain,
    //                but still within the boundary cell if (grid_coord[i] <
    //                small_end[i] /*&& grid_coord[i] >= small_end[i] - 0.5*/)
    //                nearest[i] = small_end[i];
    //
    //                // point lies beyond the "big end" of the whole domain,
    //                but still within the boundary cell else if (grid_coord[i]
    //                > big_end[i] /*&& grid_coord[i] <= big_end[i] + 0.5*/)
    //                nearest[i] = big_end[i];
    //
    //                // otherwise we round to nearest grid point
    //                else nearest[i] = (int) ceil(grid_coord[i] - 0.5);
    //            }

    //            if (box.contains(nearest))
    //            {
    //                interp_layout.rank[point_idx] = rank;
    //                interp_layout.level_idx[point_idx] = level_idx;
    //                interp_layout.box_idx[point_idx] = box_idx;

    //                if (m_verbosity >= 2)
    //                {
    //                    _pout << "    Found (";
    //                    for (int i = 0; i < CH_SPACEDIM; ++i)
    //                    {
    //                        _pout << query.m_coords[i][point_idx];
    //                        if (i < CH_SPACEDIM - 1) _pout << ",";
    //                    }
    //                    _pout << ") in level " << level_idx << " box " <<
    //                    box_idx << " (rank " << rank << ")" << endl;
    //                }

    //                points_found += 1;

    //                if (points_found == query.m_num_points)
    //                {
    //                    if (m_verbosity)
    //                    {
    //                        _pout << "    Completed findBoxes using only
    //                        memoised data" << endl;
    //                    }

    //                    // The last remaining socially acceptable use of goto.
    //                    // Cherish it.
    //                    goto found_all_points;
    //                }
    //            }
    //        }
    //    }
    //}

    // Start at the innermost level and increment points_found as we find a box
    // in which to interpolate each point.
    for (int level_idx = num_levels - 1; level_idx >= 0; --level_idx)
    {
        const AMRLevel &level = *levels[level_idx];

        const LevelData<FArrayBox> &level_data =
            dynamic_cast<const InterpSource<> &>(level).getLevelData(
                VariableType::evolution);
        const DisjointBoxLayout &box_layout = level_data.disjointBoxLayout();
        const Box &domain_box = level.problemDomain().domainBox();

        CH_assert(box_layout.isClosed() && box_layout.isSorted());

        const IntVect &small_end = domain_box.smallEnd();
        const IntVect &big_end = domain_box.bigEnd();

        const LayoutIterator &layout_it = box_layout.layoutIterator();

        for (int box_idx = 0; box_idx < box_layout.size(); ++box_idx)
        {
            const LayoutIndex &layout_idx = layout_it[box_idx];

            const Box &box = box_layout[layout_idx];
            int rank = box_layout.procID(layout_idx);

            bool new_box = false;

            for (int point_idx = 0; point_idx < query.m_num_points; ++point_idx)
            {
                // Skip points that have already been found in another box
                if (interp_layout.level_idx[point_idx] > -1)
                {
                    continue;
                }
                else
                {
                    // Calculate "grid coordinates" for current point
                    for (int i = 0; i < CH_SPACEDIM; ++i)
                    {
                        double coord =
                            apply_reflective_BC_on_coord(query, i, point_idx);
                        grid_coord[i] = (coord - m_origin[level_idx][i]) /
                                        m_dx[level_idx][i];

                        // point lies beyond the "small end" of the whole
                        // domain, but still within the boundary cell
                        if (grid_coord[i] <
                            small_end
                                [i] /*&& grid_coord[i] >= small_end[i] - 0.5*/)
                            nearest[i] = small_end[i];

                        // point lies beyond the "big end" of the whole domain,
                        // but still within the boundary cell
                        else if (
                            grid_coord[i] >
                            big_end[i] /*&& grid_coord[i] <= big_end[i] + 0.5*/)
                            nearest[i] = big_end[i];

                        // otherwise we round to nearest grid point
                        else
                            nearest[i] = (int)ceil(grid_coord[i] - 0.5);
                    }

                    if (box.contains(nearest))
                    {
                        interp_layout.rank[point_idx] = rank;
                        interp_layout.level_idx[point_idx] = level_idx;
                        interp_layout.box_idx[point_idx] = box_idx;

                        if (m_verbosity >= 2)
                        {
                            _pout << "    Found (";
                            for (int i = 0; i < CH_SPACEDIM; ++i)
                            {
                                _pout << query.m_coords[i][point_idx];
                                if (i < CH_SPACEDIM - 1)
                                    _pout << ",";
                            }
                            _pout << ") in level " << level_idx << " box "
                                  << box_idx << " (rank " << rank << ")"
                                  << endl;
                        }

                        points_found += 1;
                        new_box = true;

                        if (points_found == query.m_num_points)
                        {
                            m_mem_level.push_back(level_idx);
                            m_mem_box.push_back(box_idx);

                            // The last remaining socially acceptable use of
                            // goto. Cherish it.
                            goto found_all_points;
                        }
                    }
                }
            }

            if (new_box)
            {
                m_mem_level.push_back(level_idx);
                m_mem_box.push_back(box_idx);
            }
        }
    }

found_all_points:
    CH_assert(points_found == query.m_num_points);
    if (m_verbosity)
    {
        _pout << "    All points have been found" << endl;
        _pout << TAG << "Leaving findBoxes" << endl;
    }

    return interp_layout;
}

template <typename InterpAlgo>
void AMRInterpolator<InterpAlgo>::prepareMPI(InterpolationQuery &query,
                                             const InterpolationLayout &layout)
{
    CH_TIME("AMRInterpolator::prepareMPI");

    ostream &_pout = pout();

    if (m_verbosity)
    {
        _pout << TAG << "Entering prepareMPI" << endl;
    }

    // Count the number of points queried to each rank
    m_mpi.clearQueryCounts();
    for (int point_idx = 0; point_idx < query.m_num_points; ++point_idx)
    {
        int rank = layout.rank[point_idx];
        CH_assert(rank > -1);
        m_mpi.incrementQueryCount(rank);
    }

    // Resize MPI 'query' buffers
    m_query_level.resize(query.m_num_points);
    m_query_box.resize(query.m_num_points);
    for (int i = 0; i < CH_SPACEDIM; ++i)
        m_query_coords[i].resize(query.m_num_points);

    m_query_data.resize(query.numComps());
    for (int comp = 0; comp < query.numComps(); ++comp)
        m_query_data[comp].resize(query.m_num_points);

    // Reorder query data for MPI_Ialltoallv
    m_mpi_mapping.resize(query.m_num_points);
    std::vector<int> rank_counter(m_mpi.m_num_process, 0);
    for (int point_idx = 0; point_idx < query.m_num_points; ++point_idx)
    {
        int rank = layout.rank[point_idx];
        int idx = m_mpi.queryDispl(rank) + rank_counter[rank];

        m_query_level[idx] = layout.level_idx[point_idx];
        m_query_box[idx] = layout.box_idx[point_idx];

        for (int i = 0; i < CH_SPACEDIM; ++i)
        {
            // m_query_coords[i][idx] = query.m_coords[i][point_idx];
            m_query_coords[i][idx] =
                apply_reflective_BC_on_coord(query, i, point_idx);
        }

        m_mpi_mapping[point_idx] = idx;

        rank_counter[rank] += 1;
    }

    m_mpi.exchangeLayout();

    // Resize MPI 'answer' buffers
    int num_answers = m_mpi.totalAnswerCount();

    m_answer_level.resize(num_answers);
    m_answer_box.resize(num_answers);
    for (int i = 0; i < CH_SPACEDIM; ++i)
        m_answer_coords[i].resize(num_answers);

    m_answer_data.resize(query.numComps());
    for (int comp = 0; comp < query.numComps(); ++comp)
        m_answer_data[comp].resize(num_answers);

    if (m_verbosity >= 2)
    {
        _pout << "    Number of points that needs to be answered back:" << endl;
        for (int rank = 0; rank < m_mpi.m_num_process; ++rank)
        {
            _pout << "    Rank " << rank << "\t= " << m_mpi.queryCount(rank)
                  << endl;
        }
        _pout << TAG << "Leaving prepareMPI" << endl;
    }
}

template <typename InterpAlgo>
void AMRInterpolator<InterpAlgo>::exchangeMPIQuery()
{
    CH_TIME("AMRInterpolator::exchangeMPIQuery");

    ostream &_pout = pout();

    if (m_verbosity)
    {
        _pout << TAG << "Entering exchangeMPIQuery" << endl;
    }

#ifdef CH_MPI // TODO: it would be nicer if this ifdef were moved into
              // MPIContext ... the only issue is the MPI datatype
    m_mpi.asyncBegin();

    m_mpi.asyncExchangeQuery(&m_query_level[0], &m_answer_level[0], MPI_INT);
    m_mpi.asyncExchangeQuery(&m_query_box[0], &m_answer_box[0], MPI_INT);
    for (int i = 0; i < CH_SPACEDIM; ++i)
    {
        m_mpi.asyncExchangeQuery(&m_query_coords[i][0], &m_answer_coords[i][0],
                                 MPI_DOUBLE);
    }

    m_mpi.asyncEnd();
#else
    m_answer_level = m_query_level;
    m_answer_box = m_query_box;
    for (int i = 0; i < CH_SPACEDIM; ++i)
        m_answer_coords[i] = m_query_coords[i];
#endif

    if (m_verbosity)
    {
        _pout << TAG << "Entering exchangeMPIQuery" << endl;
    }
}

template <typename InterpAlgo>
void AMRInterpolator<InterpAlgo>::calculateAnswers(InterpolationQuery &query)
{
    CH_TIME("AMRInterpolator::calculateAnswers");

    ostream &_pout = pout();

    if (m_verbosity)
    {
        _pout << TAG << "Entering calculateAnswer" << endl;
    }

    const Vector<AMRLevel *> &levels =
        const_cast<GRAMR &>(m_gr_amr).getAMRLevels();
    // const int num_levels = levels.size();
    // const int num_comps = query.numComps();
    const int num_answers = m_mpi.totalAnswerCount();

    std::array<double, CH_SPACEDIM> grid_coord;

    for (int answer_idx = 0; answer_idx < num_answers; ++answer_idx)
    {
        const int box_idx = m_answer_box[answer_idx];
        const int level_idx = m_answer_level[answer_idx];

        const AMRLevel &level = *levels[level_idx];
        const InterpSource<> &source =
            dynamic_cast<const InterpSource<> &>(level);
        const LevelData<FArrayBox> *const evolution_level_data_ptr =
            &source.getLevelData(VariableType::evolution);
        const LevelData<FArrayBox> *diagnostics_level_data_ptr;
        if (NUM_DIAGNOSTIC_VARS > 0)
        {
            diagnostics_level_data_ptr =
                &source.getLevelData(VariableType::diagnostic);
        }
        const DisjointBoxLayout *const evolution_box_layout_ptr =
            &evolution_level_data_ptr->disjointBoxLayout();
        const DisjointBoxLayout *diagnostics_box_layout_ptr;
        if (NUM_DIAGNOSTIC_VARS > 0)
        {
            diagnostics_box_layout_ptr =
                &diagnostics_level_data_ptr->disjointBoxLayout();
        }

        // Convert the LayoutIndex to DataIndex
        const DataIndex evolution_data_idx(
            evolution_box_layout_ptr->layoutIterator()[box_idx]);
        DataIndex diagnostics_data_idx;
        if (NUM_DIAGNOSTIC_VARS > 0)
        {
            diagnostics_data_idx = DataIndex(
                diagnostics_box_layout_ptr->layoutIterator()[box_idx]);
        }

        const Box &box = (*evolution_box_layout_ptr)[evolution_data_idx];
        const FArrayBox *const evolution_fab_ptr =
            &((*evolution_level_data_ptr)[evolution_data_idx]);
        const FArrayBox *diagnostics_fab_ptr;
        if (NUM_DIAGNOSTIC_VARS > 0)
        {
            diagnostics_fab_ptr =
                &((*diagnostics_level_data_ptr)[diagnostics_data_idx]);
        }

        for (int i = 0; i < CH_SPACEDIM; ++i)
        {
            grid_coord[i] =
                (m_answer_coords[i][answer_idx] - m_origin[level_idx][i]) /
                m_dx[level_idx][i];

            if (!(grid_coord[i] >= box.smallEnd()[i] - 0.5 &&
                  grid_coord[i] <= box.bigEnd()[i] + 0.5))
            {
                std::ostringstream s;
                s << "grid_coord[" << i << "] = " << grid_coord[i]
                  << " is out of range [" << (box.smallEnd()[i] - 0.5) << ","
                  << (box.bigEnd()[i] + 0.5) << "]";
                MayDay::Abort(s.str().c_str());
            }
        }

        if (m_verbosity >= 2)
        {
            _pout << "    Interpolating (";
            for (int i = 0; i < CH_SPACEDIM; ++i)
            {
                _pout << m_answer_coords[i][answer_idx];
                if (i < CH_SPACEDIM - 1)
                    _pout << ",";
            }
            _pout << ") in level " << level_idx << " box " << box_idx << endl;
        }

        InterpAlgo algo(source);
        int comp_idx = 0;

        for (typename InterpolationQuery::iterator deriv_it =
                 query.compsBegin();
             deriv_it != query.compsEnd(); ++deriv_it)
        {
            const Derivative deriv = deriv_it->first;

            typedef std::vector<typename InterpolationQuery::out_t> comps_t;
            comps_t &comps = deriv_it->second;
            algo.setup(deriv, grid_coord);

            for (typename comps_t::iterator it = comps.begin();
                 it != comps.end(); ++it)
            {
                int comp = std::get<0>(*it);
                const FArrayBox *fab_ptr;
                if (NUM_DIAGNOSTIC_VARS > 0)
                {
                    VariableType var_type = std::get<2>(*it);
                    fab_ptr = (var_type == VariableType::evolution)
                                  ? evolution_fab_ptr
                                  : diagnostics_fab_ptr;
                }
                else
                {
                    fab_ptr = evolution_fab_ptr;
                }
                double out_val = algo.interpData(*fab_ptr, comp);
                m_answer_data[comp_idx++][answer_idx] = out_val;
            }
        }
    }
}

template <typename InterpAlgo>
void AMRInterpolator<InterpAlgo>::exchangeMPIAnswer()
{
    CH_TIME("AMRInterpolator::exchangeMPIAnswer");

    if (m_verbosity)
    {
        pout() << TAG << "Entering exchangeMPIAnswer" << endl;
    }

#ifdef CH_MPI // TODO: it would be nicer if this ifdef were moved into
              // MPIContext ... the only issue is the MPI datatype
    m_mpi.asyncBegin();

    m_mpi.asyncExchangeAnswer(&m_answer_level[0], &m_query_level[0], MPI_INT);
    m_mpi.asyncExchangeAnswer(&m_answer_box[0], &m_query_box[0], MPI_INT);
    for (int comp = 0; comp < m_answer_data.size(); ++comp)
    {
        m_mpi.asyncExchangeAnswer(&m_answer_data[comp][0],
                                  &m_query_data[comp][0], MPI_DOUBLE);
    }

    m_mpi.asyncEnd();
#else
    m_query_level = m_answer_level;
    m_query_box = m_answer_box;
    m_query_data = m_answer_data;
#endif

    if (m_verbosity)
    {
        pout() << TAG << "Leaving exchangeMPIAnswer" << endl;
    }
}

template <typename InterpAlgo>
void AMRInterpolator<InterpAlgo>::set_reflective_BC()
{
    const IntVect &big_end = const_cast<GRAMR &>(m_gr_amr)
                                 .getAMRLevels()[0]
                                 ->problemDomain()
                                 .domainBox()
                                 .bigEnd();

    FOR(i)
    {
        m_upper_corner[i] = (big_end[i] + 1) * m_coarsest_dx[i];

        m_lo_boundary_reflective[i] =
            (m_bc_params.lo_boundary[i] == BoundaryConditions::REFLECTIVE_BC);
        m_hi_boundary_reflective[i] =
            (m_bc_params.hi_boundary[i] == BoundaryConditions::REFLECTIVE_BC);
    }
}

template <typename InterpAlgo>
int AMRInterpolator<InterpAlgo>::get_var_parity(int comp,
                                                const VariableType type,
                                                int point_idx,
                                                const InterpolationQuery &query,
                                                const Derivative &deriv) const
{
    // check done every time because only one of the variables may be of
    // diagnostic type
    if (type == VariableType::diagnostic &&
        m_bc_params.vars_parity_diagnostic[comp] ==
            BoundaryConditions::UNDEFINED &&
        m_bc_params.reflective_boundaries_exist)
        MayDay::Error("Please provide parameter 'vars_parity_diagnostic' if "
                      "extracting diagnostic variables with reflective BC");

    int parity = 1;
    FOR(dir)
    {
        double coord = query.m_coords[dir][point_idx];
        if ((m_lo_boundary_reflective[dir] && coord < 0.) ||
            (m_hi_boundary_reflective[dir] && coord > m_upper_corner[dir]))
        {

            parity *= BoundaryConditions::get_var_parity(comp, dir, m_bc_params,
                                                         type);
            if (deriv[dir] == 1) // invert parity to first derivatives
                parity *= -1;
        }
    }
    return parity;
}

template <typename InterpAlgo>
double AMRInterpolator<InterpAlgo>::apply_reflective_BC_on_coord(
    const InterpolationQuery &query, double dir, int point_idx) const
{
    double coord = query.m_coords[dir][point_idx];
    if (m_lo_boundary_reflective[dir] && coord < 0.)
        coord = -coord;
    else if (m_hi_boundary_reflective[dir] && coord > m_upper_corner[dir])
        coord = 2. * m_upper_corner[dir] - coord;
    return coord;
}

#endif /* AMRINTERPOLATOR_IMPL_HPP_ */

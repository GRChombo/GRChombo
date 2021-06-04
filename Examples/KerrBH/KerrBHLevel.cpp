/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "KerrBHLevel.hpp"
#include "BoxLoops.hpp"
#include "CCZ4RHS.hpp"
#include "ChiExtractionTaggingCriterion.hpp"
#include "ComputePack.hpp"
#include "KerrBHLevel.hpp"
#include "NanCheck.hpp"
#include "NewConstraints.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SetValue.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"

// Initial data
#include "GammaCalculator.hpp"
#include "KerrBH.hpp"

#include "ADMQuantities.hpp"
#include "ADMQuantitiesExtraction.hpp"

void KerrBHLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

void KerrBHLevel::initialData()
{
    CH_TIME("KerrBHLevel::initialData");
    if (m_verbosity)
        pout() << "KerrBHLevel::initialData " << m_level << endl;

    // First set everything to zero then calculate initial data
    // Get the Kerr solution in the variables, then calculate the \tilde\Gamma^i
    // numerically as these are non zero and not calculated in the Kerr ICs
    BoxLoops::loop(
        make_compute_pack(SetValue(0.), KerrBH(m_p.kerr_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

#ifdef CH_USE_HDF5
void KerrBHLevel::prePlotLevel()
{
    fillAllGhosts();
    BoxLoops::loop(Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}
#endif /* CH_USE_HDF5 */

void KerrBHLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                  const double a_time)
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side
    if (m_p.max_spatial_derivative_order == 4)
    {
        BoxLoops::loop(CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives>(
                           m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation),
                       a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        BoxLoops::loop(CCZ4RHS<MovingPunctureGauge, SixthOrderDerivatives>(
                           m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation),
                       a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

void KerrBHLevel::specificUpdateODE(GRLevelData &a_soln,
                                    const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void KerrBHLevel::preTagCells()
{
    // We only use chi in the tagging criterion so only fill the ghosts for chi
    fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
}

void KerrBHLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                          const FArrayBox &current_state)
{
    BoxLoops::loop(ChiExtractionTaggingCriterion(m_dx, m_level, m_p.max_level,
                                                 m_p.extraction_params,
                                                 m_p.activate_extraction),
                   current_state, tagging_criterion);
}

void KerrBHLevel::specificPostTimeStep()
{
    CH_TIME("KerrBHLevel::specificPostTimeStep");
    // Do the extraction on the min extraction level
    if (m_p.activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_adm = at_level_timestep_multiple(min_level);
        if (calculate_adm)
        {
            // Populate the ADM Mass and Spin values on the grid
            fillAllGhosts();
            BoxLoops::loop(ADMQuantities(m_p.extraction_params.center, m_dx,
                                         c_Madm, c_Jadm),
                           m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);

            if (m_level == min_level)
            {
                CH_TIME("ADMExtraction");
                // Now refresh the interpolator and do the interpolation
                m_gr_amr.m_interpolator->refresh();
                ADMQuantitiesExtraction my_extraction(
                    m_p.extraction_params, m_dt, m_time, m_restart_time, c_Madm,
                    c_Jadm);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }
}
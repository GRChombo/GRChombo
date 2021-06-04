/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMQUANTITIESLEVEL_HPP_
#define ADMQUANTITIESLEVEL_HPP_

#include "ADMQuantities.hpp"
#include "BoxLoops.hpp"
#include "DiagnosticVariables.hpp"
#include "GRAMRLevel.hpp"
#include "SetValue.hpp"

class ADMQuantitiesTestLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<ADMQuantitiesTestLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // initialize data
    virtual void initialData()
    {
        // First set everything to zero then calculate initial data
        // Get the Kerr solution in the variables, then calculate the
        // \tilde\Gamma^i numerically as these are non zero and not calculated
        // in the Kerr ICs
        BoxLoops::loop(
            make_compute_pack(SetValue(0.), KerrBH(m_p.kerr_params, m_dx)),
            m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

        fillAllGhosts();
        // we don't need the Gamma's for the ADM mass calculation
        // BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
        // EXCLUDE_GHOST_CELLS);

        // set ghosts to 0 for extraction
        BoxLoops::loop(SetValue(0.), m_state_diagnostics, m_state_diagnostics,
                       INCLUDE_GHOST_CELLS);
        BoxLoops::loop(
            ADMQuantities(m_p.extraction_params.center, m_dx, c_Madm, c_Jadm),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    }

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time)
    {
    }

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state)
    {
        tagging_criterion.setVal(0.);
    };
};

#endif /* ADMQUANTITIESLEVEL_HPP_ */

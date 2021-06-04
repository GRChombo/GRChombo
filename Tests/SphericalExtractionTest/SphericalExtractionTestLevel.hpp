/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHERICALEXTRACTIONTESTLEVEL_HPP_
#define SPHERICALEXTRACTIONTESTLEVEL_HPP_

#include "BoxLoops.hpp"
#include "GRAMRLevel.hpp"
#include "SetHarmonic.hpp"
#include "SetValue.hpp"
#include "UserVariables.hpp"

class SphericalExtractionTestLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<SphericalExtractionTestLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // initialize data
    virtual void initialData()
    {
        BoxLoops::loop(SetHarmonic(c_phi_Re, c_phi_Im, m_p.es, m_p.el, m_p.em,
                                   m_p.center, m_dx),
                       m_state_new, m_state_new, FILL_GHOST_CELLS);
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

#endif /* SPHERICALEXTRACTIONTESTLEVEL_HPP_ */

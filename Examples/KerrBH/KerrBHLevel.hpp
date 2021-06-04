/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef KERRBHLEVEL_HPP_
#define KERRBHLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"

class KerrBHLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<KerrBHLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    virtual void specificAdvance() override;

    /// Initial data calculation
    virtual void initialData() override;

#ifdef CH_USE_HDF5
    /// Things to do before writing a plot file
    virtual void prePlotLevel() override;
#endif /* CH_USE_HDF5 */

    /// Calculation of the right hand side for the time stepping
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override;

    /// Things to do after dt*rhs has been added to the solution
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    /// Things to do before tagging cells (i.e. filling ghosts)
    virtual void preTagCells() override;

    // to do post each time step on every level
    virtual void specificPostTimeStep() override;

    virtual void
    computeTaggingCriterion(FArrayBox &tagging_criterion,
                            const FArrayBox &current_state) override;
};

#endif /* KERRBHLEVEL_HPP_ */

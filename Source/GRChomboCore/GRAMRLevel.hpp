/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRAMRLEVEL_HPP_
#define GRAMRLEVEL_HPP_

// Chombo includes
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "FourthOrderFillPatch.H"
#include "LevelFluxRegister.H" //We don't actually use flux conservation but Chombo assumes we do
#include "LevelRK4.H"
#include "LoadBalance.H"

// Other includes
#include "BoundaryConditions.hpp"
#include "GRAMR.hpp"
#include "GRLevelData.hpp"
#include "InterpSource.hpp"
#include "SimulationParameters.hpp"
#include "UserVariables.hpp" // need NUM_VARS
#include <fstream>
#include <limits>
#include <sys/time.h>

// Chombo namespace
#include "UsingNamespace.H"

class GRAMRLevel : public AMRLevel, public InterpSource<>
{
  public:
    GRAMRLevel(GRAMR &gr_amr, const SimulationParameters &a_p, int a_verbosity);

    virtual ~GRAMRLevel();

  public:
    /// Do casting from AMRLevel to GRAMRLevel and stop if this isn't possible
    static const GRAMRLevel *gr_cast(const AMRLevel *const amr_level_ptr);
    static GRAMRLevel *gr_cast(AMRLevel *const amr_level_ptr);

    const GRLevelData &
    getLevelData(const VariableType var_type = VariableType::evolution) const;

    bool contains(const std::array<double, CH_SPACEDIM> &point) const;

  private:
    // define
    virtual void define(AMRLevel *a_coarser_level_ptr,
                        const Box &a_problem_domain, int a_level,
                        int a_ref_ratio);

    // define
    virtual void define(AMRLevel *a_coarser_level_ptr,
                        const ProblemDomain &a_problem_domain, int a_level,
                        int a_ref_ratio);

    /// advance by one timestep
    virtual Real advance();

    /// things to do after a timestep
    virtual void postTimeStep();

    /// things to do before tagging cells (e.g. filling ghosts)
    virtual void preTagCells();

    /// tag cells that need to be refined
    virtual void tagCells(IntVectSet &a_tags);

    /// create tags at initialization
    virtual void tagCellsInit(IntVectSet &a_tags);

    /// regrid
    virtual void regrid(const Vector<Box> &a_new_grids);

    /// things to do after regridding
    virtual void postRegrid(int a_base_level);

    /// initialize grids
    virtual void initialGrid(const Vector<Box> &a_new_grids);

    /// things to do after initialization
    virtual void postInitialize();

    /// compute the size of the timestep
    virtual Real computeDt();

    /// compute the size of the initial timestep
    virtual Real computeInitialDt();

    DisjointBoxLayout loadBalance(const Vector<Box> &a_grids);

#ifdef CH_USE_HDF5
    virtual void writeCheckpointHeader(HDF5Handle &a_handle) const;

    virtual void writeCheckpointLevel(HDF5Handle &a_handle) const;

    virtual void readCheckpointHeader(HDF5Handle &a_handle);

    virtual void readCheckpointLevel(HDF5Handle &a_handle);

    virtual void writePlotHeader(HDF5Handle &a_handle) const;

    virtual void writePlotLevel(HDF5Handle &a_handle) const;
#endif

  public:
    /// evaluate d(soln)/dt at current time based on soln
    void evalRHS(GRLevelData &rhs,          //!< d(soln)/dt based on soln
                 GRLevelData &soln,         //!< soln at current time
                 LevelFluxRegister &fineFR, //!< flux register w/ finer level
                 LevelFluxRegister &crseFR, //!< flux register w/ crse level
                 const GRLevelData &oldCrseSoln, //!< old-time crse solution
                 Real oldCrseTime,               //!< old crse time
                 const GRLevelData &newCrseSoln, //!< new-time crse solution
                 Real newCrseTime,               //!< new crse time
                 Real time,      //!< current time centering of soln
                 Real fluxWeight //!< weight to apply to fluxRegister updates
    );

    /// implements soln += dt*rhs
    void updateODE(GRLevelData &soln, const GRLevelData &rhs, Real dt);

    /// define data holder newSoln based on existingSoln, including ghost cell
    /// specification
    void defineSolnData(GRLevelData &newSoln, const GRLevelData &existingSoln);

    /// define data holder for RHS based on existingSoln including ghost cell
    /// specification (which in most cases is no ghost cells)
    void defineRHSData(GRLevelData &newRHS, const GRLevelData &existingSoln);

    /// copy data from src into dest
    void copySolnData(GRLevelData &dest, const GRLevelData &src);

    /// Virtual function for the problem specific parts of Advance
    virtual void specificAdvance() {}

    /// Virtual function for the problem specific parts of postTimeStep
    virtual void specificPostTimeStep() {}

    /// (Pure) virtual function for the initial data calculation
    virtual void initialData() = 0;

    /// Computes which cells have insufficient resolution and should be tagged
    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state) = 0;

#ifdef CH_USE_HDF5
    /// Things to do immediately before checkpointing
    virtual void preCheckpointLevel() {}

    /// Things to do immediately before writing plot files
    virtual void prePlotLevel() {}

    /// Things to do immediately after restart from checkpoint
    virtual void postRestart() {}
#endif

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) = 0;

    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs, Real a_dt)
    {
    }

    // direction irrelevant, but relevant for InterpSource
    ALWAYS_INLINE double get_dx(int dir = 0) const { return m_dx; };

    /// Returns true if m_time is the same as the time at the end of the current
    /// timestep on level a_level and false otherwise
    /// Useful to check whether to calculate something in postTimeStep (which
    /// might only be needed at the end of a_level's timestep)
    bool at_level_timestep_multiple(int a_level) const;

    /// Fill all [either] evolution or diagnostic ghost cells
    virtual void fillAllGhosts(
        const VariableType var_type = VariableType::evolution,
        const Interval &a_comps = Interval(0, std::numeric_limits<int>::max()));

  protected:
    /// Fill all evolution ghosts cells (i.e. those in m_state_new)
    virtual void
    fillAllEvolutionGhosts(const Interval &a_comps = Interval(0, NUM_VARS - 1));

    /// Fill all diagnostics ghost cells (i.e. those in m_state_diagnostics)
    virtual void fillAllDiagnosticsGhosts(
        const Interval &a_comps = Interval(0, NUM_DIAGNOSTIC_VARS - 1));

    /// Fill ghosts cells from boxes on this level only. Do not interpolate
    /// between levels.
    virtual void
    fillIntralevelGhosts(const Interval &a_comps = Interval(0, NUM_VARS - 1));

    /// This function is used to fill ghost cells outside the domain
    /// (for non-periodic boundary conditions, where values depend on state)
    virtual void fillBdyGhosts(GRLevelData &a_state,
                               const Interval &a_comps = Interval(0, NUM_VARS -
                                                                         1));

    /// This function is used to copy ghost cells outside the domain
    /// (for non-periodic boundary conditions, where boundaries evolve via rhs)
    virtual void copyBdyGhosts(const GRLevelData &a_src, GRLevelData &a_dest);

    /// This function is used to define the exchange copiers required for
    /// copying ghost cells between boxes
    virtual void defineExchangeCopier(const DisjointBoxLayout &a_level_domain);

    void printProgress(const std::string &from) const;

    BoundaryConditions m_boundaries; // the class for implementing BCs

    GRLevelData m_state_old; //!< the solution at the old time
    GRLevelData m_state_new; //!< the solution at the new time
    GRLevelData m_state_diagnostics;
    Real m_dx; //!< grid spacing
    double m_restart_time;

    GRAMR &m_gr_amr; //!< The GRAMR object containing this GRAMRLevel

    // params
    SimulationParameters m_p; //!< Holds parameters necessary for the simulation
    int m_verbosity;          //!< Level of verbosity of the output

    Copier m_exchange_copier; //!< copier (for ghost cells on same level)

    CoarseAverage m_coarse_average; //!< Averages from fine to coarse level

    FourthOrderFillPatch m_patcher; //!< Organises interpolation from coarse to
                                    //!< fine levels of ghosts
    FourthOrderFillPatch
        m_patcher_diagnostics; //!< Organises interpolation from coarse to
                               //!< fine levels of ghosts for diagnostics
    FourthOrderFineInterp m_fine_interp; //!< executes the interpolation from
                                         //!< coarse to fine when regridding

    DisjointBoxLayout m_grids;       //!< Holds grid setup (the layout of boxes)
    DisjointBoxLayout m_grown_grids; //!< Holds grown grid setup (for
                                     //!< Sommerfeld BCs)

  public:
    const int m_num_ghosts; //!< Number of ghost cells
};

#endif /* GRAMRLEVEL_HPP_ */

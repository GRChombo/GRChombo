/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef AMRINTERPOLATOR_HPP_
#define AMRINTERPOLATOR_HPP_

// system includes

#include <limits>

// Chombo includes
#include "AMRLevel.H"
#include "LoHiSide.H"

// Our includes
#include "BoundaryConditions.hpp"
#include "GRAMR.hpp"
#include "InterpolationQuery.hpp"
#include "MPIContext.hpp"
#include "UserVariables.hpp"

// Chombo namespace
#include "UsingNamespace.H"

// End include

template <typename InterpAlgo> class AMRInterpolator
{
    struct InterpolationLayout
    {
        std::vector<int> rank;
        std::vector<int> level_idx;
        std::vector<int> box_idx;

        InterpolationLayout(int num_points)
            : rank(num_points, -1), level_idx(num_points, -1),
              box_idx(num_points, -1)
        {
        }
    };

  public:
    // constructor for backward compatibility
    // (adds an artificial BC with only periodic BC)
    AMRInterpolator(const GRAMR &amr,
                    const std::array<double, CH_SPACEDIM> &coarsest_origin,
                    const std::array<double, CH_SPACEDIM> &coarsest_dx,
                    int verbosity = 0);
    AMRInterpolator(const GRAMR &amr,
                    const std::array<double, CH_SPACEDIM> &coarsest_origin,
                    const std::array<double, CH_SPACEDIM> &coarsest_dx,
                    const BoundaryConditions::params_t &a_bc_params,
                    int verbosity = 0);

    void refresh(const bool a_fill_ghosts = true);

    // if not filling ghosts in refresh, call this explicitly for required vars
    void fill_multilevel_ghosts(
        const VariableType a_var_type,
        const Interval &a_comps = Interval(0, std::numeric_limits<int>::max()),
        const int a_min_level = 0,
        const int a_max_level = std::numeric_limits<int>::max());

    void limit_num_levels(unsigned int num_levels);
    void interp(InterpolationQuery &query);
    const AMR &getAMR() const;
    const std::array<double, CH_SPACEDIM> &get_coarsest_dx() const;
    const std::array<double, CH_SPACEDIM> &get_coarsest_origin() const;
    bool get_boundary_reflective(Side::LoHiSide a_side, int a_dir) const;
    bool get_boundary_periodic(int a_dir) const;

  private:
    void computeLevelLayouts();
    InterpolationLayout findBoxes(InterpolationQuery &query);

    void prepareMPI(InterpolationQuery &query,
                    const InterpolationLayout &layout);
    void exchangeMPIQuery();
    void calculateAnswers(InterpolationQuery &query);
    void exchangeMPIAnswer();

    /// set values of member 'm_lo_boundary_reflective' and
    /// 'm_hi_boundary_reflective'
    void set_reflective_BC();
    int get_var_parity(int comp, const VariableType type, int point_idx,
                       const InterpolationQuery &query,
                       const Derivative &deriv) const;
    /// reflect coordinates if BC set to reflective in that direction
    double apply_reflective_BC_on_coord(const InterpolationQuery &query,
                                        double dir, int point_idx) const;

    const GRAMR &m_gr_amr;

    // Coordinates of the point represented by IntVect::Zero in coarsest grid
    const std::array<double, CH_SPACEDIM> m_coarsest_origin;

    // Grid spacing in each direction
    const std::array<double, CH_SPACEDIM> m_coarsest_dx;

    int m_num_levels;
    const int m_verbosity;

    std::vector<std::array<double, CH_SPACEDIM>> m_origin;
    std::vector<std::array<double, CH_SPACEDIM>> m_dx;

    MPIContext m_mpi;
    std::vector<int> m_mpi_mapping;

    // Memoisation of boxes previously found
    std::vector<int> m_mem_level;
    std::vector<int> m_mem_box;

    std::vector<int> m_query_level;
    std::vector<int> m_query_box;
    std::vector<double> m_query_coords[CH_SPACEDIM];
    std::vector<std::vector<double>> m_query_data;

    std::vector<int> m_answer_level;
    std::vector<int> m_answer_box;
    std::vector<double> m_answer_coords[CH_SPACEDIM];
    std::vector<std::vector<double>> m_answer_data;

    // A bit of Android-ism here, but it's really useful!
    // Identifies the printout as originating from this class.
    const static string TAG;

    // Variables for reflective BC
    // m_bc_params can't be a 'const' reference as we need a
    // constructor with backward compatibility that builds an artificial
    // 'BoundaryConditions::params_t'
    BoundaryConditions::params_t m_bc_params;
    /// simplified bools saying whether or not boundary has
    /// a reflective condition in a given direction
    std::array<bool, CH_SPACEDIM> m_lo_boundary_reflective,
        m_hi_boundary_reflective;
    std::array<double, CH_SPACEDIM> m_upper_corner;
};

#include "AMRInterpolator.impl.hpp"

#endif /* AMRINTERPOLATOR_HPP_ */

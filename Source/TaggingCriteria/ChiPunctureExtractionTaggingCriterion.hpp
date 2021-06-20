/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHIPUNCTUREEXTRACTIONTAGGINGCRITERION_HPP_
#define CHIPUNCTUREEXTRACTIONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

//! This class tags cells based on three criteria - the
//! value of the second derivs, the extraction regions
//! and the puncture horizons (which must be covered to
//! a given level
class ChiPunctureExtractionTaggingCriterion
{
  protected:
    const double m_dx;
    const int m_level;
    const int m_max_level;
    const bool m_track_punctures;
    const bool m_activate_extraction;
    const FourthOrderDerivatives m_deriv;
    const spherical_extraction_params_t m_params;
    const std::vector<double> m_puncture_masses;
    const std::vector<std::array<double, CH_SPACEDIM>> &m_puncture_coords;

  public:
    template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
        }
    };

    // The constructor
    ChiPunctureExtractionTaggingCriterion(
        const double dx, const int a_level, const int a_max_level,
        const spherical_extraction_params_t a_params,
        const std::vector<std::array<double, CH_SPACEDIM>> &a_puncture_coords,
        const bool activate_extraction = false,
        const bool track_punctures = false,
        const std::vector<double> a_puncture_masses = {1.0, 1.0})
        : m_dx(dx), m_level(a_level), m_max_level(a_max_level),
          m_track_punctures(track_punctures),
          m_activate_extraction(activate_extraction), m_deriv(dx),
          m_params(a_params), m_puncture_masses(a_puncture_masses),
          m_puncture_coords(a_puncture_coords)
    {
        // check that the number of punctures is consistent
        CH_assert(m_puncture_masses.size() == m_puncture_coords.size());
    };

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // first test the gradients for regions of high curvature
        const auto d2 = m_deriv.template diff2<Vars>(current_cell);
        data_t mod_d2_chi = 0;
        FOR(idir, jdir)
        {
            mod_d2_chi += d2.chi[idir][jdir] * d2.chi[idir][jdir];
        }
        data_t criterion = m_dx * sqrt(mod_d2_chi);

        // if extracting weyl data at a given radius, enforce a given resolution
        // there
        if (m_activate_extraction)
        {
            for (int iradius = 0; iradius < m_params.num_extraction_radii;
                 ++iradius)
            {
                // regrid if within extraction level and not at required
                // refinement
                if (m_level < m_params.extraction_levels[iradius])
                {
                    const Coordinates<data_t> coords(current_cell, m_dx,
                                                     m_params.center);
                    const data_t r = coords.get_radius();
                    // add a 20% buffer to extraction zone so not too near to
                    // boundary
                    auto regrid = simd_compare_lt(
                        r, 1.2 * m_params.extraction_radii[iradius]);
                    criterion = simd_conditional(regrid, 100.0, criterion);
                }
            }
        }

        // ensure that the horizons of the punctures are covered
        // by the max level - for this we need
        // only check the puncture locations on the top 2 levels
        // which regrid (ie, max_level - 1 to max_level - 2)
        // (just the top level would be ok, but doing two ensures
        // the top levels are well spaced)
        if ((m_level > (m_max_level - 3)) && (m_track_punctures == 1))
        {
            // we want each level to be double the innermost one in size
            const double factor = pow(2.0, m_max_level - m_level - 1);
            // loop over puncture masses
            for (int ipuncture = 0; ipuncture < m_puncture_masses.size();
                 ++ipuncture)
            {
                // where am i?
                const Coordinates<data_t> coords(current_cell, m_dx,
                                                 m_puncture_coords[ipuncture]);
                const data_t r = coords.get_radius();
                // decide whether to tag based on distance to horizon
                // plus a fudge factor of 1.5
                auto regrid = simd_compare_lt(
                    r, 1.5 * factor * m_puncture_masses[ipuncture]);
                criterion = simd_conditional(regrid, 100.0, criterion);
            }
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* CHIPUNCTUREEXTRACTIONTAGGINGCRITERION_HPP_ */

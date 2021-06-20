/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHIEXTRACTIONTAGGINGCRITERION_HPP_
#define CHIEXTRACTIONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

//! This class tags cells based on two criteria - the
//! value of the second derivs and the extraction regions
class ChiExtractionTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const spherical_extraction_params_t m_params;
    const int m_level;
    const bool m_activate_extraction;

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
    ChiExtractionTaggingCriterion(const double dx, const int a_level,
                                  const spherical_extraction_params_t a_params,
                                  const bool activate_extraction = false)
        : m_dx(dx), m_deriv(dx), m_params(a_params), m_level(a_level),
          m_activate_extraction(activate_extraction){};
    // below is a constructor for backward compatibility
    ChiExtractionTaggingCriterion(const double dx, const int a_level,
                                  const int a_max_level,
                                  const spherical_extraction_params_t a_params,
                                  const bool activate_extraction = false)
        : ChiExtractionTaggingCriterion(dx, a_level, a_params,
                                        activate_extraction){};

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

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* CHIEXTRACTIONTAGGINGCRITERION_HPP_ */

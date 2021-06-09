/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef LAGRANGE_HPP_
#define LAGRANGE_HPP_

#include "InterpSource.hpp"
#include <utility>

template <int Order, int N_DIMS = CH_SPACEDIM> class Lagrange
{
    const InterpSource<N_DIMS> &m_source;
    bool m_verbosity;

    struct Stencil;

    struct Stencil
    {
        int m_width;
        int m_deriv;
        double m_dx;
        double m_point_offset;

        std::vector<double> m_weights;

        Stencil(int width, int deriv, double point_offset, double dx);
        inline bool operator==(const Stencil &rhs) const;
        inline bool isSameAs(int width, int deriv, double point_offset,
                             double dx) const;

        inline const double &operator[](unsigned int i) const;
    };

    typedef std::vector<Stencil> stencil_collection_t;
    stencil_collection_t m_memoized_stencils;

    Stencil getStencil(int width, int deriv, double point_offset, double dx);

    // Helper function to generate tensor product weights
    // Argument 'dim' is used for recursion over dimensions.
    pair<std::vector<IntVect>, std::vector<double>>
    generateStencil(const std::array<int, N_DIMS> &deriv,
                    const std::array<double, N_DIMS> &dx,
                    const std::array<double, N_DIMS> &eval_index,
                    int dim = N_DIMS - 1);

    std::vector<IntVect> m_interp_points;
    std::vector<double> m_interp_weights;

    // We are adding 216+ numbers at roughly the same magnitudes but alternating
    // signs. Let's keep track of positive and negative terms separately to make
    // sure we don't run into trouble.
    multiset<double> m_interp_neg;
    multiset<double> m_interp_pos;

  public:
    Lagrange(const InterpSource<N_DIMS> &source, bool verbosity = false);

    // eval_index is in 'index' coordinates, not physical coordinates
    void setup(const std::array<int, N_DIMS> &deriv,
               const std::array<double, N_DIMS> &eval_index);

    // any class with a method:
    // Real get(const IntVect &a_iv, int a_comp) const
    template <class GeneralArrayBox>
    double interpData(const GeneralArrayBox &fab, int comp = 0);

    const static string TAG;
};

#include "Lagrange.impl.hpp"

#endif /* LAGRANGE_HPP_ */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef QUINTICCONVOLUTION_HPP_
#define QUINTICCONVOLUTION_HPP_

#include "InterpSource.hpp"
#include <utility>

template <int N_DIMS> class QuinticConvolution
{
    const InterpSource<N_DIMS> &m_source;
    bool m_verbosity;

    std::vector<IntVect> m_interp_points;
    std::vector<double> m_interp_weights;

  public:
    QuinticConvolution(const InterpSource<N_DIMS> &source,
                       bool verbosity = false);

    // eval_index is in 'index' coordinates, not physical coordinates
    void setup(const std::array<int, N_DIMS> &deriv,
               const std::array<double, N_DIMS> &eval_index);

    // any class with a method:
    // Real get(const IntVect &a_iv, int a_comp) const
    template <class GeneralArrayBox>
    double interpData(const GeneralArrayBox &fab, int comp = 0);

    const static string TAG;
};

#include "QuinticConvolution.impl.hpp"

#endif /* QUINTICCONVOLUTION_HPP_ */

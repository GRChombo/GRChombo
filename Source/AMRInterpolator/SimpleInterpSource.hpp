/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMPLEINTERPSOURCE_H_
#define SIMPLEINTERPSOURCE_H_

#include "InterpSource.hpp"

// Simple interpolation source where all grid points are in the same "Box"
// Relevant to be able to use Interpolation classes by themselves
template <int N_DIMS> class SimpleInterpSource : public InterpSource<N_DIMS>
{
    std::array<int, N_DIMS> m_points_per_dir;
    std::array<bool, N_DIMS> m_is_periodic;
    std::array<double, N_DIMS> m_dxs;

    LevelData<FArrayBox> *m_fake;

  public:
    SimpleInterpSource() {}
    SimpleInterpSource(std::array<int, N_DIMS> a_points_per_dir,
                       std::array<double, N_DIMS> a_dxs,
                       std::array<bool, N_DIMS> a_is_periodic = {false})
    {
        set(a_points_per_dir, a_dxs, a_is_periodic);
    }

    void set(std::array<int, N_DIMS> a_points_per_dir,
             std::array<double, N_DIMS> a_dxs,
             std::array<bool, N_DIMS> a_is_periodic = {false})
    {
        m_points_per_dir = a_points_per_dir;
        m_dxs = a_dxs;
        m_is_periodic = a_is_periodic;
    }

    const LevelData<FArrayBox> &
    getLevelData(const VariableType var_type = VariableType::evolution) const
    {
        CH_assert("This should not be needed.");
        return *m_fake;
    }

    bool contains(const std::array<double, N_DIMS> &point) const
    {
        bool in = true;
        for (int idir = 0; idir < N_DIMS; ++idir)
        {
            if (!m_is_periodic[idir] &&
                (point[idir] < 0. || point[idir] > m_points_per_dir[idir] - 1))
                in = false;
        }
        return in;
    };

    double get_dx(int dir) const { return m_dxs[dir]; };
    std::array<double, N_DIMS> get_dxs() const { return m_dxs; }
};

#endif /* SIMPLEINTERPSOURCE_H_ */

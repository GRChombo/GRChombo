/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMPLEARRAYBOX
#define SIMPLEARRAYBOX

// Chombo includes
#include "CH_assert.H"
#include "MayDay.H"

// Chombo namespace
#include "UsingNamespace.H"

// Global Array Box class to interpolate N_DIM dimensional data
// where 'm_f' is a flat array over the dimensions
// Created to replace an FArrayBox in interpolation ::get methods (e.g.
// Lagrange)
template <int N_DIMS> class SimpleArrayBox
{
    // coordinates (u,v,f)
    std::array<int, N_DIMS> m_points_per_dir;
    std::vector<double> m_f;
    std::array<bool, N_DIMS> m_is_periodic;

    int apply_periodic(const IntVect &a_iv, int dir) const
    {
        int idx = a_iv[dir];
        if (idx < 0 || idx >= m_points_per_dir[dir])
        {
            if (m_is_periodic[dir])
            {
                while (idx < 0)
                    idx += m_points_per_dir[dir];
                while (idx >= m_points_per_dir[dir])
                    idx -= m_points_per_dir[dir];
            }
            else
            {
                pout() << "Direction " << dir << " invalid in IntVect " << a_iv
                       << std::endl;
                MayDay::Error("[SimpleArrayBox] Trying to access index out of "
                              "valid range.");
            }
        }
        return idx;
    }

  public:
    // IntVect will be CH_SPACEDIM dimension, but ignore the lasts components if
    // N_DIMS < CH_SPACEDIM
    // the 'comp' argument doesn't matter, always assume 'm_f' is what matters
    Real get(const IntVect &a_iv, int a_comp) const
    {
        int global_idx = apply_periodic(a_iv, N_DIMS - 1);

        for (int i = N_DIMS - 2; i >= 0; --i)
            global_idx =
                global_idx * m_points_per_dir[i] + apply_periodic(a_iv, i);

        return m_f[global_idx];
    }

    SimpleArrayBox() {}

    SimpleArrayBox(std::array<int, N_DIMS> a_points_per_dir,
                   std::vector<double> a_f,
                   std::array<bool, N_DIMS> a_is_periodic = {false})
    {
        set(a_points_per_dir, a_f, a_is_periodic);
    }

    void set(std::array<int, N_DIMS> a_points_per_dir, std::vector<double> a_f,
             std::array<bool, N_DIMS> a_is_periodic = {false})
    {
        m_points_per_dir = a_points_per_dir;
        m_f = a_f;
        m_is_periodic = a_is_periodic;

        CH_assert(N_DIMS <= CH_SPACEDIM);
        int total = 1.;
        for (int i = 0; i < N_DIMS; ++i)
            total *= m_points_per_dir[i];
        CH_assert(m_f.size() == total);
    }
};

#endif /* SIMPLEARRAYBOX */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTERPSOURCE_H_
#define INTERPSOURCE_H_

// Chombo includes
#include "LevelData.H"

// Other inclues
#include "VariableType.hpp"
#include <array>

// Chombo namespace
#include "UsingNamespace.H"

// Abstrace base class to identify points in a grid structure
template <int N_DIMS = CH_SPACEDIM> class InterpSource
{
  public:
    virtual const LevelData<FArrayBox> &getLevelData(
        const VariableType var_type = VariableType::evolution) const = 0;
    virtual bool contains(const std::array<double, N_DIMS> &point) const = 0;
    virtual double get_dx(int dir) const = 0;
    virtual std::array<double, N_DIMS> get_dxs() const
    {
        std::array<double, N_DIMS> dxs;
        for (int i = 0; i < N_DIMS; ++i)
            dxs[i] = get_dx(i);
        return dxs;
    }
};

#endif /* INTERPSOURCE_H_ */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"

class SimulationParameters : public ChomboParameters
{
  public:
    SimulationParameters(GRParmParse &pp) : ChomboParameters(pp) {}
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

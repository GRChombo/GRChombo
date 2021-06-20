/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTEGRATIONMETHODSETUP_HPP_
#define INTEGRATIONMETHODSETUP_HPP_

#include "IntegrationMethod.hpp"

// define the static IntegrationMethods here
// Closed methods:
const IntegrationMethod IntegrationMethod::trapezium({0.5});
const IntegrationMethod IntegrationMethod::simpson({0.3333333333333333,
                                                    1.3333333333333333});
const IntegrationMethod IntegrationMethod::simpson38({0.375, 1.125, 1.125});
const IntegrationMethod
    IntegrationMethod::boole({0.3111111111111111, 1.4222222222222222,
                              0.53333333333333, 1.4222222222222222});

// Open Methods:
const IntegrationMethod IntegrationMethod::midpoint({1.0}, false);
const IntegrationMethod IntegrationMethod::open_2nd_order({1.125, 0.75, 1.125},
                                                          false);
const IntegrationMethod
    IntegrationMethod::open_3rd_order({1.0833333333333333, 0.9166666666666666,
                                       0.9166666666666666, 1.0833333333333333},
                                      false);
const IntegrationMethod IntegrationMethod::open_4th_order(
    {1.1935763888888888, 0.4340277777777778, 1.7447916666666667,
     0.4340277777777778, 1.1935763888888888},
    false);

#endif /* INTEGRATIONMETHODSETUP_HPP_ */

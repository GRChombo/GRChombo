/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "EmptyDiagnosticVariables.hpp"
#include <array>
#include <string>

// assign enum to each variable
enum
{
    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names = {};
}

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */

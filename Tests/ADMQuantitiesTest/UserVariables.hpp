/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4UserVariables.hpp"
#include "DiagnosticVariables.hpp"

/// This enum gives the index of every variable stored in the grid
enum
{
    // Note that it is important that the first enum value is set to 1 more than
    // the last CCZ4 var enum
    NUM_VARS = NUM_CCZ4_VARS,
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names =
    ccz4_variable_names;
} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */

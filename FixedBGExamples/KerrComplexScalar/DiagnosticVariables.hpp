/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_chi,    // the conformal factor of the metric, not evolved
    c_rho,    // the energy density of the SF
    c_rhoJ,   // the energy density of the SF
    c_Edot,   // the energy density of the SF
    c_Jdot,   // the energy density of the SF
    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi", "rho", "rhoJ", "Edot", "Jdot"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */

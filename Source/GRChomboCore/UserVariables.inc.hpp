/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(USERVARIABLES_HPP)
#error "This file should only be included through UserVariables.hpp"
#endif

#ifndef USERVARIABLES_INC_HPP
#define USERVARIABLES_INC_HPP

#include "GRParmParse.hpp"
#include "VariableType.hpp"
#include "parstream.H"
#include <algorithm>
#include <array>
#include <string>

// This file must be included at the end of UserVariables.hpp

namespace UserVariables
{
/// Takes a string and returns the variable enum number if the string
/// matches one of those in UserVariables::variable_names, or returns -1
/// otherwise
static int variable_name_to_enum(const std::string &a_var_name)
{
    const auto var_name_it =
        std::find(variable_names.begin(), variable_names.end(), a_var_name);

    int var = std::distance(variable_names.begin(), var_name_it);
    if (var != NUM_VARS)
        return var;
    else
    {
        return -1;
    }
}

} // namespace UserVariables

namespace DiagnosticVariables
{
/// Takes a string and returns the variable enum number if the string
/// matches one of those in UserVariables::variable_names, or returns -1
/// otherwise
static int variable_name_to_enum(const std::string &a_var_name)
{
    const auto var_name_it =
        std::find(variable_names.begin(), variable_names.end(), a_var_name);

    int var = std::distance(variable_names.begin(), var_name_it);
    if (var != NUM_DIAGNOSTIC_VARS)
        return var;
    else
    {
        return -1;
    }
}

} // namespace DiagnosticVariables

namespace UserVariables
{
// where one has read in a subset of variables with some feature
// this reads in a set of associated values and assigns it into a full
// array of all NUM_VARS vars (setting other values to a default value)
template <class T>
void load_values_to_array(
    GRParmParse &pp, const char *a_values_vector_string,
    const std::vector<std::pair<int, VariableType>> &a_vars_vector,
    std::array<T, NUM_VARS> &a_values_array, const T a_default_value)
{
    // how many values do I need to get?
    int num_values = a_vars_vector.size();
    // make a container for them, and load
    std::vector<T> vars_values(num_values, a_default_value);
    pp.load(a_values_vector_string, vars_values, num_values, vars_values);

    // populate the values_array for the NUM_VARS values with those read in
    a_values_array.fill(a_default_value);
    for (int i = 0; i < num_values; i++)
    {
        int icomp = a_vars_vector[i].first;
        CH_assert(a_vars_vector[i].second == VariableType::evolution);
        a_values_array[icomp] = vars_values[i];
    }
}

// function to create a vector of enums of vars by reading in their
// names as strings from the params file and converting it to the enums
inline void
load_vars_to_vector(GRParmParse &pp, const char *a_vars_vector_string,
                    const char *a_vector_size_string,
                    std::vector<std::pair<int, VariableType>> &a_vars_vector,
                    int &a_vars_vector_size)
{
    int num_values;
    pp.load(a_vector_size_string, num_values, -1);
    // only set a_vars_vector and a_var_vector_size if a_vector_size_string
    // found
    if (num_values >= 0)
    {
        std::vector<std::string> var_names(num_values, "");
        pp.load(a_vars_vector_string, var_names, num_values, var_names);
        for (std::string var_name : var_names)
        {
            // first assume the variable is a normal evolution var
            int var = UserVariables::variable_name_to_enum(var_name);
            VariableType var_type = VariableType::evolution;
            if (var < 0)
            {
                // if not an evolution var check if it's a diagnostic var
                var = DiagnosticVariables::variable_name_to_enum(var_name);
                if (var < 0)
                {
                    // it's neither :(
                    pout() << "Variable with name " << var_name << " not found."
                           << endl;
                }
                else
                {
                    var_type = VariableType::diagnostic;
                }
            }
            if (var >= 0)
            {
                a_vars_vector.emplace_back(var, var_type);
            }
        }
        // overwrites read in value if entries have been ignored
        a_vars_vector_size = a_vars_vector.size();
    }
}

} // namespace UserVariables

#endif /* USERVARIABLES_INC_HPP_ */

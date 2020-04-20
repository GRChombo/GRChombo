/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERRESTART_HPP_
#define USERRESTART_HPP_

#include "SPMD.H"
#include <unistd.h>

// This class includes functions that enables the user to trigger a clean
// restart mid evolution

class UserRestart
{
  public:
    // checks whether the restart trigger file exists and updates s_activate
    static void check(std::string a_restart_trigger_file);

    // removes the restart trigger file
    static void remove_trigger(std::string a_restart_trigger_file);

    static bool activate();

  private:
    static bool s_activate;
};

#endif /* USERRESTART_HPP_ */

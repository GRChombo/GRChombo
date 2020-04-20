/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "UserRestart.hpp"

// This file includes functions that enables the user to trigger a clean restart
// mid evolution

bool UserRestart::s_activate = false;

void UserRestart::check(std::string a_restart_trigger_file)
{
    int restart = 0;
    int master_rank = 0;

    // check if the restart trigger file exists on rank 0
    if (procID() == master_rank)
    {
        restart = static_cast<int>(
            access(a_restart_trigger_file.c_str(), F_OK) != -1);
    }

    // broadcast check to all ranks
    broadcast(restart, master_rank);
    s_activate = static_cast<bool>(restart);
}

void UserRestart::remove_trigger(std::string a_restart_trigger_file)
{
    if (procID() == 0)
        ;
    {
        bool trigger_file_exists =
            (access(a_restart_trigger_file.c_str(), F_OK) != -1);
        if (trigger_file_exists)
        {
            std::remove(a_restart_trigger_file.c_str());
        }
    }
    s_activate = false;
}

bool UserRestart::activate() { return s_activate; }

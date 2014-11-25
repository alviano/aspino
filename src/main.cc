/*
 *  Copyright (C) 2014  Mario Alviano (mario@alviano.net)
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 */


#include "Solver.h"
#include "SatSolver.h"
#include "PseudoBooleanSolver.h"
#include "MaxSatSolver.h"
#include "utils/trace.h"

#include <errno.h>
#include <signal.h>
#include <zlib.h>

#include <gflags/gflags.h>
#include <string>
#include <iostream>

using namespace aspino;
using namespace std;

static bool validate_mode(const char* name, const string& value) {
    if(value == "asp") return true;
    if(value == "sat") return true;
    if(value == "maxsat") return true;
    if(value == "pbs") return true;
    cerr << "Invalid value for --" << name << ": " << value << "\n";
    return false;
}
DEFINE_string(mode, "asp", "How to interpret input. Valid values: asp, sat, maxsat, pbs.");
DEFINE_int32(n, 1, "Number of desired solutions. Non-positive integers are interpreted as unbounded.");

extern bool validate_maxsat_strat(const char* name, const string& value);

static aspino::Solver* solver;

static void SIGINT_interrupt(int) { solver->interrupt(); }

static void SIGINT_exit(int) { solver->exit(1); }

int main(int argc, char** argv)
{
    // Use signal handlers that forcibly quit until the solver will be able to respond to
    // interrupts:
    signal(SIGINT, SIGINT_exit);
    signal(SIGXCPU,SIGINT_exit);

    // Change to signal-handlers that will only notify the solver and allow it to terminate
    // voluntarily:
    signal(SIGINT, SIGINT_interrupt);
    signal(SIGXCPU,SIGINT_interrupt);

    google::SetUsageMessage(string()
        + "Solve ASP or SAT problems read from STDIN.\n"
        + "usage: " + argv[0] + " [flags]");

    google::RegisterFlagValidator(&FLAGS_mode, &validate_mode);
    google::RegisterFlagValidator(&FLAGS_maxsat_strat, &validate_maxsat_strat);
    google::ParseCommandLineFlags(&argc, &argv, true);
    
    if(FLAGS_mode == "asp")
//        solver = new AspSolver();
        { cerr << "ASP is not currently supported." << endl; exit(-1); }
    else if(FLAGS_mode == "sat")
        solver = new SatSolver();
    else if(FLAGS_mode == "maxsat")
        solver = new MaxSatSolver();
    else if(FLAGS_mode == "pbs")
        solver = new PseudoBooleanSolver();
    else
        assert(0);
    
    gzFile in = gzdopen(0, "rb");
    solver->parse(in);
    gzclose(in);

    solver->eliminate(true);
    if(!solver->okay()) {
        cout << "UNSATISFIABLE" << endl;
        solver->exit(20);
    }
    
    lbool ret = solver->solve(FLAGS_n);
    
#ifdef NDEBUG
    solver->exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
#else
    delete solver;
    return (ret == l_True ? 10 : ret == l_False ? 20 : 0);
#endif
}
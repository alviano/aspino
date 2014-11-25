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

int main(int argc, char** argv)
{
    signal(SIGINT, SIGINT_interrupt);
    signal(SIGTERM, SIGINT_interrupt);

    google::SetUsageMessage(string()
        + "Solve ASP or SAT problems read from STDIN or provided as command-line argument.\n"
        + "usage: " + argv[0] + " [flags] [input-file]");

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
    
    if(argc > 2) {
        cerr << "Extra argument: " << argv[2] << endl;
        solver->exit(-1);
    }

    gzFile in = argc == 1 ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
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
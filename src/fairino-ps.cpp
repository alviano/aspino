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


#include "AbstractSolver.h"
#include "SatSolver.h"
#include "PseudoBooleanSolver.h"
#include "MaxSatSolver.h"
#include "FairSatSolver.h"
#include "AspSolver.h"
#include "utils/trace.h"

#include <utils/Options.h>

#include <errno.h>
#include <signal.h>
#include <zlib.h>

#include <string>
#include <iostream>
#if defined(__linux__)
#include <fpu_control.h>
#endif

using namespace aspino;
using namespace std;

static Glucose::EnumOption option_mode("MAIN", "mode", "How to interpret input.\n", "asp|sat|maxsat|pbs|fairsat");
Glucose::IntOption option_n("MAIN", "n", "Number of desired solutions. Non-positive integers are interpreted as unbounded.\n", 1, Glucose::IntRange(0, INT32_MAX));

extern Glucose::EnumOption option_fairsat_alg;

static aspino::AbstractSolver* solver;

static void SIGINT_interrupt(int) { solver->interrupt(); }

int main(int argc, char** argv)
{
    signal(SIGINT, SIGINT_interrupt);
    signal(SIGTERM, SIGINT_interrupt);

    Glucose::setUsageHelp(
        "Solve ASP or SAT problems read from STDIN or provided as command-line argument.\n\n"
        "usage: %s [flags] [input-file]\n");

#if defined(__linux__)
    fpu_control_t oldcw, newcw;
    _FPU_GETCW(oldcw); newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(newcw);
    //printf("c WARNING: for repeatability, setting FPU to use double precision\n");
#endif

    Glucose::IntOption cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", INT32_MAX, Glucose::IntRange(0, INT32_MAX));
    Glucose::IntOption mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", INT32_MAX, Glucose::IntRange(0, INT32_MAX));

//    Glucose::parseOptions(argc, argv, true);
    option_mode = "fairsat";
    option_fairsat_alg = "progression";

    if(strcmp(option_mode, "asp") == 0)
        solver = new AspSolver();
    else if(strcmp(option_mode, "sat") == 0)
        solver = new SatSolver();
    else if(strcmp(option_mode, "maxsat") == 0)
        solver = new MaxSatSolver();
    else if(strcmp(option_mode, "pbs") == 0)
        solver = new PseudoBooleanSolver();
    else if(strcmp(option_mode, "fairsat") == 0)
        solver = new FairSatSolver();
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
    
    lbool ret = solver->solve(option_n);
    
#ifndef SAFE_EXIT
    solver->exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
#else
    delete solver;
    return (ret == l_True ? 10 : ret == l_False ? 20 : 0);
#endif
}
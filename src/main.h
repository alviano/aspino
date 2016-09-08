/*
 *  Copyright (C) 2016  Mario Alviano (mario@alviano.net)
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

#ifndef __main_h__
#define __main_h__

#include "AbstractSolver.h"
#include "SatSolver.h"
#include "PseudoBooleanSolver.h"
#include "MaxSatSolver.h"
#include "AspSolver.h"
#include "FairSatSolver.h"
#include "LTLSolver.h"
#include "TGDsSolver.h"
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

extern Glucose::EnumOption option_mode;
extern Glucose::IntOption option_n;

extern Glucose::IntOption cpu_lim;
extern Glucose::IntOption mem_lim;

extern aspino::AbstractSolver* solver;

void SIGINT_interrupt(int);

void premain();
int postmain(int argc, char** argv);


#endif //__main_h__
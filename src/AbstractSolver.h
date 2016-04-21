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

#ifndef __AbstractSolver_h__
#define __AbstractSolver_h__

#include <zlib.h>
#include <ostream>

using namespace std;

#include "utils/assert.h"
#include "utils/trace.h"

#include "core/Solver.h"
#include "core/SolverTypes.h"

namespace aspino {

using Glucose::lbool;
using Glucose::Lit;
using Glucose::Var;
using Glucose::mkLit;
using Glucose::vec;
using Glucose::CRef;
using Glucose::CRef_Undef;
using Glucose::lit_Undef;
using Glucose::ClauseAllocator;
using Glucose::Clause;

ostream& operator<<(ostream& o, lbool l);
ostream& operator<<(ostream& o, const Lit& l);
ostream& operator<<(ostream& o, const Clause& c);
template <class T> ostream& operator<<(ostream& o, const vec<T>& v) {
    o << "[ ";
    for(int i = 0; i < v.size(); i++)
        o << v[i] << " ";
    return o << "]";
}

class AbstractSolver {
public:
    AbstractSolver() {}
    virtual ~AbstractSolver() {}
    
    virtual void interrupt() = 0;
    virtual void exit(int code) { _exit(code); }

    virtual void parse(gzFile in) = 0;
    virtual bool simplify() = 0;
    virtual bool eliminate(bool turn_off_elim) = 0;
    virtual bool okay() const = 0;
    virtual lbool solve() = 0;
    virtual lbool solve(int n) = 0;
};
    
} // namespace aspino

#endif
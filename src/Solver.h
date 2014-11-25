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

#ifndef __Solver_h__
#define __Solver_h__

#include <zlib.h>
#include <ostream>

using namespace std;

#include "utils/assert.h"
#include "utils/trace.h"

#include "core/Solver.h"
#include "core/SolverTypes.h"

namespace aspino {

#define Minisat Glucose
    
using Minisat::lbool;
using Minisat::Lit;
using Minisat::Var;
using Minisat::mkLit;
using Minisat::vec;
using Minisat::CRef;
using Minisat::CRef_Undef;
using Minisat::lit_Undef;
using Minisat::ClauseAllocator;
using Minisat::Clause;

ostream& operator<<(ostream& o, const Lit& l);
template <class T> ostream& operator<<(ostream& o, const vec<T>& v) {
    o << "[ ";
    for(int i = 0; i < v.size(); i++)
        o << v[i] << " ";
    return o << "]";
}

class Solver {
public:
    Solver() {}
    virtual ~Solver() {}
    
    virtual void interrupt() = 0;
    virtual void exit(int code) { _exit(code); }

    virtual void parse(gzFile in) = 0;
    virtual bool simplify() = 0;
    virtual bool eliminate(bool turn_off_elim) = 0;
    virtual bool okay() const = 0;
    virtual lbool solve() = 0;
    virtual lbool solve(int n) = 0;
};
    
} // namespace miniasp

#endif
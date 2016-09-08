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

#ifndef __TGDsSolver_h__
#define __TGDsSolver_h__

#include "SatSolver.h"

namespace aspino {

class TGDsSolver : public SatSolver {
public:
    class Listener {
    public:
        virtual ~Listener() {}
        virtual void onTrueCandidate(Lit candidate, int index);
        virtual void onFalseCandidate(Lit candidate, int index);
    };

    TGDsSolver();
    virtual ~TGDsSolver() {}
    
    void setListener(Listener* value) { assert(value != NULL); listener = value; }
    
    void addCandidate(Lit lit) { candidates.push(lit); }
    
    void clearAssumptions() { assumptions.clear(); }
    void addAssumption(Lit lit) { assumptions.push(lit); }

    void disprove();
    
    virtual void parse(gzFile in);
    virtual lbool solve() { return SatSolver::solve(); }
    virtual lbool solve(int) { return l_True; }

private:
    vec<Lit> candidates;
    Listener* listener;
    
    static Listener defaultListener;
};

} // namespace miniasp

#endif

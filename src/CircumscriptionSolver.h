/*
 *  Copyright (C) 2017  Mario Alviano (mario@alviano.net)
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

#ifndef __CircumscriptionSolver_h__
#define __CircumscriptionSolver_h__

#include "PseudoBooleanSolver.h"
#include "DepGraph.h"

namespace aspino {

class CircumscriptionSolver : public PseudoBooleanSolver {
public:
    CircumscriptionSolver();
    virtual ~CircumscriptionSolver();
    
    virtual void interrupt();
    
    virtual void parse(gzFile in);
    virtual bool eliminate(bool) { return SatSolver::eliminate(false); } 
    
    virtual lbool solve() { return solve(1); }
    virtual lbool solve(int);
    
    inline int nOLits() const { return objectLits.size(); }
    inline Lit getOLit(int idx) const { assert(idx >= 0 && idx < objectLits.size()); return objectLits[idx]; }
    
    inline int nCLits() const { return careLits.size(); }
    inline Lit getCLit(int idx) const { assert(idx >= 0 && idx < careLits.size()); return careLits[idx]; }

private:
    bool checkConsistencyAndRemoveDuplicates(vec<Lit>& lits);
    bool disjoint(const vec<Lit>& a, const vec<Lit>& b);
    void onComputedModel();
    void extendModel();
    
    void enumerate();
    void enumerateByBlockingClauses();
    void trim();
    void shrink();
    void one();
    void kdyn();
    void pmres();
    
    lbool sat() const;
    lbool unsat() const;

    int requiredModels;
    int requiredWitnesses;
    int computedModels;

    vec<Lit> objectLits;
    vec<Lit> careLits;
    vec<Lit> visible;
    vec<char*> visibleValue;
    
    vec<Lit> softLits;
    vec<bool> removed;
    vec<Lit> nonInputSoftLits;
};


} // namespace aspino

#endif
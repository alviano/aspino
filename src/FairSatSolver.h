/*
 *  Copyright (C) 2015  Mario Alviano (mario@alviano.net)
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

#ifndef __FairSatSolver_h__
#define __FairSatSolver_h__

#include "MaxSatSolver.h"

namespace aspino {

class ObjectFunction {
public:
    inline ObjectFunction() {}
//    inline ObjectFunction(ObjectFunction& init) : sumOfInCoeffs(init.sumOfInCoeffs) { init.lits.moveTo(lits); init.coeffs.moveTo(coeffs); init.selectorVars.moveTo(selectorVars); }
    ~ObjectFunction() { lits.clear(); coeffs.clear(); selectorVars.clear(); }
    
    void init(vec<Lit>& lits, vec<int64_t>& coeffs);
    
    
//private:
    vec<Lit> lits;
    vec<int64_t> coeffs;
    int64_t sumOfInCoeffs;
    vec<Var> selectorVars;
    int64_t modelValue;
    
private:
    ObjectFunction(ObjectFunction&) {}
};
    
class FairSatSolver : public PseudoBooleanSolver {
public:
    FairSatSolver();
    virtual ~FairSatSolver();
    
    virtual void interrupt() { cout << "s UNKNOWN" << endl; this->exit(1); }
    
    virtual void parse(gzFile in);
    
    void addObjectFunction(vec<Lit>& lits, vec<int64_t>& coeffs);
    int64_t processObjectFunctions();
    
    void setMinObjectFunction(int64_t min);
    void updateLowerBound();
    
    virtual lbool solve();
    virtual lbool solve(int) { return solve(); }
    
private:
    vec<ObjectFunction*> objectFunctions;
    int64_t lowerbound;
    int64_t upperbound;
    
    void (FairSatSolver::*search_alg)();
    void search_alg_bb();
    void search_alg_binary();
    void search_alg_progression();
};


} // namespace aspino

#endif
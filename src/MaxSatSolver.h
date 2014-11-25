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

#ifndef __MaxSatSolver_h__
#define __MaxSatSolver_h__

#include "PseudoBooleanSolver.h"

DECLARE_string(maxsat_strat);
DECLARE_bool(maxsat_reitereted_disjoint_cores);

namespace aspino {

class MaxSatSolver : public PseudoBooleanSolver {
public:
    MaxSatSolver();
    virtual ~MaxSatSolver();
    
    virtual void parse(gzFile in);
    
    void addWeightedClause(vec<Lit>& lits, long weight);
    
    virtual lbool solve();
    virtual lbool solve(int) { return solve(); }
    
    
private:
    long upperbound;
    long lowerbound;
    vec<Lit> softLiterals;
    vec<long> weights;
    int lastSoftLiteral;
    
    void sort();
    void detectLevels();
    long setAssumptions(long limit);
    lbool solve_(long limit);
    
    void (MaxSatSolver::*corestrat)(long);
    void corestrat_one(long limit);
    void corestrat_pmres(long limit);
    void corestrat_pmres_split_conj(long limit);
    void corestrat_pmreslog(long limit);
};



} // namespace miniasp

#endif
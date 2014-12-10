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
DECLARE_string(maxsat_disjcores);
DECLARE_bool(maxsat_saturate);

namespace aspino {

class MaxSatSolver : public PseudoBooleanSolver {
public:
    MaxSatSolver();
    virtual ~MaxSatSolver();
    
    virtual void interrupt() { cout << "s UNKNOWN" << endl; this->exit(1); }
    
    virtual void parse(gzFile in);
    
    void addWeightedClause(vec<Lit>& lits, long weight);
    
    virtual lbool solve();
    virtual lbool solve(int) { return solve(); }
    
    
private:
    long upperbound;
    long lowerbound;
    
    vec<Lit> softLiterals;
    vec<long> weights;
    
    vec<vec<Lit>*> levels;
    vec<long> weightOfPreviousLevel;
    
    int lastSoftLiteral;
    long firstLimit;
    
    void sameSoftVar(Lit soft, long weight);
    
    void detectLevels();
    long setAssumptions(long limit);
    lbool solveCurrentLevel();
    lbool solve_();
    
    void updateUpperBound();
    
    void trim();
    
    void (MaxSatSolver::*corestrat)(long);
    void corestrat_one(long limit);
    void corestrat_pmres(long limit);
    void corestrat_pmres_split_conj(long limit);
    void corestrat_pmreslog(long limit);

    enum DisjunctCores {NO = 0, PRE, ALL};
    DisjunctCores disjcores;
    
    bool saturate;
};



} // namespace miniasp

#endif
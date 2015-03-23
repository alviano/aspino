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

namespace aspino {

class MaxSatSolver : public PseudoBooleanSolver {
public:
    MaxSatSolver();
    virtual ~MaxSatSolver();
    
    virtual void interrupt() { cout << "s UNKNOWN" << endl; this->exit(1); }
    
    virtual void parse(gzFile in);
    
    void addWeightedClause(vec<Lit>& lits, int64_t weight);
    
    virtual lbool solve();
    virtual lbool solve(int) { return solve(); }
    
    
private:
    int64_t upperbound;
    int64_t lowerbound;
    
    vec<Lit> softLiterals;
    vec<int64_t> weights;
    
    vec<vec<Lit>*> levels;
    vec<int64_t> weightOfPreviousLevel;
    
    int lastSoftLiteral;
    int64_t firstLimit;
    
    void sameSoftVar(Lit soft, int64_t weight);
    
    void removeSoftLiteralsAtLevelZero();
    void detectLevels();
    void preprocess();
    int64_t setAssumptions(int64_t limit);
    void solveCurrentLevel();
    void solve_();
    
    void updateUpperBound();
    
    void trim();
    void minimize();
    
    void (MaxSatSolver::*corestrat)(int64_t);
    void corestrat_one(int64_t limit);
    void corestrat_one_2(int64_t limit);
    void corestrat_one_wc(int64_t limit);
    void corestrat_one_neg(int64_t limit);
    void corestrat_one_neg_wc(int64_t limit);
    void corestrat_one_pmres(int64_t limit);
    void corestrat_one_pmres_2(int64_t limit);
    void corestrat_pmres(int64_t limit);
    void corestrat_pmres_reverse(int64_t limit);
    void corestrat_pmres_split_conj(int64_t limit);
    void corestrat_pmreslog(int64_t limit);

    enum DisjunctCores {NO = 0, PRE, ALL};
    DisjunctCores disjcores;
    
    bool saturate;
    
    int numberOfCores;
    int sizeOfCores;
    vec<vec<Lit>*> cores;
};



} // namespace aspino

#endif
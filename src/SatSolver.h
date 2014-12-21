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

#ifndef __SatSolver_h__
#define __SatSolver_h__

#include "Solver.h"
#include "simp/SimpSolver.h"

namespace aspino {

class SatSolver : public Solver, protected Minisat::SimpSolver {
public:
    SatSolver() : inVars(0), status(l_Undef) {}
    virtual ~SatSolver() {}
    
    virtual void interrupt() { Minisat::SimpSolver::interrupt(); }
    
    virtual void parse(gzFile in);
    virtual bool simplify() { return Minisat::SimpSolver::simplify(); }
    virtual bool eliminate(bool turn_off_elim) { return Minisat::SimpSolver::eliminate(turn_off_elim); }
    virtual bool okay() const { return Minisat::SimpSolver::okay(); }
    virtual lbool solve();
    virtual lbool solve(int n);
    void printStatus() const;
    void copyModel();
    void printModel() const;
    
    inline void nInVars(int value) { inVars = value; }
    inline int nInVars() const { return inVars; }
    
    inline int nextToPropagateByUnit() const { return qhead; }
    inline Lit mainTrail(int idx) const { return trail[idx]; }
    inline int mainTrailLimit(int level) const { return trail_lim[level]; }

    virtual inline void newVar() { Minisat::SimpSolver::newVar(); }
    using Minisat::SimpSolver::nVars;
    using Minisat::SimpSolver::addClause;

protected:
    virtual inline CRef morePropagate() { return CRef_Undef; }
    virtual inline void onCancel() {}
    
    void solve_();
    lbool search(int nof_conflicts);
    inline void cancelUntil(int level) { Minisat::SimpSolver::cancelUntil(level); onCancel(); }
    static double luby(double y, int x);
    void learnClauseFromModel();
    
    int inVars;
    int inClauses;
    lbool status;
};

} // namespace miniasp

#endif
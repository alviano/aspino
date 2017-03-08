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

#ifndef __QBFSolver_h__
#define __QBFSolver_h__

#include "PseudoBooleanSolver.h"
#include "DepGraph.h"

namespace aspino {

class QBFSolver : public PseudoBooleanSolver {
public:
    QBFSolver();
    virtual ~QBFSolver();
    
    virtual void interrupt();
    
    virtual void parse(gzFile in);
    
    virtual bool eliminate(bool) { return true; } 
        
    virtual lbool solve();
    virtual lbool solve(int) { return solve(); }
    
    Var getFirstCheckVar() const { return nInVars(); }
    Var getClauseVar(int index) const { return nInVars() + index + 1; }
    Var getClausePVar(int index) const { return nInVars() + clauses.size() + index + 1; }
    
private:
    PseudoBooleanSolver outSolver;
    PseudoBooleanSolver inSolver;
    vec<bool> universal;
    vec<Lit> selectors;
    UndirectedWeightedGraph g;
};


} // namespace aspino

#endif
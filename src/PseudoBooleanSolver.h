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

#ifndef __PseudoBooleanSolver_h__
#define __PseudoBooleanSolver_h__

#include "SatSolver.h"

namespace aspino {

class WeightConstraint {
    friend ostream& operator<<(ostream& out, const WeightConstraint& wc) { return out << wc.toString(); }
public:
    inline WeightConstraint() : bound(0) {}
    inline WeightConstraint(WeightConstraint& wc) : bound(wc.bound), loosable(wc.loosable) { wc.lits.moveTo(lits); wc.coeffs.moveTo(coeffs); }
    ~WeightConstraint() { lits.clear(); coeffs.clear(); trail.clear(); }
    
//    void free(ClauseAllocator& ca);
    
    inline void clear() { lits.clear(); coeffs.clear(); bound = 0; }
    string toString() const;
    
    inline int size() const { assert(lits.size() == coeffs.size()); return lits.size(); }
    inline void shrink_(int n) { lits.shrink_(n); coeffs.shrink_(n); }
    inline void pop() { lits.pop(); coeffs.pop(); }
    
    vec<Lit> lits;
    vec<int> coeffs;
    int bound;
    int loosable;
    vec<int> trail;
    int first;
};
    
class PseudoBooleanSolver : public SatSolver {
public:
    PseudoBooleanSolver() : nextToPropagate(0) {}
    virtual ~PseudoBooleanSolver();
    
    virtual void newVar();
    
    virtual void parse(gzFile in);

    bool addConstraint(WeightConstraint& wc);
    bool addEquality(WeightConstraint& wc);
    void attach(WeightConstraint& wc);
    
protected:
    virtual CRef morePropagate();
    CRef morePropagate(Lit lit);
    void restore(WeightConstraint& wc);
    CRef checkConflict(Lit lit, WeightConstraint& wc, int pos);
    CRef checkInference(Lit lit, WeightConstraint& wc, int pos);
    
    virtual void moreReason(Lit lit);
    virtual void onCancel();
    virtual void relocAll(ClauseAllocator& to);
    
    vec<WeightConstraint*> wconstraints;
    vec< vec<WeightConstraint*> > propagators[2];
    vec< vec<int> > positions[2];

    vec<WeightConstraint*> moreReasonWC;
    vec<int> moreReasonTrailSize;
    vec<Var> moreReasonVars;
    vec<CRef> moreReasonClauses;
    
    vec<bool> propagated;
    int nextToPropagate;
    
private:
    bool readConstraint(Minisat::StreamBuffer& in, WeightConstraint& wc);
    static int gcd(int u, int v);
};

} // namespace aspino

#endif

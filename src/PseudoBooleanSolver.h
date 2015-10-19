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
    
    inline void clear() { lits.clear(); coeffs.clear(); bound = 0; }
    string toString() const;
    
    inline int size() const { assert(lits.size() == coeffs.size()); return lits.size(); }
    inline void shrink_(int n) { lits.shrink_(n); coeffs.shrink_(n); }
    inline void pop() { lits.pop(); coeffs.pop(); }
    
    vec<Lit> lits;
    vec<int64_t> coeffs;
    int bound;
    int loosable;
    vec<int> trail;
    int first;
};

class CardinalityConstraint {
    friend ostream& operator<<(ostream& out, const CardinalityConstraint& cc) { return out << cc.toString(); }
public:
    inline CardinalityConstraint() : bound(0) {}
    inline CardinalityConstraint(CardinalityConstraint& cc) : bound(cc.bound), loosable(cc.loosable) { cc.lits.moveTo(lits); }
    ~CardinalityConstraint() { lits.clear(); trail.clear(); }
    
    inline void clear() { lits.clear(); bound = 0; }
    string toString() const;
    
    inline int size() const { return lits.size(); }
    inline void shrink_(int n) { lits.shrink_(n); }
    inline void pop() { lits.pop(); }
    
    vec<Lit> lits;
    int bound;
    int loosable;
    vec<Lit> trail;
};
    
class PseudoBooleanSolver : public SatSolver {
public:
    static const CRef CRef_MoreConflict;
    
    PseudoBooleanSolver() : moreConflictWC(NULL), moreConflictCC(NULL), nextToPropagate(0) {}
    virtual ~PseudoBooleanSolver();
    
    virtual void newVar();
    
    virtual void parse(gzFile in);

    bool addConstraint(WeightConstraint& wc);
    bool addEquality(WeightConstraint& wc);
    void attach(WeightConstraint& wc);

    bool addConstraint(CardinalityConstraint& wc);
    bool addEquality(CardinalityConstraint& wc);
    void attach(CardinalityConstraint& wc);
    
    void freeze();
    
protected:
    virtual CRef morePropagate();
    CRef morePropagate(Lit lit);
    void restore(WeightConstraint& wc);
    void restore(CardinalityConstraint& wc);
    CRef checkConflict(Lit lit, WeightConstraint& wc, int pos);
    CRef checkConflict(Lit lit, CardinalityConstraint& wc);
    CRef checkInference(Lit lit, WeightConstraint& wc, int pos);
    CRef checkInference(Lit lit, CardinalityConstraint& wc);
    
    virtual bool moreReason(Lit lit, vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC);
    void _moreReasonCC(Lit lit, vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC);
    void _moreReasonWC(Lit lit, vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC);
    virtual bool moreReason(Lit lit);
    void _moreReasonCC(Lit lit);
    void _moreReasonWC(Lit lit);
    virtual bool moreConflict(vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC);
    void _moreConflictCC(vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC);
    void _moreConflictWC(vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC);
    virtual void onCancel();
    
    vec<CardinalityConstraint*> cconstraints;
    vec< vec<CardinalityConstraint*> > cpropagators[2];
    
    vec<WeightConstraint*> wconstraints;
    vec< vec<WeightConstraint*> > wpropagators[2];
    vec< vec<int> > positions[2];

    vec<WeightConstraint*> moreReasonWC;
    vec<CardinalityConstraint*> moreReasonCC;
    vec<int> moreReasonTrailSize;
    vec<Var> moreReasonVars;
    WeightConstraint* moreConflictWC;
    CardinalityConstraint* moreConflictCC;
    Lit moreConflictLit;
    
    vec<bool> propagated;
    int nextToPropagate;
    
private:
    bool readConstraint(Glucose::StreamBuffer& in, WeightConstraint& wc);
    static int gcd(int u, int v);
};

} // namespace aspino

#endif

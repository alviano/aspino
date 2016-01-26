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

#include "PseudoBooleanSolver.h"

#include "utils/ParseUtils.h"
#include "core/Dimacs.h"

#include <iostream>
#include <cmath>
#include <sstream>

using namespace std;

namespace aspino {

const CRef PseudoBooleanSolver::CRef_MoreConflict = CRef_Undef - 1;
 

template<class B>
int64_t parseLong(B& in) {
    int64_t    val = 0;
    bool    neg = false;
    skipWhitespace(in);
    if      (*in == '-') neg = true, ++in;
    else if (*in == '+') ++in;
    if (*in < '0' || *in > '9') fprintf(stderr, "PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    while (*in >= '0' && *in <= '9')
        val = val*10 + (*in - '0'),
        ++in;
    return neg ? -val : val; }
   
//void WeightConstraint::free(Glucose::ClauseAllocator& ca) {
//    for(int i = 0; i < reasons.size(); i++)
//        if(reasons[i] != CRef_Undef)
//            ca.free(reasons[i]);
//}
    
string WeightConstraint::toString() const {
    stringstream ss;
    if(lits.size() == 0)
        ss << "0 ";
    else {
        for(int i = 0; i < lits.size(); i++)
            ss << (coeffs[i] >= 0 ? "+" : "") << coeffs[i] << "x" << lits[i] << " ";
    }
    ss << ">= " << bound;
    return ss.str();
}

string CardinalityConstraint::toString() const {
    stringstream ss;
    if(lits.size() == 0)
        ss << "0 ";
    else {
        for(int i = 0; i < lits.size(); i++)
            ss << "+x" << lits[i] << " ";
    }
    ss << ">= " << bound;
    return ss.str();
}

PseudoBooleanSolver::~PseudoBooleanSolver() {
    cpropagators[0].clear();
    cpropagators[1].clear();
    wpropagators[0].clear();
    wpropagators[1].clear();
    
    for(int i = 0; i < cconstraints.size(); i++) delete  cconstraints[i];
    cconstraints.clear();

    for(int i = 0; i < wconstraints.size(); i++) delete  wconstraints[i];
    wconstraints.clear();
}

bool PseudoBooleanSolver::readConstraint(Glucose::StreamBuffer& in, WeightConstraint& wc) {
    int var;
    bool ret = false;
    wc.clear();
    for (;;){
        if(*in == '>' || *in == '=') {
            if(*in == '=') ret = true;
            else {
                ++in;
                if(*in != '=') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            }
            ++in;
            skipWhitespace(in);
            wc.bound = parseLong(in);
            skipWhitespace(in);
            if(*in != ';') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            ++in;
            return ret;
        }
        
        wc.coeffs.push(parseLong(in));
        skipWhitespace(in);
        if(*in != 'x') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
        ++in;
        var = parseInt(in);
        if(var-- <= 0) cerr << "PARSE ERROR! Variable with non-positive index: x" << var+1 << endl, exit(3);
        while(var >= nVars()) newVar();
        wc.lits.push(mkLit(var));
        skipWhitespace(in);
    }
}
    
void PseudoBooleanSolver::parse(gzFile in_) {
    Glucose::StreamBuffer in(in_);

    WeightConstraint wc;
    for(;;) {
        skipWhitespace(in);
        if(*in == EOF) break;
        if(*in == '*')
            skipLine(in);
        else {
            if(readConstraint(in, wc))
                addEquality(wc);
            else
                addConstraint(wc);
            morePropagate();
        }
    }
    inVars = nVars();
    freeze();
}

void PseudoBooleanSolver::freeze() {
    for(int i = 0; i < cconstraints.size(); i++)
        for(int j = 0; j < cconstraints[i]->size(); j++)
            setFrozen(var(cconstraints[i]->lits[j]), true);
    for(int i = 0; i < wconstraints.size(); i++)
        for(int j = 0; j < wconstraints[i]->size(); j++)
            setFrozen(var(wconstraints[i]->lits[j]), true);
}

bool PseudoBooleanSolver::addConstraint(WeightConstraint& wc) {
    assert(decisionLevel() == 0);
    assert(wc.size() > 0);

    trace(pbs, 2, "Processing constraint: " << wc);
    
    // remove literals with zero coefficient and literals at level zero, and transform negative coefficients
    int j = 0;
    for(int i = 0; i < wc.size(); i++) {
        if(wc.coeffs[i] == 0) continue;
        if(value(wc.lits[i]) != l_Undef) {
            assert(level(var(wc.lits[i])) == 0);
            if(value(wc.lits[i]) == l_True)
                wc.bound -= wc.coeffs[i];
            continue;
        }
        if(wc.coeffs[i] < 0) {
            wc.coeffs[j] = -wc.coeffs[i];
            wc.bound += wc.coeffs[j];
            wc.lits[j++] = ~wc.lits[i];
        }
        else {
            wc.coeffs[j] = wc.coeffs[i];
            wc.lits[j++] = wc.lits[i];
        }
    }
    wc.shrink_(wc.size()-j);
    trace(pbs, 10, "After removing zero coefficients and literals at level zero: " << wc);
    
    // trivially satisfied
    if(wc.bound <= 0) {
        trace(pbs, 4, "Trivially satisfied constraint: " << wc);
        return true;
    }
    
    // flattening coefficients greater than bound and compute gcd and sum
    int64_t d = wc.bound;
    int64_t sum = 0;
    j = 0;
    for(int i = 0; i < wc.size(); i++) {
        wc.coeffs[j] = wc.coeffs[i] > wc.bound ? wc.bound : wc.coeffs[i];
        d = gcd(wc.coeffs[j], d);
        sum += wc.coeffs[j];
        wc.lits[j++] = wc.lits[i];
    }
    trace(pbs, 10, "After flattening bounds: " << wc);
    
    // cannot be satisfied
    if(sum < wc.bound) {
        trace(pbs, 4, "Unsatisfiable constraint: " << wc);
        return ok = false;
    }
    
    if(d > 1) {
        wc.bound /= d;
        sum /= d;
        for(int i = 0; i < wc.size(); i++)
            wc.coeffs[i] /= d;
    }
    trace(pbs, 10, "After dividing bounds for " << d << ": " << wc);
    
    // sort
    int n = wc.size();
    while(n > 0) {
       int newn = 0;
       for(int i = 1; i < n; i++) {
          if(wc.coeffs[i-1] > wc.coeffs[i]) {
             int64_t ctmp = wc.coeffs[i-1];
             wc.coeffs[i-1] = wc.coeffs[i];
             wc.coeffs[i] = ctmp;
             
             Lit ltmp = wc.lits[i-1];
             wc.lits[i-1] = wc.lits[i];
             wc.lits[i] = ltmp;
             
             newn = i;
          }
       }
       n = newn;
    }
    trace(pbs, 10, "After sorting bounds: " << wc);
    
    // find required literals
    wc.loosable = sum - wc.bound;
    assert(wc.loosable >= 0);
    while(wc.size() > 0) {
        if(wc.coeffs.last() <= wc.loosable) break;
        trace(pbs, 4, "Required literal " << wc.lits.last() << " in constraint " << wc);
        if(!addClause(wc.lits.last())) return ok = false;
        wc.bound -= wc.coeffs.last();
        wc.pop();
    }
    if(wc.bound <= 0) {
        trace(pbs, 4, "Trivially satisfied constraint: " << wc);
        return true;
    }
    assert(wc.bound > 0);
    //if(lits.size() == 0) return true;
    assert(wc.size() > 1);
    trace(pbs, 10, "After removing required literals: " << wc);
    
    // check if it is an at-least-one
    if(wc.bound <= wc.coeffs[0]) {
        trace(pbs, 4, "At-least-one constraint: " << wc);
        return addClause(wc.lits);
    }
    
    attach(wc);
    
    return true;
}

bool PseudoBooleanSolver::addEquality(WeightConstraint& wc) {
    WeightConstraint wc2;
    wc2.lits.growTo(wc.size());
    wc2.coeffs.growTo(wc.size());
    for(int i = 0; i < wc.size(); i++) {
        wc2.lits[i] = wc.lits[i];
        wc2.coeffs[i] = -wc.coeffs[i];
    }
    wc2.bound = -wc.bound;

    return addConstraint(wc) && addConstraint(wc2);
}

void PseudoBooleanSolver::attach(WeightConstraint& wc) {
    trace(pbs, 4, "Adding constraint: " << wc);

    WeightConstraint* ref = new WeightConstraint(wc);
    ref->first = ref->size() - 1;
    wconstraints.push(ref);
    
    for(int i = 0; i < ref->size(); i++) {
        wpropagators[1-sign(ref->lits[i])][var(ref->lits[i])].push(ref);
        positions[1-sign(ref->lits[i])][var(ref->lits[i])].push(i);
    }
}

bool PseudoBooleanSolver::addConstraint(CardinalityConstraint& cc) {
    assert(decisionLevel() == 0);
    assert(cc.size() > 0);

    trace(pbs, 2, "Processing constraint: " << cc);
    
    // remove literals at level zero, and transform negative coefficients
    int j = 0;
    for(int i = 0; i < cc.size(); i++) {
        if(value(cc.lits[i]) != l_Undef) {
            assert(level(var(cc.lits[i])) == 0);
            if(value(cc.lits[i]) == l_True)
                cc.bound--;
            continue;
        }
        cc.lits[j++] = cc.lits[i];
    }
    cc.shrink_(cc.size()-j);
    trace(pbs, 10, "After removing literals at level zero: " << cc);
    
    // trivially satisfied
    if(cc.bound <= 0) {
        trace(pbs, 4, "Trivially satisfied constraint: " << cc);
        return true;
    }
    
    // cannot be satisfied
    if(cc.size() < cc.bound) {
        trace(pbs, 4, "Unsatisfiable constraint: " << cc);
        return ok = false;
    }
    
    // find required literals
    cc.loosable = cc.size() - cc.bound;
    assert(cc.loosable >= 0);
    while(cc.size() > 0) {
        if(1 <= cc.loosable) break;
        trace(pbs, 4, "Required literal " << cc.lits.last() << " in constraint " << cc);
        if(!addClause(cc.lits.last())) return ok = false;
        cc.bound--;
        cc.pop();
    }
    if(cc.bound <= 0) {
        trace(pbs, 4, "Trivially satisfied constraint: " << cc);
        return true;
    }
    assert(cc.bound > 0);
    //if(lits.size() == 0) return true;
    assert(cc.size() > 1);
    trace(pbs, 10, "After removing required literals: " << cc);
    
    // check if it is an at-least-one
    if(cc.bound <= 1) {
        trace(pbs, 4, "At-least-one constraint: " << cc);
        return addClause(cc.lits);
    }
    
    attach(cc);
    
    return true;
}

bool PseudoBooleanSolver::addEquality(CardinalityConstraint& cc) {
    CardinalityConstraint cc2;
    cc2.lits.growTo(cc.size());
    for(int i = 0; i < cc.size(); i++) cc2.lits[i] = ~cc.lits[i];
    cc2.bound = -cc.bound + cc2.size();
    return addConstraint(cc) && addConstraint(cc2);
}

void PseudoBooleanSolver::attach(CardinalityConstraint& cc) {
    trace(pbs, 4, "Adding constraint: " << cc);

    CardinalityConstraint* ref = new CardinalityConstraint(cc);
    cconstraints.push(ref);
    
    for(int i = 0; i < ref->size(); i++) cpropagators[1-sign(ref->lits[i])][var(ref->lits[i])].push(ref);
}

CRef PseudoBooleanSolver::morePropagate() {
    while(nextToPropagate < nextToPropagateByUnit()) {
        Lit lit = mainTrail(nextToPropagate++);
        CRef ret = morePropagate(lit);
        if(ret != CRef_Undef) return ret;
    }
        
    return CRef_Undef;
}

CRef PseudoBooleanSolver::morePropagate(Lit lit) {
    trace(pbs, 20, "Propagating " << lit << "@" << level(var(lit)));

    vec<CardinalityConstraint*>& c = cpropagators[sign(lit)][var(lit)];
    for(int i = 0; i < c.size(); i++) {
        CardinalityConstraint& cc = *c[i];

        trace(pbs, 6, "Processing " << cc);
        restore(cc);
        CRef ret = checkConflict(lit, cc);
        if(ret != CRef_Undef) return ret;
        ret = checkInference(lit, cc);
        if(ret != CRef_Undef) return ret;
    }

    vec<WeightConstraint*>& p = wpropagators[sign(lit)][var(lit)];
    vec<int>& ppos = positions[sign(lit)][var(lit)];
    for(int i = 0; i < p.size(); i++) {
        WeightConstraint& wc = *p[i];
        int pos = ppos[i];

        trace(pbs, 6, "Processing " << wc);
        restore(wc);
        CRef ret = checkConflict(lit, wc, pos);
        if(ret != CRef_Undef) return ret;
        ret = checkInference(lit, wc, pos);
        if(ret != CRef_Undef) return ret;
    }
        
    propagated[var(lit)] = true;
    return CRef_Undef;
}

void PseudoBooleanSolver::restore(WeightConstraint& wc) {
    trace(pbs, 10, "Restoring status of " << wc << " (first: " << wc.first << "; loosable: " << wc.loosable << "; trail: " << wc.trail << ")");
    for(;;) {
        if(wc.trail.size() == 0) break;
        int idx = wc.trail.last();
        if(propagated[var(wc.lits[idx])] && value(wc.lits[idx]) == l_False) break;
        trace(pbs, 15, "Removing literal " << wc.lits[idx] << " from the trail: loosable was " << wc.loosable << " and now is " << wc.loosable + wc.coeffs[idx]);
        wc.trail.pop();
        wc.loosable += wc.coeffs[idx];
        if(wc.first < idx) {
            trace(pbs, 17, "First was " << wc.first << " and now is " << idx);
            wc.first = idx;
        }
    }
}

void PseudoBooleanSolver::restore(CardinalityConstraint& cc) {
    trace(pbs, 10, "Restoring status of " << cc << " (loosable: " << cc.loosable << "; trail: " << cc.trail << ")");
    for(;;) {
        if(cc.trail.size() == 0) break;
        Lit last = cc.trail.last();
        if(propagated[var(last)] && value(last) == l_False) break;
        trace(pbs, 15, "Removing literal " << last << " from the trail: loosable was " << cc.loosable << " and now is " << cc.loosable + 1);
        cc.trail.pop();
        cc.loosable++;
    }
}

CRef PseudoBooleanSolver::checkConflict(Lit lit, WeightConstraint& wc, int pos) {
    if(wc.loosable < wc.coeffs[pos]) {
        trace(pbs, 4, "Conflict on " << wc);
        moreConflictWC = &wc;
        moreConflictLit = ~lit;
        return CRef_MoreConflict;
    }
    return CRef_Undef;
}

CRef PseudoBooleanSolver::checkConflict(Lit lit, CardinalityConstraint& cc) {
    if(cc.loosable < 1) {
        trace(pbs, 4, "Conflict on " << cc);
        moreConflictCC = &cc;
        moreConflictLit = ~lit;
        return CRef_MoreConflict;
    }
    return CRef_Undef;
}

CRef PseudoBooleanSolver::checkInference(Lit lit, WeightConstraint& wc, int pos) {
    trace(pbs, 10, "Adding literal " << wc.lits[pos] << " to the trail of " << wc << ": loosable was " << wc.loosable << " and now is " << wc.loosable - wc.coeffs[pos]);
    wc.trail.push(pos);
    wc.loosable -= wc.coeffs[pos];
    for(int j = wc.first; j >= 0; j--) {
        if(wc.coeffs[j] <= wc.loosable) break;
        Lit wlit = wc.lits[j];
        if(value(wlit) == l_True) continue;
        if(value(wlit) == l_Undef) {
            trace(pbs, 20, "Inferring " << wlit << "@" << decisionLevel() << " from " << wc);
            assert(moreReasonWC[var(wlit)] == NULL);
            moreReasonWC[var(wlit)] = &wc;
            moreReasonTrailSize[var(wlit)] = wc.trail.size();
            moreReasonVars.push(var(wlit));
            uncheckedEnqueue(wlit);
        }
        else if(wlit != ~lit && !propagated[var(wlit)]) {
            trace(pbs, 4, "Conflict on literal " << wlit << " in " << wc);
            moreConflictWC = &wc;
            moreConflictLit = wlit;
            return CRef_MoreConflict;
        }
    }
    while(value(wc.lits[wc.first]) == l_False) if(--wc.first < 0) break;
    
    return CRef_Undef;
}

CRef PseudoBooleanSolver::checkInference(Lit lit, CardinalityConstraint& cc) {
    trace(pbs, 10, "Adding literal " << ~lit << " to the trail of " << cc << ": loosable was " << cc.loosable << " and now is " << cc.loosable - 1);
    cc.trail.push(~lit);
    cc.loosable--;
    assert(cc.loosable >= 0);
    
    if(cc.loosable == 0) {
        for(int i = 0; i < cc.size(); i++) {
            Lit clit = cc.lits[i];
            if(value(clit) == l_True) continue;
            if(value(clit) == l_Undef) {
                trace(pbs, 20, "Inferring " << clit << "@" << decisionLevel() << " from " << cc);
                assert(moreReasonCC[var(clit)] == NULL);
                moreReasonCC[var(clit)] = &cc;
                moreReasonTrailSize[var(clit)] = cc.trail.size();
                moreReasonVars.push(var(clit));
                uncheckedEnqueue(clit);
            }
            else if(clit != ~lit && !propagated[var(clit)]) {
                trace(pbs, 4, "Conflict on literal " << clit << " in " << cc);
                moreConflictCC = &cc;
                moreConflictLit = clit;
                return CRef_MoreConflict;
            }
        }
    }
    
    return CRef_Undef;
}

void PseudoBooleanSolver::newVar() {
    SatSolver::newVar();
    cpropagators[0].push();
    cpropagators[1].push();
    wpropagators[0].push();
    wpropagators[1].push();
    positions[0].push();
    positions[1].push();
    moreReasonWC.push(NULL);
    moreReasonCC.push(NULL);
    moreReasonTrailSize.push();
    propagated.push(false);
}

bool PseudoBooleanSolver::moreReason(Lit lit, vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC) {
    if(moreReasonCC[var(lit)] != NULL) { _moreReasonCC(lit, out_learnt, selectors, pathC); return true; }
    if(moreReasonWC[var(lit)] != NULL) { _moreReasonWC(lit, out_learnt, selectors, pathC); return true; }
    return false;
}

void PseudoBooleanSolver::_moreReasonCC(Lit lit, vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC) {
    assert(decisionLevel() != 0);
    assert(reason(var(lit)) == CRef_Undef);
    CardinalityConstraint& cc = *moreReasonCC[var(lit)];
    int trailSize = moreReasonTrailSize[var(lit)];
    assert(trailSize >= 0);
    assert(trailSize <= cc.trail.size());
    while(--trailSize >= 0) {
        Lit q = cc.trail[trailSize];
        assert(value(q) == l_False);
        assert(level(var(q)) <= level(var(lit)));
        
        if(seen[var(q)]) continue;
        if(level(var(q)) == 0) continue;
        
        if(!isSelector(var(q)))
            varBumpActivity(var(q));
        
        seen[var(q)] = 1;
        
        if(level(var(q)) >= decisionLevel()) {
            pathC++;
            // UPDATEVARACTIVITY trick (see competition'09 companion paper)
            if(!isSelector(var(q)) && (reason(var(q)) != CRef_Undef) && ca[reason(var(q))].learnt())
                lastDecisionLevel.push(q);
        }
        else {
            if(isSelector(var(q))) {
                assert(value(q) == l_False);
                selectors.push(q);
            }
            else 
                out_learnt.push(q);
        }
    }
}

void PseudoBooleanSolver::_moreReasonWC(Lit lit, vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC) {
    assert(decisionLevel() != 0);
    assert(reason(var(lit)) == CRef_Undef);
    WeightConstraint& wc = *moreReasonWC[var(lit)];
    int trailSize = moreReasonTrailSize[var(lit)];
    assert(trailSize >= 0);
    assert(trailSize <= wc.trail.size());
    while(--trailSize >= 0) {
        Lit q = wc.lits[wc.trail[trailSize]];
        assert(value(q) == l_False);
        assert(level(var(q)) <= level(var(lit)));
        
        if(seen[var(q)]) continue;
        if(level(var(q)) == 0) continue;
        
        if(!isSelector(var(q)))
            varBumpActivity(var(q));
        
        seen[var(q)] = 1;
        
        if(level(var(q)) >= decisionLevel()) {
            pathC++;
            // UPDATEVARACTIVITY trick (see competition'09 companion paper)
            if(!isSelector(var(q)) && (reason(var(q)) != CRef_Undef) && ca[reason(var(q))].learnt())
                lastDecisionLevel.push(q);
        }
        else {
            if(isSelector(var(q))) {
                assert(value(q) == l_False);
                selectors.push(q);
            }
            else 
                out_learnt.push(q);
        }
    }
}

bool PseudoBooleanSolver::moreReason(Lit lit) {
    if(moreReasonCC[var(lit)] != NULL) { _moreReasonCC(lit); return true; }
    if(moreReasonWC[var(lit)] != NULL) { _moreReasonWC(lit); return true; }
    return false;
}

void PseudoBooleanSolver::_moreReasonCC(Lit lit) {
    assert(decisionLevel() != 0);
    assert(reason(var(lit)) == CRef_Undef);
    CardinalityConstraint& cc = *moreReasonCC[var(lit)];
    int trailSize = moreReasonTrailSize[var(lit)];
    assert(trailSize >= 0);
    assert(trailSize <= cc.trail.size());
    while(--trailSize >= 0) {
        Lit l = cc.trail[trailSize];
        assert(value(l) == l_False);
        assert(level(var(l)) <= level(var(lit)));
        if(level(var(l)) == 0) continue;
        seen[var(l)] = 1;
    }
}

void PseudoBooleanSolver::_moreReasonWC(Lit lit) {
    assert(decisionLevel() != 0);
    assert(reason(var(lit)) == CRef_Undef);
    WeightConstraint& wc = *moreReasonWC[var(lit)];
    int trailSize = moreReasonTrailSize[var(lit)];
    assert(trailSize >= 0);
    assert(trailSize <= wc.trail.size());
    while(--trailSize >= 0) {
        Lit l = wc.lits[wc.trail[trailSize]];
        assert(value(l) == l_False);
        assert(level(var(l)) <= level(var(lit)));
        if(level(var(l)) == 0) continue;
        seen[var(l)] = 1;
    }
}

bool PseudoBooleanSolver::moreConflict(vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC) {
    if(moreConflictCC != NULL) { _moreConflictCC(out_learnt, selectors, pathC); return true; }
    if(moreConflictWC != NULL) { _moreConflictWC(out_learnt, selectors, pathC); return true; }
    return false;
}
    
void PseudoBooleanSolver::_moreConflictCC(vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC) {
    assert(decisionLevel() != 0);
    
    if(!seen[var(moreConflictLit)] && level(var(moreConflictLit)) > 0) {
        if(!isSelector(var(moreConflictLit))) varBumpActivity(var(moreConflictLit));
        seen[var(moreConflictLit)] = 1;
        assert(level(var(moreConflictLit)) == decisionLevel());
        pathC++;
        // UPDATEVARACTIVITY trick (see competition'09 companion paper)
        if(!isSelector(var(moreConflictLit)) && (reason(var(moreConflictLit)) != CRef_Undef) && ca[reason(var(moreConflictLit))].learnt())
            lastDecisionLevel.push(moreConflictLit);
    }
    
    CardinalityConstraint& cc = *moreConflictCC;
    moreConflictCC = NULL;
    int trailSize = cc.trail.size();
    assert(trailSize >= 0);
    while(--trailSize >= 0) {
        Lit q = cc.trail[trailSize];
        assert(value(q) == l_False);
        
        if(seen[var(q)]) continue;
        if(level(var(q)) == 0) continue;
        
        if(!isSelector(var(q)))
            varBumpActivity(var(q));
        
        seen[var(q)] = 1;
        
        if(level(var(q)) >= decisionLevel()) {
            pathC++;
            // UPDATEVARACTIVITY trick (see competition'09 companion paper)
            if(!isSelector(var(q)) && (reason(var(q)) != CRef_Undef) && ca[reason(var(q))].learnt())
                lastDecisionLevel.push(q);
        }
        else {
            if(isSelector(var(q))) {
                assert(value(q) == l_False);
                selectors.push(q);
            }
            else 
                out_learnt.push(q);
        }
    }
}

void PseudoBooleanSolver::_moreConflictWC(vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC) {
    assert(decisionLevel() != 0);
    
    if(!seen[var(moreConflictLit)] && level(var(moreConflictLit)) > 0) {
        if(!isSelector(var(moreConflictLit))) varBumpActivity(var(moreConflictLit));
        seen[var(moreConflictLit)] = 1;
        assert(level(var(moreConflictLit)) == decisionLevel());
        pathC++;
        // UPDATEVARACTIVITY trick (see competition'09 companion paper)
        if(!isSelector(var(moreConflictLit)) && (reason(var(moreConflictLit)) != CRef_Undef) && ca[reason(var(moreConflictLit))].learnt())
            lastDecisionLevel.push(moreConflictLit);
    }
    
    WeightConstraint& wc = *moreConflictWC;
    moreConflictWC = NULL;
    int trailSize = wc.trail.size();
    assert(trailSize >= 0);
    while(--trailSize >= 0) {
        Lit q = wc.lits[wc.trail[trailSize]];
        assert(value(q) == l_False);
        
        if(seen[var(q)]) continue;
        if(level(var(q)) == 0) continue;
        
        if(!isSelector(var(q)))
            varBumpActivity(var(q));
        
        seen[var(q)] = 1;
        
        if(level(var(q)) >= decisionLevel()) {
            pathC++;
            // UPDATEVARACTIVITY trick (see competition'09 companion paper)
            if(!isSelector(var(q)) && (reason(var(q)) != CRef_Undef) && ca[reason(var(q))].learnt())
                lastDecisionLevel.push(q);
        }
        else {
            if(isSelector(var(q))) {
                assert(value(q) == l_False);
                selectors.push(q);
            }
            else 
                out_learnt.push(q);
        }
    }
}

void PseudoBooleanSolver::onCancel() {
    trace(pbs, 2, "Cancel until level " << decisionLevel());
    while(nextToPropagate > nextToPropagateByUnit()) { 
        int v = var(mainTrail(--nextToPropagate));
        propagated[v] = false;
    }

    while(moreReasonVars.size() > 0) {
        Var v = moreReasonVars.last();
        if(level(v) <= decisionLevel()) break;
        moreReasonVars.pop();
        moreReasonWC[v] = NULL;
        moreReasonCC[v] = NULL;
        moreReasonTrailSize[v] = -1;
    }
}

int64_t PseudoBooleanSolver::gcd(int64_t a, int64_t b) {
    assert(a > 0);
    assert(b > 0);
    
    if(a <= b)
        b = b % a;

    int64_t tmp;
    while(b > 0) {
        tmp = a;
        a = b;
        b = tmp % b;
    }
    return a;
}

} // namespace aspino

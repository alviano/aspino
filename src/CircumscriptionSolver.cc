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

#include "CircumscriptionSolver.h"

#include "utils/algorithm.h"
#include <core/Dimacs.h>

extern Glucose::BoolOption option_print_model;

Glucose::IntOption option_circumscription_witnessess("CIRCUMSCRIPTION", "circ-wit", "Number of desired witnesses for each solution. Non-positive integers are interpreted as unbounded.\n", 1, Glucose::IntRange(0, INT32_MAX));

namespace aspino {

CircumscriptionSolver::CircumscriptionSolver() : PseudoBooleanSolver() {
    setIncrementalMode();
}

CircumscriptionSolver::~CircumscriptionSolver() {
    for(int i = 0; i < visibleValue.size(); i++) delete visibleValue[i];
    visibleValue.clear();
}

void CircumscriptionSolver::parse(gzFile in_) {
    cerr << "c parsing start" << endl;
    Glucose::StreamBuffer in(in_);

    char buff[1024];
    vec<Lit> lits;
    bool min = true;
    for(;;) {
        skipWhitespace(in);
        if(*in == EOF) break;
        if(*in == 'p') {
            if(!eagerMatch(in, "p ccnf")) cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            skipWhitespace(in);
            if(*in == '+') min = false;
            else if(*in != '-') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            skipLine(in);
        }
        else if(*in == 'o') {
            if(objectLits.size() != 0) cerr << "PARSE ERROR! Only one line can start with o: " << static_cast<char>(*in) << endl, exit(3);
            ++in;
            readClause(in, *this, objectLits);
            if(!checkConsistencyAndRemoveDuplicates(objectLits)) cerr << "PARSE ERROR! Inconsistency in o literals!" << endl, exit(3);
            for(int i = 0; i < nOLits(); i++) setFrozen(var(getOLit(i)), true);
        }
        else if(*in == 'C') {
            if(careLits.size() != 0) cerr << "PARSE ERROR! Only one line can start with C: " << static_cast<char>(*in) << endl, exit(3);
            ++in;
            readClause(in, *this, careLits);
            if(!checkConsistencyAndRemoveDuplicates(careLits)) cerr << "PARSE ERROR! Inconsistency in C literals!" << endl, exit(3);
            for(int i = 0; i < nCLits(); i++) setFrozen(var(getCLit(i)), true);
        }
        else if(*in == 'v') {
            ++in;
            int parsed_lit = parseInt(in);
            Var var = abs(parsed_lit)-1;
            while (var >= nVars()) newVar();
            Lit lit = (parsed_lit > 0) ? mkLit(var) : ~mkLit(var);
            
            ++in;
            int count = 0;
            while(*in != EOF && *in != '\n' && static_cast<unsigned>(count) < sizeof(buff)) { buff[count++] = *in; ++in; }
            if(static_cast<unsigned>(count) >= sizeof(buff)) cerr << "PARSE ERROR! Value of visible literal " << lit << " is too long!" << endl, exit(3);
            buff[count] = '\0';
            
            visible.push(lit);
            char* value = new char[count+1];
            strcpy(value, buff);
            visibleValue.push(value);
            
            if(option_circumscription_witnessess != 1) setFrozen(var, true);
        }
        else if(*in == '+') {
            ++in;
            CardinalityConstraint cc;
            cc.bound = parseLong(in);
            readClause(in, *this, cc.lits);
            addConstraint(cc);
        }
        else if(*in == 'c')
            skipLine(in);
        else {
            readClause(in, *this, lits);
            addClause(lits); 
        }
    }
    nInVars(nVars());
    
    if(!disjoint(objectLits, careLits)) cerr << "Object literals and care literals must be disjoint!" << endl, exit(3);
    
    if(min) for(int i = 0; i < objectLits.size(); i++) objectLits[i] = ~objectLits[i];

    cerr << "c parsing end" << endl;
}

void CircumscriptionSolver::interrupt() { 
    cout << "s UNKNOWN" << endl;
    this->exit(1);
}

void CircumscriptionSolver::enumerate() {
    cancelUntil(0);
    
    for(int i = 0; i < nonInputSoftLits.size(); i++) assumptions.push(~nonInputSoftLits[i]);
    int nAss = assumptions.size();
    
    vec<Lit> branchingLits;
    vec<bool> branchingFlags;
    
    int witnesses = 0;
    while(true) {
        assumptions.shrink_(assumptions.size() - nAss);
        assert(assumptions.size() == nAss);
        for(int i = 0; i < branchingLits.size(); i++) assumptions.push(branchingLits[i]);
        
        PseudoBooleanSolver::solve();
        
        if(status == l_True) {
            for(int i = assumptions.size(); i < decisionLevel(); i++) {
                branchingLits.push(mainTrail(mainTrailLimit(i)));
                branchingFlags.push(false);
            }
            
            onComputedModel();
            if(computedModels == requiredModels) break;
            if(++witnessess == requiredWitnesses) break;
        }
        else if(status == l_False) {
            // backjump
            Var v = conflict.size() > 0 ? var(conflict[0]) : -1;
            while(branchingLits.size() > 0 && v != var(branchingLits.last())) {
                branchingLits.pop();
                branchingFlags.pop();
            }
        }
    
        // remove flipped branching literals
        while(branchingLits.size() > 0 && branchingFlags.last()) {
            branchingLits.pop();
            branchingFlags.pop();
        }
    
        if(branchingLits.size() == 0) break;

        // flip top branching literal
        branchingLits[branchingLits.size()-1] = ~branchingLits[branchingLits.size()-1];
        branchingFlags[branchingFlags.size()-1] = true;
        
        cancelUntil(branchingLits.size()-1);
    }
}

void CircumscriptionSolver::enumerateByBlockingClauses() {
    cancelUntil(0);
    
    for(int i = 0; i < nonInputSoftLits.size(); i++) assumptions.push(~nonInputSoftLits[i]);
    
    vec<Lit> blockingClause;
    int witnesses = 0;
    while(true) {
        PseudoBooleanSolver::solve();
        
        if(status == l_False) break;
 
        onComputedModel();
        if(computedModels == requiredModels) break;
        if(++witnessess == requiredWitnesses) break;
        
        assert(blockingClause.size() == 0);
        for(int i = 0; i < nCLits(); i++) blockingClause.push(value(getCLit(i)) == l_True ? ~getCLit(i) : getCLit(i));
        for(int i = 0; i < nOLits(); i++) if(value(getOLit(i)) == l_False) blockingClause.push(getOLit(i));
        cancelUntil(0);
        addClause(blockingClause);
        blockingClause.clear();
    }
}

void CircumscriptionSolver::trim() {
    int counter = 0;
    do{
        counter++;
        assumptions.clear();
        for(int i = 0; i < conflict.size(); i++) assumptions.push(~conflict[i]);
        PseudoBooleanSolver::solve();
        assert(status == l_False);
        cancelUntil(0);
        if(conflict.size() <= 1) return;
    }while(assumptions.size() > conflict.size());
    if(counter % 2 == 1) for(int i = 0; i < assumptions.size(); i++) conflict[i] = ~assumptions[i];
    assert(conflict.size() > 1);
}

void CircumscriptionSolver::shrink() {
    assert(decisionLevel() == 0);
    if(conflict.size() <= 1) return;
    
    trim();

    vec<Lit> core;
    conflict.moveTo(core);
    
    vec<Lit> allAssumptions;
    for(int i = 0; i < core.size(); i++) allAssumptions.push(~core[i]);
    
    assumptions.clear();
    const int progressionFrom = 1;
    int progression = progressionFrom;
    int fixed = 0;
    while(true) {
        if(fixed + progression >= allAssumptions.size()) {
            if(progression == progressionFrom) break;
            progression = progressionFrom;
            fixed = assumptions.size();
            continue;
        }

        trace(maxsat, 15, "Minimize: progress to " << progression << "; fixed = " << fixed);
        
        int prec = assumptions.size();
        for(int i = assumptions.size(); i < fixed + progression; i++) {
            assert(i < allAssumptions.size());
            assumptions.push(allAssumptions[i]);
        }
        
        PseudoBooleanSolver::solve();

        if(status == l_False) {
            trace(maxsat, 10, "Minimize: reduce to size " << conflict.size());
            progression = progressionFrom;
            
            assumptions.moveTo(core);
            cancelUntil(0);
            trim();
            core.moveTo(assumptions);
            conflict.moveTo(core);
            
            int j = 0;
            for(int i = 0, k = core.size() - 1; i < prec; i++) {
                if(k < 0) break;
                if(assumptions[i] != ~core[k]) continue;
                assumptions[j++] = assumptions[i];
                k--;
            }
            assumptions.shrink_(assumptions.size() - j);
            fixed = assumptions.size();
            
            j = 0;
            for(int i = 0, k = core.size() - 1; i < allAssumptions.size(); i++) {
                if(k < 0) break;
                if(allAssumptions[i] != ~core[k]) continue;
                allAssumptions[j++] = allAssumptions[i];
                k--;
            }
            allAssumptions.shrink_(allAssumptions.size() - j);
        }
        else {
            trace(maxsat, 20, (status == l_True ? "SAT!" : "UNDEF"));
            progression *= 2;
        }
        cancelUntil(0);
    }
    core.moveTo(conflict);
}

void CircumscriptionSolver::one() {
    assert(decisionLevel() == 0);

    CardinalityConstraint cc;
    cc.bound = conflict.size() - 1;
    while(conflict.size() > 0) {
        assert(!removed[var(conflict.last())]);
        removed[var(conflict.last())] = true;
        cc.lits.push(~conflict.last());
        conflict.pop();
    }
    assert(conflict.size() == 0);
    
    int newLen = 0;
    for(int i = 0; i < softLits.size(); i++) {
        softLits[newLen] = softLits[i];
        if(!removed[var(softLits[i])]) newLen++;
    }
    softLits.shrink(softLits.size()-newLen);
    
    if(cc.bound == 0) return;   // core was singleton
    
    for(int i = 0; i < cc.bound; i++) {
        newVar(); setFrozen(nVars()-1, true);
        removed.push(false);
        
        if(i != 0) addClause(~softLits.last(), mkLit(nVars()-1));
        
        softLits.push(mkLit(nVars()-1));
        nonInputSoftLits.push(softLits.last());
        cc.lits.push(~softLits.last());
    }
    
    assert(cc.size() > 1);
    addConstraint(cc);
}

void CircumscriptionSolver::pmres() {
    assert(decisionLevel() == 0);

    Lit prec = lit_Undef;
    vec<Lit> lits;
    while(conflict.size() > 0) {
        assert(!removed[var(conflict.last())]);
        removed[var(conflict.last())] = true;
        
        if(prec == lit_Undef) prec = ~conflict.last();
        else {
            // disjunction
            newVar(); setFrozen(nVars()-1, true);
            removed.push(false);
            softLits.push(mkLit(nVars()-1));
            nonInputSoftLits.push(softLits.last());
            lits.push(prec);
            lits.push(~conflict.last());
            lits.push(mkLit(nVars()-1, true));
            addClause(lits);
            lits.clear();
            addClause(~prec, mkLit(nVars()-1));
            addClause(conflict.last(), mkLit(nVars()-1));
            
            if(conflict.size() > 1) {
                // conjunction
                newVar();
                removed.push(false);
                addClause(prec, mkLit(nVars()-1, true));
                addClause(~conflict.last(), mkLit(nVars()-1, true));
                lits.push(~prec);
                lits.push(conflict.last());
                lits.push(mkLit(nVars()-1));
                addClause(lits);
                lits.clear();
                
                prec = mkLit(nVars()-1);
            }
        }
        
        conflict.pop();
    }
    assert(conflict.size() == 0);
    
    int newLen = 0;
    for(int i = 0; i < softLits.size(); i++) {
        softLits[newLen] = softLits[i];
        if(!removed[var(softLits[i])]) newLen++;
    }
    softLits.shrink(softLits.size()-newLen);
}

void CircumscriptionSolver::kdyn() {
    assert(decisionLevel() == 0);
    
    const int b = conflict.size() <= 2 ? 8 : ceil(log10(conflict.size()) * 16);
    const int m = ceil(2.0 * conflict.size() / (b-2.0));
    const int N = ceil(
            (
                conflict.size()         // literals in the core
                + conflict.size() - 1   // new soft literals
                + 2 * (m-1)             // new connectors
            ) / (m * 2.0)
        );

    Lit prec = lit_Undef;
    while(conflict.size() > 0) {
        CardinalityConstraint cc;
        
        int i = N;
        if(prec != lit_Undef) { cc.lits.push(prec); i--; }
        for(; i > 0; i--) {
            if(conflict.size() == 0) break;
            assert(!removed[var(conflict.last())]);
            removed[var(conflict.last())] = true;
            cc.lits.push(~conflict.last());
            conflict.pop();
        }
        assert(cc.size() > 0);
        cc.bound = cc.size()-1;
        
        if(conflict.size() > 0) cc.bound++;

        for(i = 0; i < cc.bound; i++) {
            newVar();
            setFrozen(nVars()-1, true);
            removed.push(false);
            cc.lits.push(~mkLit(nVars()-1));
            if(i != 0) addClause(~mkLit(nVars()-2), mkLit(nVars()-1)); // symmetry breaker
            if(i == 0 && conflict.size() > 0) { 
                prec = mkLit(nVars()-1);
            }
            else {
                softLits.push(mkLit(nVars()-1));
                nonInputSoftLits.push(softLits.last());
            }
        }
        
        addConstraint(cc);
    }
    assert(conflict.size() == 0);
    
    int newLen = 0;
    for(int i = 0; i < softLits.size(); i++) {
        softLits[newLen] = softLits[i];
        if(!removed[var(softLits[i])]) newLen++;
    }
    softLits.shrink(softLits.size()-newLen);
}

lbool CircumscriptionSolver::sat() const {
    cerr << "c solving end" << endl;
    return l_True;
}

lbool CircumscriptionSolver::unsat() const {
    cout << "UNSATISFIABLE" << endl;
    cerr << "c solving end" << endl;
    return l_False;
}

lbool CircumscriptionSolver::solve(int n) {
    cerr << "c solving start" << endl;
    
    requiredModels = n;
    requiredWitnesses = option_circumscription_witnessess;
    
    assert(softLits.size() == 0);
    objectLits.copyTo(softLits);
    removed.growTo(nVars(), false);
    
    vec<Lit> blockingClause;
    
    while(true) {
        assert(decisionLevel() == 0);
        assumptions.clear();
        softLits.copyTo(assumptions);
        PseudoBooleanSolver::solve();
        if(status == l_False) {
            if(conflict.size() == 0) break;
            cancelUntil(0);
            shrink();
            one();
            //pmres();
            //kdyn();
            continue;
        }
        
        assert(status == l_True);

        assert(blockingClause.size() == 0);
        for(int i = 0; i < nCLits(); i++) blockingClause.push(value(getCLit(i)) == l_True ? ~getCLit(i) : getCLit(i));
        for(int i = 0; i < nOLits(); i++) if(value(getOLit(i)) == l_False) blockingClause.push(getOLit(i));

        if(requiredWitnesses == 1) onComputedModel();
        else {
            assumptions.clear();
            for(int i = 0; i < nCLits(); i++) assumptions.push(value(getCLit(i)) == l_True ? getCLit(i) : ~getCLit(i));
            for(int i = 0; i < nOLits(); i++) assumptions.push(value(getOLit(i)) == l_True ? getOLit(i) : ~getOLit(i));
            cerr << "c enumerate start" << endl;
            enumerate();
            cerr << "c enumerate end" << endl;
        }
        if(computedModels == requiredModels) break;
        
        cancelUntil(0);
        //cout << "BC " << blockingClause << endl;
        addClause(blockingClause);
        blockingClause.clear();
    }
    
    return computedModels > 0 ? sat() : unsat();
}

bool CircumscriptionSolver::checkConsistencyAndRemoveDuplicates(vec<Lit>& lits) {
    if(lits.size() == 0) return true;
    vec<int> tmp;
    tmp.growTo(nVars(), 0);
    int newLen = 0;
    for(int i = 0; i < lits.size(); i++) {
        lits[newLen] = lits[i];
        if(tmp[var(lits[i])] == 0) {
            newLen++;
            tmp[var(lits[i])] = sign(lits[i]) ? -1 : 1;
        }
        else if(tmp[var(lits[i])] != (sign(lits[i]) ? -1 : 1))
            return false;
    }
    lits.shrink_(lits.size()-newLen);
    return true;
}

bool CircumscriptionSolver::disjoint(const vec<Lit>& a, const vec<Lit>& b) {
    if(a.size() == 0 || b.size() == 0) return true;
    vec<bool> in(nVars(), false);
    for(int i = 0; i < a.size(); i++) in[var(a[i])] = true;
    for(int i = 0; i < b.size(); i++) if(in[var(b[i])]) return false;
    return true;
}

void CircumscriptionSolver::onComputedModel() {
    computedModels++;
    cerr << "c model " << computedModels << endl;
    if(!option_print_model) return;
    extendModel();
    cout << "v";
    for(int i = 0; i < visible.size(); i++) if(value(visible[i]) == l_True) cout << " " << visibleValue[i];
    cout << endl;
}

void CircumscriptionSolver::extendModel() {
    int i, j;
    Lit x;
    for(i = elimclauses.size()-1; i > 0; i -= j) {
        for(j = elimclauses[i--]; j > 1; j--, i--) if(value(Glucose::toLit(elimclauses[i])) != l_False) goto next;
        x = Glucose::toLit(elimclauses[i]);
        assert(var(x) < assigns.size());
        assigns[var(x)] = lbool(!sign(x));
    next:;
    }
}

} // namespace aspino

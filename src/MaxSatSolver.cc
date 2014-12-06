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

#include "MaxSatSolver.h"

#include <mtl/Map.h>
#include <core/Dimacs.h>

#include "utils/algorithm.h"

using Glucose::Map;

bool validate_maxsat_strat(const char* name, const string& value) {
    if(value == "one") return true;
    if(value == "pmres") return true;
    if(value == "pmres_split_conj") return true;
    if(value == "pmres_log") return true;
    cerr << "Invalid value for --" << name << ": " << value << "\n";
    return false;
}
DEFINE_string(maxsat_strat, "one", "Set optimization strategy. Valid values: one, pmres, pmres_log, pmres_split_conj.");

bool validate_maxsat_disjcores(const char* name, const string& value) {
    if(value == "no") return true;
    if(value == "pre") return true;
    if(value == "all") return true;
    cerr << "Invalid value for --" << name << ": " << value << "\n";
    return false;
}
DEFINE_string(maxsat_disjcores, "pre", "Set disjunct unsatisfiable cores policy. Valid values: no, pre, all.");

DEFINE_bool(maxsat_saturate, false, "Eliminate all cores of weight W before considering any core of level smaller than W.");

namespace aspino {

template<class B>
static long parseLong(B& in) {
    long    val = 0;
    bool    neg = false;
    skipWhitespace(in);
    if      (*in == '-') neg = true, ++in;
    else if (*in == '+') ++in;
    if (*in < '0' || *in > '9') fprintf(stderr, "PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    while (*in >= '0' && *in <= '9')
        val = val*10 + (*in - '0'),
        ++in;
    return neg ? -val : val; }

    
MaxSatSolver::MaxSatSolver() : lowerbound(0) {
    if(FLAGS_maxsat_strat == "one") corestrat = &MaxSatSolver::corestrat_one;
    else if(FLAGS_maxsat_strat == "pmres") corestrat = &MaxSatSolver::corestrat_pmres;
    else if(FLAGS_maxsat_strat == "pmres_split_conj") corestrat = &MaxSatSolver::corestrat_pmres_split_conj;
    else if(FLAGS_maxsat_strat == "pmres_log") corestrat = &MaxSatSolver::corestrat_pmreslog;
    else assert(0);
    
    if(FLAGS_maxsat_disjcores == "no") disjcores = NO;
    else if(FLAGS_maxsat_disjcores == "pre") disjcores = PRE;
    else if(FLAGS_maxsat_disjcores == "all") disjcores = ALL;
    else assert(0);
    
    saturate = FLAGS_maxsat_saturate;
}

MaxSatSolver::~MaxSatSolver() {
}

void MaxSatSolver::sameSoftVar(Lit soft, long weight) {
    assert(weights[var(soft)] != 0);
    assert(decisionLevel() == 0);
    int pos = 0;
    for(int i = 0; i < softLiterals.size(); i++, pos++) if(var(softLiterals[i]) == var(soft)) break;
    
    if(softLiterals[pos] == soft) {
        weights[var(soft)] += weight;
        upperbound += weight;
        return;
    }
        
    if(weights[var(soft)] == weight) {
        lowerbound += weight;
        setFrozen(var(soft), false);
        int j = 0;
        for(int i = 0; i < softLiterals.size(); i++) {
            if(softLiterals[i] == ~soft) continue;
            softLiterals[j++] = softLiterals[i];
        }
        assert(softLiterals.size() == j+1);
        softLiterals.shrink_(1);
    }else if(weights[var(soft)] < weight) {
        lowerbound += weights[var(soft)];
        for(int i = 0; i < softLiterals.size(); i++) if(softLiterals[i] == ~soft) { softLiterals[i] = soft; break; }
        weights[var(soft)] = weight - weights[var(soft)];
        upperbound += weights[var(soft)];
    }
    else {
        assert(weights[var(soft)] > weight);
        lowerbound += weight;
        weights[var(soft)] -= weight;
        upperbound += weights[var(soft)];
    }
}

void MaxSatSolver::addWeightedClause(vec<Lit>& lits, long weight) {
    Lit soft;
    if(lits.size() == 1)
        soft = lits[0];
    else {
        newVar();
        weights.push(0);
        soft = mkLit(nVars()-1);
        lits.push(~soft);
        addClause_(lits);
    }

    assert(weights.size() == nVars());
    if(weights[var(soft)] != 0) {
        sameSoftVar(soft, weight);
        return;
    }
    
    softLiterals.push(soft);
    assert_msg(weights[var(soft)] == 0, lits);
    weights[var(soft)] = weight;
    upperbound += weight;
    setFrozen(var(soft), true);
}
    
void MaxSatSolver::parse(gzFile in_) {
    Minisat::StreamBuffer in(in_);

    bool weighted = false;
    long top = -1;
    long weight = 1;
    
    vec<Lit> lits;
    int vars = 0;
    int clauses = 0;
    int count = 0;
    for(;;) {
        skipWhitespace(in);
        if(*in == EOF) break;
        if(*in == 'p') {
            ++in;
            if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            ++in;
            if(*in == 'w') { weighted = true; ++in; }
            
            if(eagerMatch(in, "cnf")) {
                vars = parseInt(in);
                clauses = parseInt(in);
                if(weighted && *in != '\n') top = parseInt(in);
                
                nInVars(vars);
                while(weights.size() < nInVars()) weights.push(0);
                while(nVars() < nInVars()) newVar();
            }
            else {
                cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            }
        }
        else if(*in == 'c')
            skipLine(in);
        else {
            count++;
            if(weighted) weight = parseLong(in);
            readClause(in, *this, lits);
            if(weight == top) addClause_(lits);
            else addWeightedClause(lits, weight);
        }
    }
    if(count != clauses)
        cerr << "WARNING! DIMACS header mismatch: wrong number of clauses." << endl;
}

long MaxSatSolver::setAssumptions(long limit) {
    assumptions.clear();
    long next = limit;
    int j = 0;
    for(int i = 0; i < softLiterals.size(); i++) {
        long w = weights[var(softLiterals[i])];
        if(w == 0) continue;
        if(w + lowerbound > upperbound) {
            addClause(softLiterals[i]);
            trace(maxsat, 10, "Hardening of " << softLiterals[i] << " of weight " << weights[var(softLiterals[i])]);
            weights[var(softLiterals[i])] = 0;
            continue;
        }
        softLiterals[j++] = softLiterals[i];
        if(var(softLiterals[i]) >= lastSoftLiteral) continue;
        if(w >= limit)
            assumptions.push(softLiterals[i]);
        else if(next == limit || w > next)
            next = w;
    }
    softLiterals.shrink_(softLiterals.size()-j);
    return next;
}

lbool MaxSatSolver::solve() {
//    sort();
    detectLevels();

    assert(levels.size() > 0);

    while(lowerbound < upperbound) {
        assert(levels.size() > 0);
        lbool ret = solveCurrentLevel();
        if(ret == l_False) {
            cout << "s UNSATISFIABLE" << endl;
            return l_False;
        }
        assert(ret == l_True);
        
        for(int i = 0; i < softLiterals.size(); i++) {
            assert(weights[var(softLiterals[i])] + lowerbound > upperbound);
            addClause(softLiterals[i]);
            trace(maxsat, 10, "Hardening of " << softLiterals[i] << " of weight " << weights[var(softLiterals[i])]);
            weights[var(softLiterals[i])] = 0;
        }
    }
    
    while(levels.size() > 0) { delete levels.last(); levels.pop(); }
    
    cout << "s OPTIMUM FOUND" << endl;
    copyModel();
    printModel();
    return l_True;
}

lbool MaxSatSolver::solveCurrentLevel() {
    trace(maxsat, 1, "Solve level " << levels.size());
    levels.last()->moveTo(softLiterals);
    delete levels.last();
    levels.pop();
    
    lastSoftLiteral = disjcores == NO ? INT_MAX : nVars();
    firstLimit = LONG_MAX;
    
    int iteration = 1;

    trace(maxsat, 2, "Iteration " << iteration);
    lbool ret = solve_();
    trace(maxsat, 2, "Bounds after iteration " << iteration << ": [" << lowerbound << ":" << upperbound << "]");
    
    if(ret == l_False) return l_False;
    assert(ret == l_True);
    
    while(lastSoftLiteral < nVars()) {
        iteration++;
        lastSoftLiteral = disjcores == ALL ? nVars() : INT_MAX;
        trace(maxsat, 2, "Iteration " << iteration);
        ret = solve_();
        trace(maxsat, 2, "Bounds after iteration " << iteration << ": [" << lowerbound << ":" << upperbound << "]");
        assert(ret == l_True);
    }

    return l_True;
}

lbool MaxSatSolver::solve_() {
    long limit = firstLimit;
    long nextLimit;
    bool foundCore = false;
    
    for(;;) {
        //aspino::shuffle(softLiterals);
        nextLimit = setAssumptions(limit);
        
        trace(maxsat, 2, "Solve with " << assumptions.size() << " assumptions. Current bounds: [" << lowerbound << ":" << upperbound << "]");
        trace(maxsat, 10, "Assumptions: " << assumptions);
        PseudoBooleanSolver::solve();
        
        if(status == l_True) {
            if(saturate && lastSoftLiteral < nVars()) {
                lastSoftLiteral = nVars();
                trace(maxsat, 8, "Continue on limit " << limit << " considering " << (nVars() - lastSoftLiteral) << " more literals");
                continue;
            }
            
            if(nextLimit == limit) {
                updateUpperBound();
                return l_True;
            }
            
            trace(maxsat, 4, "SAT! Decrease limit to " << nextLimit);
            limit = nextLimit;
            if(!foundCore) firstLimit = limit;

            updateUpperBound();

            continue;
        }
        
        foundCore = true;
        
        trace(maxsat, 2, "UNSAT! Conflict of size " << conflict.size());
        trace(maxsat, 10, "Conflict: " << conflict);
        
        if(conflict.size() == 0) return l_False;
        
        trim();
        
        trace(maxsat, 4, "Analyze conflict of size " << conflict.size() << " and weight " << limit);
        lowerbound += limit;
        cout << "o " << lowerbound << endl;
        (this->*corestrat)(limit);
    }
}

void MaxSatSolver::updateUpperBound() {
    long newupperbound = lowerbound;
    for(int i = 0; i < levels.size(); i++) newupperbound += levels[i]->size() * weights[var(levels[i]->last())];
    for(int i = 0; i < softLiterals.size(); i++) {
        if(value(softLiterals[i]) == l_False) newupperbound += weights[var(softLiterals[i])];
    }
    if(newupperbound < upperbound) {
        upperbound = newupperbound;
        cout << "c ub " << upperbound << endl;
        cancelUntil(0); // new hardening are possible
    }
}

void MaxSatSolver::trim() {
    cancelUntil(0);
    
    if(conflict.size() > 1) {
        addClause(conflict);
        int oldSize;
        do{
            assumptions.clear();
            for(int i = 0; i < conflict.size(); i++) assumptions.push(~conflict[i]);
            oldSize = conflict.size();
            PseudoBooleanSolver::solve();
            assert(status == l_False);
            trace(maxsat, 4, "Trim " << oldSize - conflict.size() << " literals from conflict");
            trace(maxsat, 10, "Conflict: " << conflict);
            cancelUntil(0);
        }while(oldSize > conflict.size() && conflict.size() > 1);
    }
        
//        if(conflict.size() > 1) {
//            vec<Lit> copy;
//            int n = conflict.size()-1;
//            int iter = 0;
//            do{
//                iter++;
//                conflict.moveTo(copy);
//                assumptions.clear();
//                for(int i = 0; i < copy.size(); i++) assumptions.push(~copy[i]);
//                for(int i = 0; i < assumptions.size(); i++) {
//                    int j = rand() % assumptions.size();
//                    Lit tmp = assumptions[i];
//                    assumptions[i] = assumptions[j];
//                    assumptions[j] = tmp;
//                }
//                PseudoBooleanSolver::solve();
//                assert(status == l_False);
//                trace(maxsat, 2, "UNSAT! Conflict of size (after trim)  " << conflict.size());
//                trace(maxsat, 10, "Conflict: " << conflict);
//                cancelUntil(0);
//                trace(maxsat, 2, "DIFF " << copy.size() << " - " << conflict.size() << " = " << - conflict.size() + copy.size() << " - iteration " << iter);
//                if(iter == n) break;
//                if(conflict.size() == 1) break;
//            }while(true);//}while(copy.size() > conflict.size());
//        }
}

void MaxSatSolver::corestrat_one(long limit) {
    trace(maxsat, 10, "Use algorithm one");
    WeightConstraint wc;
    wc.bound = conflict.size() - 1;
    while(conflict.size() > 0) {
        weights[var(conflict.last())] -= limit;
        wc.lits.push(~conflict.last());
        wc.coeffs.push(1);
        conflict.pop();
    }
    assert(conflict.size() == 0);
    for(int i = 0; i < wc.bound; i++) {
        newVar();
        if(i != 0) addClause(~softLiterals.last(), mkLit(nVars()-1));
        setFrozen(nVars()-1, true);
        weights.push(limit);
        softLiterals.push(mkLit(nVars()-1));
        wc.lits.push(~mkLit(nVars()-1));
        wc.coeffs.push(1);
    }
    
    if(wc.size() > 1) addConstraint(wc);
}

void MaxSatSolver::corestrat_pmres(long limit) {
    trace(maxsat, 10, "Use algorithm pmres");

    Lit prec = lit_Undef;
    vec<Lit> lits;
    while(conflict.size() > 0) {
        weights[var(conflict.last())] -= limit;
        
        if(prec == lit_Undef) prec = ~conflict.last();
        else {
            // disjunction
            newVar();
            setFrozen(nVars()-1, true);
            weights.push(limit);
            softLiterals.push(mkLit(nVars()-1));
            lits.push(prec);
            lits.push(~conflict.last());
            lits.push(mkLit(nVars()-1, true));
            addClause(lits);
            lits.clear();
            addClause(~prec, mkLit(nVars()-1));
            addClause(conflict.last(), mkLit(nVars()-1));
            
            if(conflict.size() == 0) break;

            // conjunction
            newVar();
            weights.push(0);
            addClause(prec, mkLit(nVars()-1, true));
            addClause(~conflict.last(), mkLit(nVars()-1, true));
            lits.push(~prec);
            lits.push(conflict.last());
            lits.push(mkLit(nVars()-1));
            addClause(lits);
            lits.clear();
            
            prec = mkLit(nVars()-1);
        }
        
        conflict.pop();
    }
    assert(conflict.size() == 0);
}

void MaxSatSolver::corestrat_pmres_split_conj(long limit) {
    trace(maxsat, 10, "Use algorithm pmres_split_conj");

    vec<Lit> lits;
    for(int i = 0; i < conflict.size(); i++) {
        weights[var(conflict[i])] -= limit;
        
        if(i == 0) continue;
        
        newVar();
        setFrozen(nVars()-1, true);
        weights.push(limit);
        softLiterals.push(mkLit(nVars()-1));
        
        for(int j = 0; j < i; j++) {
            lits.push(mkLit(nVars()-1, true));
            lits.push(~conflict[j]);
            lits.push(~conflict[i]);
            addClause_(lits);
            lits.clear();
        }
        
        addClause(mkLit(nVars()-1), conflict[i]);
        lits.push(mkLit(nVars()-1));
        for(int j = 0; j < i; j++) lits.push(conflict[j]);
        addClause_(lits);
        lits.clear();
    }
}

void MaxSatSolver::corestrat_pmreslog(long limit) {
    trace(maxsat, 10, "Use algorithm pmreslog");

    vec<Lit> lits;
    while(conflict.size() > 1) {
        int j = 0;
        for(int i = 0; i < conflict.size(); i += 2) {
            if(i+1 == conflict.size()) { conflict[j++] = conflict[i]; break; }
            
            weights[var(conflict[i])] -= limit;
            assert(weights[var(conflict[i])] >= 0);
            weights[var(conflict[i+1])] -= limit;
            assert(weights[var(conflict[i+1])] >= 0);
            
            // disjunction
            newVar();
            setFrozen(nVars()-1, true);
            weights.push(limit);
            softLiterals.push(mkLit(nVars()-1));
            lits.push(~conflict[i]);
            lits.push(~conflict[i+1]);
            lits.push(mkLit(nVars()-1, true));
            addClause(lits);
            lits.clear();
            addClause(conflict[i], mkLit(nVars()-1));
            addClause(conflict[i+1], mkLit(nVars()-1));
            
            if(conflict.size() <= 2) return;
            
            // conjunction
            newVar();
            weights.push(limit);
            addClause(~conflict[i], mkLit(nVars()-1, true));
            addClause(~conflict[i+1], mkLit(nVars()-1, true));
            lits.push(conflict[i]);
            lits.push(conflict[i+1]);
            lits.push(mkLit(nVars()-1));
            addClause(lits);
            lits.clear();
            
            conflict[j++] = mkLit(nVars()-1, true);
        }
        conflict.shrink_(conflict.size()-j);
    }
    assert(conflict.size() == 1);
    weights[var(conflict[0])] -= limit;
}

void MaxSatSolver::detectLevels() {
    vec<long> allWeights;
    Map<long, vec<Lit>*> wMap;
    for(int i = 0; i < softLiterals.size(); i++) {
        Lit lit = softLiterals[i];
        long w = weights[var(lit)];
        if(!wMap.has(w)) { wMap.insert(w, new vec<Lit>()); allWeights.push(w); }
        wMap[w]->push(lit);
    }
    
    int n = allWeights.size();
    while(n > 0) {
        int newn = 0;
        for(int i = 1; i < n; i++) {
            if(allWeights[i-1] > allWeights[i]) {
                long tmp = allWeights[i-1];
                allWeights[i-1] = allWeights[i];
                allWeights[i] = tmp;
                newn = i;
            }
        }
        n = newn;
    }
    
    assert(levels.size() == 0);
    
    long cumulative = 0;
    for(int i = 0; i < allWeights.size(); i++) {
        long w = allWeights[i];
        if(w > cumulative) levels.push(new vec<Lit>());
        vec<Lit>& v = *wMap[w];
        for(int j = 0; j < v.size(); j++) levels.last()->push(v[j]);
        cumulative += w * v.size();
        delete wMap[w];
    }
}

} // namespace aspino

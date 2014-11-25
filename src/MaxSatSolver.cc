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

#include "core/Dimacs.h"
#include "utils/algorithm.h"

bool validate_maxsatstrat(const char* name, const string& value) {
    if(value == "one") return true;
    if(value == "pmres") return true;
    if(value == "pmres_split_conj") return true;
    if(value == "pmreslog") return true;
    cerr << "Invalid value for --" << name << ": " << value << "\n";
    return false;
}
DEFINE_string(maxsatstrat, "one", "Set optimization strategy. Valid values: one, pmres, pmreslog.");


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
    if(FLAGS_maxsatstrat == "one") corestrat = &MaxSatSolver::corestrat_one;
    else if(FLAGS_maxsatstrat == "pmres") corestrat = &MaxSatSolver::corestrat_pmres;
    else if(FLAGS_maxsatstrat == "pmres_split_conj") corestrat = &MaxSatSolver::corestrat_pmres_split_conj;
    else if(FLAGS_maxsatstrat == "pmreslog") corestrat = &MaxSatSolver::corestrat_pmreslog;
    else assert(0);
}

MaxSatSolver::~MaxSatSolver() {
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
    softLiterals.push(soft);
    assert(weights[var(soft)] == 0);
    weights[var(soft)] = weight;
    upperbound += weight;
    setFrozen(var(soft), true);
}
    
void MaxSatSolver::parse(gzFile in_) {
    Minisat::StreamBuffer in(in_);

    bool weighted = false;
    long top = -1;
    long weight = -1;
    
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
    sort();
    lastSoftLiteral = nVars();

    lbool ret = l_False;
    for(;;) {
        trace(maxsat, 1, "Preprocessing: determining disjoint unsatisfiable cores...");
        ret = solve_(LONG_MAX);
        trace(maxsat, 1, "Bounds after preprocessing: [" << lowerbound << ":" << upperbound << "]");
        
        if(ret == l_False) {
            cout << "s UNSATISFIABLE" << endl;
            return l_False;
        }
        
        assert(ret == l_True);
        
        do{
            lastSoftLiteral = nVars(); //INT_MAX;
            long limit = 0;
            for(int i = 0; i < softLiterals.size(); i++) if(weights[var(softLiterals[i])] > limit) limit = weights[var(softLiterals[i])];
            ret = solve_(limit);
            assert(ret == l_True);
        }while(lastSoftLiteral < nVars());
        cout << "s OPTIMUM FOUND" << endl;
        copyModel();
        printModel();
        return l_True;
    }
}

lbool MaxSatSolver::solve_(long limit) {
    long nextLimit;
    
    for(;;) {
        //aspino::shuffle(softLiterals);
        nextLimit = setAssumptions(limit);
        
        trace(maxsat, 2, "Solve with " << assumptions.size() << " assumptions. Current bounds: [" << lowerbound << ":" << upperbound << "]");
        trace(maxsat, 10, "Assumptions: " << assumptions);
        PseudoBooleanSolver::solve();
        
        if(status == l_True) {
            if(nextLimit == limit) return l_True;
            
            trace(maxsat, 4, "SAT! Decrease limit to " << nextLimit);
            limit = nextLimit;

            long newupperbound = lowerbound;
            for(int i = 0; i < softLiterals.size(); i++) {
                if(value(softLiterals[i]) == l_False) newupperbound += weights[var(softLiterals[i])];
            }
            if(newupperbound < upperbound) {
                trace(maxsat, 8, "Decrease upper bound to " << newupperbound);
                upperbound = newupperbound;
                cout << "c ub " << upperbound << endl;
                cancelUntil(0);
            }

            continue;
        }
        
        trace(maxsat, 2, "UNSAT! Conflict of size (before trim) " << conflict.size());
        trace(maxsat, 10, "Conflict: " << conflict);
        
        if(conflict.size() == 0) return l_False;
        cancelUntil(0);
        
        if(conflict.size() > 1) {
            addClause(conflict);
            vec<Lit> copy;
            int iter = 0;
            do{
                iter++;
                conflict.moveTo(copy);
                assumptions.clear();
                for(int i = 0; i < copy.size(); i++) assumptions.push(~copy[i]);
                PseudoBooleanSolver::solve();
                assert(status == l_False);
                trace(maxsat, 2, "UNSAT! Conflict of size (after trim)  " << conflict.size());
                trace(maxsat, 10, "Conflict: " << conflict);
                cancelUntil(0);
                trace(maxsat, 2, "DIFF " << copy.size() << " - " << conflict.size() << " = " << - conflict.size() + copy.size() << " - iteration " << iter);
            }while(copy.size() > conflict.size() && conflict.size() > 1);
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
        
        trace(maxsat, 4, "Analyze conflict of weight " << limit);
        lowerbound += limit;
        trace(maxsat, 8, "Increase lower bound to " << lowerbound);
        cout << "o " << lowerbound << endl;
        (this->*corestrat)(limit);
    }
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

void MaxSatSolver::sort() {
    int n = softLiterals.size();
    while(n > 0) {
        int newn = 0;
        for(int i = 1; i < n; i++) {
            Lit lprec = softLiterals[i-1];
            Lit lcurr = softLiterals[i];
            if(weights[var(lprec)] < weights[var(lcurr)]) {
                softLiterals[i-1] = lcurr;
                softLiterals[i] = lprec;
                newn = i;
            }
        }
        n = newn;
    }
}

//void MaxSatSolver::detectLevels() {
//    levels.growTo(weights.size(), 0);
//    if(softLiterals.size() == 0) return;
//
//    currentLevel = 0;
//    long cumulative = 0;
//    for(int i = softLiterals.size() - 1; i >= 0; i--) {
//        Lit lit = softLiterals[i];
//        long w = weights[var(lit)];
//        if(w > cumulative) currentLevel++;
//        levels[var(lit)] = currentLevel;
//        cumulative += w;
//    }
//}

} // namespace aspino

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
#include <utils/System.h>
#include <core/Dimacs.h>

#include "utils/algorithm.h"

using Glucose::Map;

Glucose::EnumOption option_maxsat_strat("MAXSAT", "maxsat-strat", "Set optimization strategy.", "one|one-2|one-neg|one-wc|one-neg-wc|one-pmres|one-pmres-2|pmres|pmres-reverse|pmres-log|pmres-split-conj|kdyn");
Glucose::EnumOption option_maxsat_disjcores("MAXSAT", "maxsat-disjcores", "Set disjunct unsatisfiable cores policy.", "no|pre|all", 1);

Glucose::BoolOption option_maxsat_printmodel("MAXSAT", "maxsat-print-model", "Print optimal model if found.", true);

Glucose::IntOption option_maxsat_tag("MAXSAT", "maxsat-tag", "Parameter for maxsat-strat.", 16, Glucose::IntRange(2, INT32_MAX));

namespace aspino {

//#include <signal.h>
//static MaxSatSolver* solver;
//static void watchdog(int) {
//    assert(solver != NULL);
//    solver->onTick();
//    alarm(1);
//}
//void MaxSatSolver::onTick() {
//    trace(maxsat, 25, "Tick: " << Glucose::cpuTime() << "; LBD ratio: " << (sumLBD / conflictsRestarts));
//    if(timeBudget > 0.0 && Glucose::cpuTime() > timeBudget) {
//        trace(maxsat, 1, "interrupt!");
//        Solver::asynch_interrupt = true;
//    }
//}
    
MaxSatSolver::MaxSatSolver() : lowerbound(0), lastConflict(0), lastPropagation(0), lastCallCpuTime(0.0) /*, timeBudget(0.0)*/ {
    if(strcmp(option_maxsat_strat, "one") == 0) corestrat = &MaxSatSolver::corestrat_one;
    else if(strcmp(option_maxsat_strat, "one-2") == 0) corestrat = &MaxSatSolver::corestrat_one_2;
    else if(strcmp(option_maxsat_strat, "one-neg") == 0) corestrat = &MaxSatSolver::corestrat_one_neg;
    else if(strcmp(option_maxsat_strat, "one-wc") == 0) corestrat = &MaxSatSolver::corestrat_one_wc;
    else if(strcmp(option_maxsat_strat, "one-neg-wc") == 0) corestrat = &MaxSatSolver::corestrat_one_neg_wc;
    else if(strcmp(option_maxsat_strat, "one-pmres") == 0) corestrat = &MaxSatSolver::corestrat_one_pmres;
    else if(strcmp(option_maxsat_strat, "one-pmres-2") == 0) corestrat = &MaxSatSolver::corestrat_one_pmres_2;
    else if(strcmp(option_maxsat_strat, "pmres") == 0) corestrat = &MaxSatSolver::corestrat_pmres;
    else if(strcmp(option_maxsat_strat, "pmres-reverse") == 0) corestrat = &MaxSatSolver::corestrat_pmres_reverse;
    else if(strcmp(option_maxsat_strat, "pmres-split-conj") == 0) corestrat = &MaxSatSolver::corestrat_pmres_split_conj;
    else if(strcmp(option_maxsat_strat, "pmres-log") == 0) corestrat = &MaxSatSolver::corestrat_pmreslog;
    else if(strcmp(option_maxsat_strat, "kdyn") == 0) corestrat = &MaxSatSolver::corestrat_kdyn;
    else assert(0);
    
    if(strcmp(option_maxsat_disjcores, "no") == 0) disjcores = NO;
    else if(strcmp(option_maxsat_disjcores, "pre") == 0) disjcores = PRE;
    else if(strcmp(option_maxsat_disjcores, "all") == 0) disjcores = ALL;
    else assert(0);
    
    setIncrementalMode();
}

MaxSatSolver::~MaxSatSolver() {
}

void MaxSatSolver::sameSoftVar(Lit soft, int64_t weight) {
    assert(weights[var(soft)] != 0);
    assert(decisionLevel() == 0);
    int pos = 0;
    for(int i = 0; i < softLiterals.size(); i++, pos++) if(var(softLiterals[i]) == var(soft)) break;
    assert(pos < softLiterals.size());
    
    if(softLiterals[pos] == soft) {
        weights[var(soft)] += weight;
        return;
    }
        
    if(weights[var(soft)] == weight) {
        updateLowerBound(weight);
        setFrozen(var(soft), false);
        softLiterals[pos] = softLiterals[softLiterals.size()-1];
        softLiterals.shrink_(1);
        weights[var(soft)] = 0;
    }
    else if(weights[var(soft)] < weight) {
        updateLowerBound(weights[var(soft)]);
        softLiterals[pos] = soft;
        weights[var(soft)] = weight - weights[var(soft)];
    }
    else {
        assert(weights[var(soft)] > weight);
        updateLowerBound(weight);
        weights[var(soft)] -= weight;
    }
}

void MaxSatSolver::addWeightedClause(vec<Lit>& lits, int64_t weight) {
    if(weight == 0) return;

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
    setFrozen(var(soft), true);
}
    
void MaxSatSolver::parse(gzFile in_) {
    Glucose::StreamBuffer in(in_);

    bool weighted = false;
    int64_t top = -1;
    int64_t weight = 1;
    
    vec<Lit> lits;
    int vars = 0;
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
                inClauses = parseInt(in);
                if(weighted && *in != '\n') top = parseLong(in);
                
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
    if(count != inClauses)
        cerr << "WARNING! DIMACS header mismatch: wrong number of clauses." << endl, exit(3);
}

void MaxSatSolver::hardening() {
    cancelUntil(0);
    int j = 0;
    for(int i = 0; i < softLiterals.size(); i++) {
        int64_t w = weights[var(softLiterals[i])];
        if(w + lowerbound >= upperbound) {
            addClause(softLiterals[i]);
            trace(maxsat, 30, "Hardening of " << softLiterals[i] << " of weight " << w);
            weights[var(softLiterals[i])] = 0;
            continue;
        }
        softLiterals[j++] = softLiterals[i];
    }
    softLiterals.shrink_(softLiterals.size()-j);
}

void MaxSatSolver::setAssumptions(int64_t limit) {
    assumptions.clear();
    cancelUntil(0);
    int j = 0;
    for(int i = 0; i < softLiterals.size(); i++) {
        int64_t w = weights[var(softLiterals[i])];
        if(w == 0) continue;
        assert(w + lowerbound < upperbound);
        softLiterals[j++] = softLiterals[i];
        if(var(softLiterals[i]) >= lastSoftLiteral) continue;
        if(w >= limit) assumptions.push(softLiterals[i]);
    }
    softLiterals.shrink_(softLiterals.size()-j);
}

int64_t MaxSatSolver::computeNextLimit(int64_t limit) const {
    int64_t next = limit;
    for(int i = 0; i < softLiterals.size(); i++) {
        int64_t w = weights[var(softLiterals[i])];
        if(w == 0) continue;
        if(w >= limit) continue;
        if(next == limit || w > next) next = w;
    }
    return next;
}

void MaxSatSolver::preprocess() {
    assert(decisionLevel() == 0);
    if(softLiterals.size() == 0) return;
    trace(maxsat, 10, "Preprocessing");
    
    trace(maxsat, 20, "Preprocessing: cache signs of soft literals");
    vec<bool> signs(nVars());
    for(int i = 0; i < softLiterals.size(); i++) {
        if(weights[var(softLiterals[i])] != weights[var(softLiterals[0])]) {
            trace(maxsat, 10, "Preprocessing: detected weighted instance; skip preprocessing");
            return;
        }
        signs[var(softLiterals[i])] = sign(softLiterals[i]);
    }
    
    trace(maxsat, 20, "Preprocessing: partition clauses by increasing size");
    vec<vec<CRef>*> clausesPartition;
    vec<int> sizes;
    Map<int, int> sizeMap;
    for(int i = 0; i < clauses.size(); i++) {
        Clause& clause = ca[clauses[i]];
        assert(clause.size() >= 2);
        if(!sizeMap.has(clause.size())) { 
            sizes.push(clause.size());
            sizeMap.insert(clause.size(), clausesPartition.size());
            clausesPartition.push(new vec<CRef>());
        }
        clausesPartition[sizeMap[clause.size()]]->push(clauses[i]);
    }
    sort(sizes);

    trace(maxsat, 20, "Preprocessing: search for input clauses being cores");
    for(int i = 0; i < sizes.size(); i++) {
        trace(maxsat, 30, "Preprocessing: consider clauses of size " << sizes[i]);
        vec<CRef>& clauses = *clausesPartition[sizeMap[sizes[i]]];
        for(int j = 0; j < clauses.size(); j++) {
            Clause& clause = ca[clauses[j]];
            assert(clause.size() == sizes[i]);
            
            int64_t min = LONG_MAX;
            for(int k = 0; k < clause.size(); k++) {
                if(value(clause[k]) == l_False) continue;
                if(weights[var(clause[k])] == 0 || signs[var(clause[k])] == sign(clause[k])) { min = LONG_MAX; break; }
                if(weights[var(clause[k])] < min) min = weights[var(clause[k])];
            }
            if(min == LONG_MAX) continue;
            
            conflict.clear();
            for(int k = 0; k < clause.size(); k++) if(value(clause[k]) != l_False) conflict.push(clause[k]);
            trace(maxsat, 4, "Analyze conflict of size " << conflict.size() << " and weight " << min);
            updateLowerBound(min);
            (this->*corestrat)(min);
        }
    }

    trace(maxsat, 20, "Preprocessing: clean up");
    for(int i = 0; i < clausesPartition.size(); i++) delete clausesPartition[i];
}

void MaxSatSolver::removeSoftLiteralsAtLevelZero() {
    assert(decisionLevel() == 0);
    int j = 0;
    for(int i = 0; i < softLiterals.size(); i++) {
        if(value(softLiterals[i]) != l_Undef) { 
            if(value(softLiterals[i]) == l_False) {
                trace(maxsat, 20, "False soft literal at level 0: " << softLiterals[i] << " of weight " << weights[var(softLiterals[i])]);
                updateLowerBound(weights[var(softLiterals[i])]);
            }
            weights[var(softLiterals[i])] = 0;
            continue;
        }
        softLiterals[j++] = softLiterals[i];
    }
    softLiterals.shrink_(softLiterals.size() - j);
}

lbool MaxSatSolver::solve() {
//    solver = this;
//    signal(SIGALRM, watchdog);
//    alarm(1);

    inClauses = clauses.size();
    upperbound = LONG_MAX;

//    // spend some time to identify easy backbones literals
//    timeBudget = Glucose::cpuTime() + 10.0;
//    setConfBudget(1000);
    PseudoBooleanSolver::solve();
//    budgetOff();
//    clearInterrupt();
//    timeBudget = 0.0;
    if(status == l_False) { cout << "s UNSATISFIABLE" << endl; return l_False; }
    if(status == l_True) updateUpperBound();
    cancelUntil(0);

    removeSoftLiteralsAtLevelZero();
    preprocess();
    hardening();

    solve_();
    trace(maxsat, 2, "Bounds: [" << lowerbound << ":" << upperbound << "]");

    if(upperbound == LONG_MAX) { cout << "s UNSATISFIABLE" << endl; return l_False; }
        
    if(lowerbound < upperbound) {
        if(!ok) status = l_False;
        else PseudoBooleanSolver::solve();
        if(status == l_True) updateUpperBound();
        else lowerbound = upperbound;
    }
    
    assert(lowerbound == upperbound);
    assert(upperbound < LONG_MAX);
    
    cout << "o " << lowerbound << endl;
    cout << "s OPTIMUM FOUND" << endl;
    if(option_maxsat_printmodel) printModel();
    return l_True;
}

void MaxSatSolver::solve_() {
    int vars = nVars();
    lastSoftLiteral = disjcores == NO ? INT_MAX : vars;

    int64_t limit = computeNextLimit(LONG_MAX);
    int64_t nextLimit;
    
    for(;;) {
        hardening();
        setAssumptions(limit);
        lastConflict = conflicts;
        lastPropagation = propagations;
        lastCallCpuTime = Glucose::cpuTime();

        trace(maxsat, 2, "Solve with " << assumptions.size() << " assumptions. Current bounds: [" << lowerbound << ":" << upperbound << "]");
        trace(maxsat, 100, "Assumptions: " << assumptions);
        
        if(assumptions.size() == 0 && upperbound != LONG_MAX) status = l_Undef;
        else //progressionBinaryFind(limit);  
             PseudoBooleanSolver::solve(); //
        
        if(status != l_False) {
            if(status == l_True) updateUpperBound();
            
            if(lastSoftLiteral == vars && vars < nVars()) {
                lastSoftLiteral = INT_MAX; //disjcores == ALL ? nVars() : INT_MAX;
                trace(maxsat, 4, "Continue on limit " << limit);
                continue;
            }
            
            nextLimit = computeNextLimit(limit);
            if(nextLimit == limit) {
                trace(maxsat, 4, (status == l_True ? "SAT!" : "Skip!") << " No other limit to try");
                return;
            }
            
            trace(maxsat, 4, (status == l_True ? "SAT!" : "Skip!") << " Decrease limit to " << nextLimit);
            limit = nextLimit;
            lastSoftLiteral = disjcores == NO ? INT_MAX : vars;

            continue;
        }
        
        cancelUntil(0);
        trace(maxsat, 2, "UNSAT! Conflict of size " << conflict.size());
        trace(maxsat, 100, "Conflict: " << conflict);
        
        if(conflict.size() == 0) return;

        updateLowerBound(limit);
        
        assert(decisionLevel() == 0);
        progressionMinimize(limit);
//        biprogressionMinimize(limit);
//        binaryMinimize(limit);
//        progressionBinaryMinimize(limit);
        trim(); // last trim, just in case some new learned clause may help to further reduce the core

        assert(conflict.size() > 0);
        trace(maxsat, 4, "Analyze conflict of size " << conflict.size() << " and weight " << limit);
        if(conflict.size() == 1) weights[var(conflict.last())] = 0;
        else (this->*corestrat)(limit);
    }
}

void MaxSatSolver::updateLowerBound(int64_t limit) {
    lowerbound += limit;
    cout << "o " << lowerbound << endl;
}

void MaxSatSolver::updateUpperBound() {
    assert(status == l_True);
    int64_t newupperbound = lowerbound;
    for(int i = 0; i < softLiterals.size(); i++)
        if(value(softLiterals[i]) == l_False) newupperbound += weights[var(softLiterals[i])];
    if(newupperbound < upperbound) {
        upperbound = newupperbound;
        copyModel();
        trace(maxsat, 200, "Model: " << model);
        cout << "c ub " << upperbound << endl;
    }
}

void MaxSatSolver::progressionBinaryFind(int64_t limit) {
    assert(decisionLevel() == 0);
    
    vec<Lit> allAssumptions;
    assumptions.moveTo(allAssumptions);
    vec<Lit> remainingAssumptions;
    
    satisfiableSoftLiterals.growTo(nVars(), false);
    for(int i = 0; i < allAssumptions.size(); i++) {
        assert(var(allAssumptions[i]) < satisfiableSoftLiterals.size());
        if(satisfiableSoftLiterals[var(allAssumptions[i])]) 
            assumptions.push(allAssumptions[i]);
        else
            remainingAssumptions.push(allAssumptions[i]);
    }
    assumptions.copyTo(allAssumptions);
    for(int i = 0; i < remainingAssumptions.size(); i++) allAssumptions.push(remainingAssumptions[i]);
    const int progressionFrom = 1;
    int progression = progressionFrom;
    int fixed = assumptions.size();
    int sat = fixed;
//*
    while(lowerbound + limit < upperbound) {
        if(fixed + progression >= allAssumptions.size()) {
            assert(sat <= assumptions.size());
            assumptions.shrink_(assumptions.size() - sat);
            break;
        }

        trace(maxsat, 15, "Find: progress to " << progression << "; fixed = " << fixed << "; LBD ratio: " << (sumLBD / conflictsRestarts));
        
        for(int i = assumptions.size(); i < fixed + progression; i++) {
            assert(i < allAssumptions.size());
            assumptions.push(allAssumptions[i]);
        }

//        float sumLBD_ = sumLBD;
//        uint64_t conflictsRestarts_ = conflictsRestarts;
        sumLBD = 0.0;
        conflictsRestarts = 0;
//        timeBudget = Glucose::cpuTime() + 3.0;
        assert(decisionLevel() == 0);
        PseudoBooleanSolver::solve();
//        clearInterrupt();
//        timeBudget = 0.0;
        
        if(status == l_False) {
            trace(maxsat, 10, "Find: core of size " << conflict.size());
            cancelUntil(0);
            return;
        }
        else {
            trace(maxsat, 20, (status == l_True ? "SAT!" : "UNDEF"));
            if(status == l_True) {
                for(int i = sat; i < assumptions.size(); i++) satisfiableSoftLiterals[var(assumptions[i])] = true;
                sat = assumptions.size();
                updateUpperBound();
            }
//            else if(conflictsRestarts_ > 0 && ((sumLBD_*0.97 / conflictsRestarts_) > (sumLBD / conflictsRestarts))) {
//                trace(maxsat, 10, "Find: improve LBD from " << (sumLBD_ / conflictsRestarts_) << " to " << (sumLBD / conflictsRestarts));
//                continue;            
//            }
//            else {
//                trace(maxsat, 10, "Find: non improve LBD from " << (sumLBD_ / conflictsRestarts_) << " to " << (sumLBD / conflictsRestarts));
//            }
            progression *= 2;
        }
        cancelUntil(0);
    }
 //   */
//*
    fixed = assumptions.size();
    int toAdd = (allAssumptions.size() - fixed) / 2;
    while(lowerbound + limit < upperbound) {
        if(toAdd == 0) break;

        trace(maxsat, 15, "Find: try with " << fixed << "+" << toAdd << " literals; LBD ratio: " << (sumLBD / conflictsRestarts));
        
        for(int i = assumptions.size(); i < fixed + toAdd; i++) {
            assert(i < allAssumptions.size());
            assumptions.push(allAssumptions[i]);
        }
        
//        float sumLBD_ = sumLBD;
//        uint64_t conflictsRestarts_ = conflictsRestarts;
        sumLBD = 0.0;
        conflictsRestarts = 0;
//        timeBudget = Glucose::cpuTime() + 3.0;
        assert(decisionLevel() == 0);
        PseudoBooleanSolver::solve();
//        clearInterrupt();
//        timeBudget = 0.0;
        if(status == l_False) {
            trace(maxsat, 10, "Find: core of size " << conflict.size());
            return;
        }
        else {
            trace(maxsat, 20, (status == l_True ? "SAT!" : "UNDEF"));
            if(status == l_True) {
                updateUpperBound();
            }
//            else if(conflictsRestarts_ > 0 && ((sumLBD_*0.97 / conflictsRestarts_) > (sumLBD / conflictsRestarts))) {
//                trace(maxsat, 10, "Find: improve LBD from " << (sumLBD_ / conflictsRestarts_) << " to " << (sumLBD / conflictsRestarts));
//                continue;            
//            }
//            else {
//                trace(maxsat, 10, "Find: non improve LBD from " << (sumLBD_ / conflictsRestarts_) << " to " << (sumLBD / conflictsRestarts));
//            }
            fixed += toAdd;
            toAdd /= 2;
        }
        cancelUntil(0);
    }
//*/    
    sumLBD = 0.0;
    conflictsRestarts = 0;
    allAssumptions.moveTo(assumptions);
    trace(maxsat, 15, "Find: try with all " << assumptions.size() << " assumptions (no budget)");
    assert(decisionLevel() == 0);
    PseudoBooleanSolver::solve();
}

void MaxSatSolver::trim() {
    assert(decisionLevel() == 0);
    
    if(conflict.size() <= 1) return;

    int counter = 0;

    do{
        counter++;
        assumptions.clear();
        for(int i = 0; i < conflict.size(); i++) assumptions.push(~conflict[i]);
        PseudoBooleanSolver::solve();
        assert(status == l_False);
        trace(maxsat, 15, "Trim " << assumptions.size() - conflict.size() << " literals from conflict");
        trace(maxsat, 100, "Conflict: " << conflict);
        cancelUntil(0);
        if(conflict.size() <= 1) return;
    }while(assumptions.size() > conflict.size());
    
    if(counter % 2 == 1) for(int i = 0; i < assumptions.size(); i++) conflict[i] = ~assumptions[i];
    
    assert(conflict.size() > 1);
}

void MaxSatSolver::progressionMinimize(int64_t limit) {
    assert(decisionLevel() == 0);
    if(conflict.size() <= 1) return;
    
    double cpuTime = Glucose::cpuTime() - lastCallCpuTime;
    uint64_t budget = conflicts - lastConflict;
    uint64_t budgetMin = budget / cpuTime;
    if(budget > budget * 30 / cpuTime) budget = budget * 30 / cpuTime;
    trace(maxsat, 10, "Minimize core of size " << conflict.size() << " (" << (conflicts - lastConflict) << " conflicts; " << cpuTime << " seconds; each check with budget " << budget << ")");
    trim();
    if(budget == 0) return;

    vec<Lit> core;
    conflict.moveTo(core);
    
    vec<Lit> allAssumptions;
    for(int i = 0; i < core.size(); i++) allAssumptions.push(~core[i]);
    
    assumptions.clear();
    const int progressionFrom = 1;
    int progression = progressionFrom;
    int fixed = 0;
    while(lowerbound + limit < upperbound) {
        if(fixed + progression >= allAssumptions.size()) {
            if(progression == progressionFrom) break;
            progression = progressionFrom;
            fixed = assumptions.size();
            if(budget > budgetMin) budget /= 2;
            continue;
        }

        trace(maxsat, 15, "Minimize: progress to " << progression << "; fixed = " << fixed);
        
        int prec = assumptions.size();
        for(int i = assumptions.size(); i < fixed + progression; i++) {
            assert(i < allAssumptions.size());
            assumptions.push(allAssumptions[i]);
        }
        
        float sumLBD_ = sumLBD;
        uint64_t conflictsRestarts_ = conflictsRestarts;
        sumLBD = 0;
        conflictsRestarts = 0;
        setConfBudget(budget);
        PseudoBooleanSolver::solve();
        budgetOff();
        sumLBD += sumLBD_;
        conflictsRestarts += conflictsRestarts_;

        if(status == l_False) {
            trace(maxsat, 10, "Minimize: reduce to size " << conflict.size());
            progression = progressionFrom;
            if(budget > budgetMin) budget /= 2;
            
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
            assumptions.shrink(assumptions.size() - j);
            fixed = assumptions.size();
            
            j = 0;
            for(int i = 0, k = core.size() - 1; i < allAssumptions.size(); i++) {
                if(k < 0) break;
                if(allAssumptions[i] != ~core[k]) continue;
                allAssumptions[j++] = allAssumptions[i];
                k--;
            }
            allAssumptions.shrink(allAssumptions.size() - j);
        }
        else {
            trace(maxsat, 20, (status == l_True ? "SAT!" : "UNDEF"));
            progression *= 2;
            if(status == l_True) {
                updateUpperBound();
            }
        }
        cancelUntil(0);
    }
    core.moveTo(conflict);
}

void MaxSatSolver::biprogressionMinimize(int64_t limit) {
    assert(decisionLevel() == 0);
    if(conflict.size() <= 1) return;
    
    double cpuTime = Glucose::cpuTime() - lastCallCpuTime;
    uint64_t budget = conflicts - lastConflict;
    if(budget > budget * 30 / cpuTime) budget = budget * 30 / cpuTime;
    trace(maxsat, 10, "Minimize core of size " << conflict.size() << " (" << (conflicts - lastConflict) << " conflicts; " << cpuTime << " seconds; each check with budget " << budget << ")");
    trim();
    if(budget == 0) return;

    vec<Lit> core;
    conflict.moveTo(core);
    
    vec<Lit> allAssumptions;
    for(int i = 0; i < core.size(); i++) allAssumptions.push(~core[i]);
    
    assumptions.clear();
    const int progressionFrom = 1;
    int progression = progressionFrom;
    bool increasing = true;
    int fixed = 0;
    while(lowerbound + limit < upperbound) {
        if(increasing && fixed + progression >= allAssumptions.size()) {
            increasing = false;
//            if(progression == progressionFrom) break;
//            progression = progressionFrom;
            progression /= 2;
            fixed = assumptions.size();
//            if(budget > 100) budget /= 2;
            continue;
        }
        if(progression == 0) break;
        if(!increasing && fixed + progression >= allAssumptions.size()) { progression /= 2; continue; }

        trace(maxsat, 15, "Minimize: progress to " << progression << "; fixed = " << fixed);
        
        int prec = assumptions.size();
        for(int i = assumptions.size(); i < fixed + progression; i++) {
            assert(i < allAssumptions.size());
            assumptions.push(allAssumptions[i]);
        }
        
        float sumLBD_ = sumLBD;
        uint64_t conflictsRestarts_ = conflictsRestarts;
        sumLBD = 0;
        conflictsRestarts = 0;
        setConfBudget(budget);
        PseudoBooleanSolver::solve();
        budgetOff();
        sumLBD += sumLBD_;
        conflictsRestarts += conflictsRestarts_;

        if(status == l_False) {
            trace(maxsat, 10, "Minimize: reduce to size " << conflict.size());
//            progression = progressionFrom;
            increasing = false;
            progression /= 2;
            
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
            assumptions.shrink(assumptions.size() - j);
            fixed = assumptions.size();
            
            j = 0;
            for(int i = 0, k = core.size() - 1; i < allAssumptions.size(); i++) {
                if(k < 0) break;
                if(allAssumptions[i] != ~core[k]) continue;
                allAssumptions[j++] = allAssumptions[i];
                k--;
            }
            allAssumptions.shrink(allAssumptions.size() - j);
        }
        else {
            trace(maxsat, 20, (status == l_True ? "SAT!" : "UNDEF"));
            progression = increasing ? progression * 2 : progression / 2;
            if(!increasing && status == l_Undef) fixed = assumptions.size();
            if(status == l_True) updateUpperBound();
        }
        cancelUntil(0);
    }
    core.moveTo(conflict);
}

void MaxSatSolver::progressionBinaryMinimize(int64_t limit) {
    assert(decisionLevel() == 0);
    if(conflict.size() <= 1) return;
    
    double cpuTime = Glucose::cpuTime() - lastCallCpuTime;
    uint64_t budgetConflics = conflicts - lastConflict;
    if(budgetConflics > budgetConflics * 30 / cpuTime) budgetConflics = budgetConflics * 30 / cpuTime;
    trace(maxsat, 10, "Minimize core of size " << conflict.size() << " (" << (conflicts - lastConflict) << " conflicts); " << cpuTime << " seconds)");
    trim();
    if(conflicts == lastConflict) return;

    vec<Lit> core;
    conflict.moveTo(core);
    
    vec<Lit> allAssumptions;
    for(int i = 0; i < core.size(); i++) allAssumptions.push(~core[i]);
    
    assumptions.clear();
    const int progressionFrom = 1;
    int progression = progressionFrom;
    vec<Lit> sat;
    while(lowerbound + limit < upperbound) {
        if(progression >= allAssumptions.size()) break;

        trace(maxsat, 15, "Minimize: progress to " << progression << "; LBD ratio: " << (sumLBD / conflictsRestarts));
        
        int prec = assumptions.size();
        for(int i = assumptions.size(); i < progression; i++) {
            assert(i < allAssumptions.size());
            assumptions.push(allAssumptions[i]);
        }
        
        float sumLBD_ = sumLBD;
        uint64_t conflictsRestarts_ = conflictsRestarts;
        sumLBD = 0;
        conflictsRestarts = 0;
        setConfBudget(budgetConflics);
        PseudoBooleanSolver::solve();
        budgetOff();
        sumLBD += sumLBD_;
        conflictsRestarts += conflictsRestarts_;
        if(status == l_False) {
            trace(maxsat, 10, "Minimize: reduce to size " << conflict.size());
            
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
            assumptions.shrink(assumptions.size() - j);
            
            j = 0;
            for(int i = 0, k = core.size() - 1; i < allAssumptions.size(); i++) {
                if(k < 0) break;
                if(allAssumptions[i] != ~core[k]) continue;
                allAssumptions[j++] = allAssumptions[i];
                k--;
            }
            allAssumptions.shrink(allAssumptions.size() - j);

            j = 0;
            for(int i = 0, k = core.size() - 1; i < sat.size(); i++) {
                if(k < 0) break;
                if(sat[i] != ~core[k]) continue;
                sat[j++] = sat[i];
                k--;
            }
            sat.shrink(sat.size() - j);
            
            break;
        }

        trace(maxsat, 20, (status == l_True ? "SAT!" : "UNDEF"));
        progression *= 2;
        if(status == l_True) {
            while(sat.size() < assumptions.size()) sat.push(assumptions[sat.size()]);
            updateUpperBound();
        }
        cancelUntil(0);
    }

    sat.moveTo(assumptions);
    int fixed = assumptions.size();
    int toAdd = (core.size() - fixed) / 2;
    while(lowerbound + limit < upperbound) {
        if(toAdd <= core.size() / 10) break;

        trace(maxsat, 15, "Minimize: try with " << fixed << "+" << toAdd << " literals; LBD ratio: " << (sumLBD / conflictsRestarts));
        
        int prec = assumptions.size();
        for(int i = assumptions.size(); i < fixed + toAdd; i++) {
            assert_msg(i < allAssumptions.size(), "i=" << i <<"; allAssumptions.size()=" << allAssumptions.size());
            assumptions.push(allAssumptions[i]);
        }
        
        float sumLBD_ = sumLBD;
        uint64_t conflictsRestarts_ = conflictsRestarts;
        sumLBD = 0;
        conflictsRestarts = 0;
        setConfBudget(budgetConflics);
        PseudoBooleanSolver::solve();
        budgetOff();
        sumLBD += sumLBD_;
        conflictsRestarts += conflictsRestarts_;
        if(status == l_False) {
            trace(maxsat, 10, "Minimize: reduce to size " << conflict.size());
            
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
            assumptions.shrink(assumptions.size() - j);
            fixed = assumptions.size();

            j = 0;
            for(int i = 0, k = core.size() - 1; i < allAssumptions.size(); i++) {
                if(k < 0) break;
                if(allAssumptions[i] != ~core[k]) continue;
                allAssumptions[j++] = allAssumptions[i];
                k--;
            }
            allAssumptions.shrink(allAssumptions.size() - j);
            assert_msg(allAssumptions.size() == core.size(), allAssumptions << " " << core);
            
            toAdd = (core.size() - fixed) / 2;
        }
        else {
            trace(maxsat, 20, (status == l_True ? "SAT!" : "UNDEF"));
            fixed += toAdd;
            toAdd /= 2;
            if(status == l_True) updateUpperBound();
        }
        cancelUntil(0);
    }
    
    core.moveTo(conflict);
}

//void MaxSatSolver::initUpperBound() {
//    assumptions.clear();
//    vec<Lit> allAssumptions;
//    softLiterals.copyTo(allAssumptions);
//    const int budget = 10000;
//    const int progressionFrom = 1;
//    int progression = progressionFrom;
//    int fixed = 0;
//    while(true) {
//        if(fixed + progression > allAssumptions.size()) {
//            if(progression == progressionFrom) return;
//            progression = progressionFrom;
//            fixed = assumptions.size();
//            continue;
//        }
//
//        trace(maxsat, 150, "Init upper bound: progress to " << progression << "; fixed = " << fixed);
//        
//        for(int i = assumptions.size(); i < fixed + progression; i++) {
//            assert(i < allAssumptions.size());
//            assumptions.push(allAssumptions[i]);
//        }
//        
//        setConfBudget(budget);
//        PseudoBooleanSolver::solve();
//        budgetOff();
//        if(status == l_False) {
//            do{
//                trace(maxsat, 100, "Init upper bound: remove one assumption");
//                for(int i = assumptions.size(); i < allAssumptions.size(); i++) allAssumptions[i-1] = allAssumptions[i];
//                assumptions.pop();
//                allAssumptions.pop();
//                setConfBudget(budget);
//                PseudoBooleanSolver::solve();
//                budgetOff();
//            }while(status == l_False);
//            
//            progression = progressionFrom;
//            fixed = assumptions.size();
//        }
//        else {
//            trace(maxsat, 200, (status == l_True ? "SAT!" : "UNDEF"));
//            progression *= 2;
//            if(status == l_True) updateUpperBound();
//        }
//        cancelUntil(0);
//    }
//}

void MaxSatSolver::binaryMinimize(int64_t limit) {
    assert(decisionLevel() == 0);
    if(conflict.size() <= 1) return;
    
    double cpuTime = Glucose::cpuTime() - lastCallCpuTime;
    uint64_t budget = conflicts - lastConflict;
    if(budget > budget * 60 / cpuTime) budget = budget * 60 / cpuTime;
    trace(maxsat, 10, "Minimize core of size " << conflict.size() << " (" << (conflicts - lastConflict) << " conflicts; " << static_cast<int>(cpuTime) << " seconds; each check with budget " << budget << ")");
    trim();
    if(budget == 0) return;

    vec<Lit> core;
    conflict.moveTo(core);
    
    vec<Lit> allAssumptions;
    for(int i = 0; i < core.size(); i++) allAssumptions.push(~core[i]);

    assumptions.clear();
    int toAdd = core.size() / 2;
    int fixed = 0;
    while(lowerbound + limit < upperbound) {
        if(toAdd == 0) break;

        trace(maxsat, 15, "Minimize: try with " << fixed << "+" << toAdd << " literals");
        
        int prec = assumptions.size();
        for(int i = assumptions.size(); i < fixed + toAdd; i++) {
            assert(i < allAssumptions.size());
            assumptions.push(allAssumptions[i]);
        }
        
        setConfBudget(budget);
        PseudoBooleanSolver::solve();
        budgetOff();
        if(status == l_False) {
            trace(maxsat, 10, "Minimize: reduce to size " << conflict.size());
            
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
            assumptions.shrink(assumptions.size() - j);
            fixed = assumptions.size();
            
            j = 0;
            for(int i = 0, k = core.size() - 1; i < allAssumptions.size(); i++) {
                if(k < 0) break;
                if(allAssumptions[i] != ~core[k]) continue;
                allAssumptions[j++] = allAssumptions[i];
                k--;
            }
            allAssumptions.shrink(allAssumptions.size() - j);
            
            toAdd = (core.size() - fixed) / 2;
        }
        else {
            trace(maxsat, 20, (status == l_True ? "SAT!" : "UNDEF"));
            fixed += toAdd;
            toAdd /= 2;
            if(status == l_True) {
                updateUpperBound();
            }
        }
        cancelUntil(0);
    }
    core.moveTo(conflict);
}

void MaxSatSolver::minimize() {
    assert(decisionLevel() == 0);
    if(conflict.size() <= 1) return;
    
    static int64_t previousConflicts = 0;
    int64_t budget = 25 * (conflicts - previousConflicts) / 100;
    
    static vec<lbool> required;
    required.growTo(nVars(), l_False);
    for(int i = 0; i < conflict.size(); i++) required[var(conflict[i])] = l_False;
    
    vec<Lit> core;
    conflict.moveTo(core);

    int requiredAssumptions = 0;
    
    for(;;) {
        assumptions.shrink_(assumptions.size()-requiredAssumptions);
        int nUnknown = 0;
        for(int i = 0; i < core.size(); i++) {
            if(required[var(core[i])] != l_Undef) continue;
            assumptions.push(~core[i]);
            nUnknown++; 
        }
        Lit last = lit_Undef;
        for(int i = 0; i < core.size(); i++) {
            if(required[var(core[i])] != l_False) continue;
            if(last != lit_Undef) assumptions.push(~last);
            last = core[i];
        }
        if(last == lit_Undef) break;
        
        trace(maxsat, 10, "Try to minimize with " << assumptions.size() << " assumptions (" << requiredAssumptions << " required, ie. " << 100*requiredAssumptions/(assumptions.size()+1) << "% of the core; " << nUnknown << " unknowns, ie. " << 100*nUnknown/(assumptions.size()+1) << "% of the core)");
        setConfBudget(budget);
        PseudoBooleanSolver::solve();
        

        if(status == l_False) {
//            cancelUntil(0);
//            addClause(conflict);
            
            conflict.moveTo(core);
            
            reverse(core);
            
            continue;
        }
        
        required[var(last)] = status;
        if(status == l_True) {
            assumptions[requiredAssumptions++] = ~last;
        }
    }
    
    core.moveTo(conflict);
    cancelUntil(0);
    
    budgetOff();
    previousConflicts = conflicts;
}

//void MaxSatSolver::minimize() {
//    assert(decisionLevel() == 0);
//    assert(upperbound != LONG_MAX);
//    if(conflict.size() <= 1) return;
//
//    static int64_t previousConflicts = 0;
//    int64_t budget = 100 * (conflicts - previousConflicts) / 100;
//
//    static vec<lbool> required;
//    required.growTo(nVars(), l_False);
//    for(int i = 0; i < conflict.size(); i++) required[var(conflict[i])] = l_False;
//    
//    vec<Lit> core;
//    conflict.moveTo(core);
//    
//    assumptions.clear();
//    int coreSize = 0;
//    for(int i = 0; i < core.size(); i++) {
//        if(modelValue(core[i]) != l_True) { assumptions.push(~core[i]); continue; }
//        core[coreSize++] = core[i];
//    }
//    core.shrink_(core.size()-coreSize);
//    
//    int requiredAssumptions = assumptions.size();
//    
//    for(;;) {
//        if(core.size() <= 1) break;
//        
//        assumptions.shrink_(assumptions.size()-requiredAssumptions);
//        assumptions.push(~core[0]);
//        
//        trace(maxsat, 10, "Try to find a model with " << assumptions.size() << " assumptions (" << requiredAssumptions << " required, ie. " << 100*requiredAssumptions/(assumptions.size()+core.size()-1) << "% of the core)");
//        setConfBudget(budget);
//        PseudoBooleanSolver::solve();
//    
//        
//
//        if(status == l_False) {
//            cancelUntil(0);
//            trace(maxsat, 15, "UNSAT! The new conflict is minimal (size " << conflict.size() << ")");
//            return;
//        }
//        
//        assumptions.shrink_(assumptions.size()-requiredAssumptions);
//        coreSize = 0;
//        for(int i = 0; i < core.size(); i++) {
//            if(value(core[i]) != l_True) { assumptions.push(~core[i]); continue; }
//            core[coreSize++] = core[i];
//        }
//        requiredAssumptions = assumptions.size();
//        core.shrink_(core.size()-coreSize);
//    }
//    
//    assert(core.size() <= 1);
//    assert(conflict.size() == 0);
//    conflict.push(core.last());
//    for(int i = 0; i < assumptions.size(); i++) conflict.push(~assumptions[i]);
//    cancelUntil(0);
//    
//    budgetOff();
//    previousConflicts = conflicts;
//}


void MaxSatSolver::corestrat_one(int64_t limit) {
    assert(decisionLevel() == 0);
    trace(maxsat, 10, "Use algorithm one");
    CardinalityConstraint cc;
    cc.bound = conflict.size() - 1;
    while(conflict.size() > 0) {
        weights[var(conflict.last())] -= limit;
        cc.lits.push(~conflict.last());
        conflict.pop();
    }
    assert(conflict.size() == 0);
    for(int i = 0; i < cc.bound; i++) {
        newVar();
        if(i != 0) addClause(~softLiterals.last(), mkLit(nVars()-1));
        setFrozen(nVars()-1, true);
        weights.push(limit);
        softLiterals.push(mkLit(nVars()-1));
        cc.lits.push(~mkLit(nVars()-1));
    }
    
    assert(cc.size() > 1);
    addConstraint(cc);
}

void MaxSatSolver::corestrat_one_2(int64_t limit) {
    assert(decisionLevel() == 0);
    trace(maxsat, 10, "Use algorithm one-2");
    CardinalityConstraint cc;
    cc.bound = conflict.size() - 1;
    while(conflict.size() > 0) {
        weights[var(conflict.last())] -= limit;
        cc.lits.push(~conflict.last());
        conflict.pop();
    }
    assert(conflict.size() == 0);
    for(int i = 0; i < cc.bound; i++) {
        newVar();
        if(i != 0) addClause(~softLiterals.last(), mkLit(nVars()-1));
        setFrozen(nVars()-1, true);
        weights.push(limit);
        softLiterals.push(mkLit(nVars()-1));
        cc.lits.push(~mkLit(nVars()-1));
    }
    
    assert(cc.size() > 1);
    addEquality(cc);
}

void MaxSatSolver::corestrat_one_wc(int64_t limit) {
    assert(decisionLevel() == 0);
    trace(maxsat, 10, "Use algorithm one (wc)");
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
    
    assert(wc.size() > 1);
    addConstraint(wc);
}

void MaxSatSolver::corestrat_one_neg(int64_t limit) {
    assert(decisionLevel() == 0);
    trace(maxsat, 10, "Use algorithm one (neg)");
    CardinalityConstraint cc;
    cc.bound = conflict.size() - 1;
    while(conflict.size() > 0) {
        weights[var(conflict.last())] -= limit;
        cc.lits.push(~conflict.last());
        conflict.pop();
    }
    assert(conflict.size() == 0);
    for(int i = 0; i < cc.bound; i++) {
        newVar();
        if(i != 0) addClause(~softLiterals.last(), ~mkLit(nVars()-1));
        setFrozen(nVars()-1, true);
        weights.push(limit);
        softLiterals.push(~mkLit(nVars()-1));
        cc.lits.push(mkLit(nVars()-1));
    }
    
    assert(cc.size() > 1);
    addConstraint(cc);
}

void MaxSatSolver::corestrat_one_neg_wc(int64_t limit) {
    assert(decisionLevel() == 0);
    trace(maxsat, 10, "Use algorithm one (neg; wc)");
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
        if(i != 0) addClause(~softLiterals.last(), ~mkLit(nVars()-1));
        setFrozen(nVars()-1, true);
        weights.push(limit);
        softLiterals.push(~mkLit(nVars()-1));
        wc.lits.push(mkLit(nVars()-1));
        wc.coeffs.push(1);
    }
    
    assert(wc.size() > 1);
    addConstraint(wc);
}

void MaxSatSolver::corestrat_one_pmres(int64_t limit) {
    assert(decisionLevel() == 0);
    
    trace(maxsat, 10, "Use algorithm one-pmres");
    
    const int b = option_maxsat_tag * 2;
    const int m = ceil(2.0 * conflict.size() / (b-2.0));
    const int N = ceil(
            (
                conflict.size()         // literals in the core
                + conflict.size() - 1   // new soft literals
                + 2 * (m-1)             // new connectors
            ) / (m * 2.0)
        );
    // ceil((conflict.size() + m) / static_cast<double>(m));
    trace(maxsat, 15, "At most " << N*2 << " elements in " << m << " new constraints");

    Lit prec = lit_Undef;
    for(;;) {
        assert(conflict.size() > 0);
        
        CardinalityConstraint cc;
        
        int i = N;
        if(prec != lit_Undef) { cc.lits.push(prec); i--; }
        for(; i > 0; i--) {
            if(conflict.size() == 0) break;
            weights[var(conflict.last())] -= limit;
            cc.lits.push(~conflict.last());
            conflict.pop();
        }
        assert(cc.size() > 0);
        cc.bound = cc.size()-1;
        
        if(conflict.size() > 0) cc.bound++;

        for(i = 0; i < cc.bound; i++) {
            newVar();
            setFrozen(nVars()-1, true);
            cc.lits.push(~mkLit(nVars()-1));
            if(i != 0) addClause(~mkLit(nVars()-2), mkLit(nVars()-1)); // symmetry breaker
            if(i == 0 && conflict.size() > 0) { 
                weights.push(0);
                prec = mkLit(nVars()-1);
                
//                vec<Lit> lits;
//                for(int j = 0; j < cc.size()-1; j++) lits.push(~cc.lits[j]);
//                lits.push(prec);
//                addClause(lits);
            }
            else {
                weights.push(limit);
                softLiterals.push(mkLit(nVars()-1));
            }
        }
        
        trace(maxsat, 25, "Add constraint of size " << cc.size());
        addConstraint(cc);
        
        if(conflict.size() == 0) break;
    }
    
    assert(conflict.size() == 0);
}

void MaxSatSolver::corestrat_one_pmres_2(int64_t limit) {
    assert(decisionLevel() == 0);
    trace(maxsat, 10, "Use algorithm one-pmres-2");
    
    const int N = option_maxsat_tag;

    Lit prec = lit_Undef;
    vec<Lit> lits;
    for(;;) {
        assert(conflict.size() > 0);
        
        CardinalityConstraint cc;
        
        int i = N;
        if(prec != lit_Undef) { cc.lits.push(prec); lits.push(~prec); i--; }
        for(; i > 0; i--) {
            if(conflict.size() == 0) break;
            weights[var(conflict.last())] -= limit;
            cc.lits.push(~conflict.last());
            lits.push(conflict.last());
            conflict.pop();
        }
        assert(cc.size() > 0);
        cc.bound = cc.size()-1;
        
        for(i = 0; i < cc.bound; i++) {
            newVar();
            cc.lits.push(~mkLit(nVars()-1));
            if(i != 0) addClause(~mkLit(nVars()-2), mkLit(nVars()-1)); // symmetry breaker
            setFrozen(nVars()-1, true);
            weights.push(limit);
            softLiterals.push(mkLit(nVars()-1));
        }
        
        addConstraint(cc);
        
        if(conflict.size() == 0) break;

        newVar();
        weights.push(0);
        prec = mkLit(nVars()-1);
        for(int j = 0; j < lits.size(); j++) addClause(~lits[j], ~prec);
        lits.push(prec);
        addClause(lits);
        lits.clear();
    }
    
    assert(conflict.size() == 0);
}

void MaxSatSolver::corestrat_pmres(int64_t limit) {
    assert(decisionLevel() == 0);
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
            
            if(conflict.size() > 1) {
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
        }
        
        conflict.pop();
    }
    assert(conflict.size() == 0);
}

void MaxSatSolver::corestrat_pmres_reverse(int64_t limit) {
    assert(decisionLevel() == 0);
    trace(maxsat, 10, "Use algorithm pmres (reverse)");

    Lit prec = lit_Undef;
    softLiterals.growTo(softLiterals.size() + conflict.size() - 1);
    int pos = softLiterals.size();
    vec<Lit> lits;
    while(conflict.size() > 0) {
        weights[var(conflict.last())] -= limit;
        
        if(prec == lit_Undef) prec = ~conflict.last();
        else {
            // disjunction
            newVar();
            setFrozen(nVars()-1, true);
            weights.push(limit);
            softLiterals[--pos] = mkLit(nVars()-1);
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
        }
        
        conflict.pop();
    }
    assert(conflict.size() == 0);
}

void MaxSatSolver::corestrat_pmres_split_conj(int64_t limit) {
    assert(decisionLevel() == 0);
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

void MaxSatSolver::corestrat_pmreslog(int64_t limit) {
    assert(decisionLevel() == 0);
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

void MaxSatSolver::corestrat_kdyn(int64_t limit) {
    assert(decisionLevel() == 0);
    trace(maxsat, 10, "Use algorithm kdyn");
    
    const int b = conflict.size() <= 2 ? 8 : ceil(log10(conflict.size()) * 16);
    const int m = ceil(2.0 * conflict.size() / (b-2.0));
    const int N = ceil(
            (
                conflict.size()         // literals in the core
                + conflict.size() - 1   // new soft literals
                + 2 * (m-1)             // new connectors
            ) / (m * 2.0)
        );
    // ceil((conflict.size() + m) / static_cast<double>(m));
    trace(maxsat, 15, "At most " << N*2 << " elements in " << m << " new constraints");

    Lit prec = lit_Undef;
    for(;;) {
        assert(conflict.size() > 0);
        
        CardinalityConstraint cc;
        
        int i = N;
        if(prec != lit_Undef) { cc.lits.push(prec); i--; }
        for(; i > 0; i--) {
            if(conflict.size() == 0) break;
            weights[var(conflict.last())] -= limit;
            cc.lits.push(~conflict.last());
            conflict.pop();
        }
        assert(cc.size() > 0);
        cc.bound = cc.size()-1;
        
        if(conflict.size() > 0) cc.bound++;

        for(i = 0; i < cc.bound; i++) {
            newVar();
            setFrozen(nVars()-1, true);
            cc.lits.push(~mkLit(nVars()-1));
            if(i != 0) addClause(~mkLit(nVars()-2), mkLit(nVars()-1)); // symmetry breaker
            if(i == 0 && conflict.size() > 0) { 
                weights.push(0);
                prec = mkLit(nVars()-1);
            }
            else {
                weights.push(limit);
                softLiterals.push(mkLit(nVars()-1));
            }
        }
        
        trace(maxsat, 25, "Add constraint of size " << cc.size());
        addConstraint(cc);
        
        if(conflict.size() == 0) break;
    }
    
    assert(conflict.size() == 0);
}

} // namespace aspino

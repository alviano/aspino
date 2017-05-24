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

#include "PMaxSatSolver.h"

#include <mtl/Map.h>
#include <utils/System.h>
#include <core/Dimacs.h>

#include "utils/algorithm.h"

using Glucose::Map;

extern Glucose::BoolOption option_print_model;

Glucose::IntOption option_threads("MAIN", "t", "Number of threads (used by PMaxSat).", 4, Glucose::IntRange(0, INT32_MAX));

namespace aspino {

int PMaxSatSolver::Task::currLevel = 0;
    
static const int buffer_size = 1048576;

static std::mutex trace_lock;

#define trace_(id, level, msg) { \
    std::unique_lock<std::mutex> locker(trace_lock); \
    trace(maxsat, level, "[solver " << id << "] " << msg); \
    locker.unlock(); \
}


class LineStream {
private:    
    int           pos;

public:
    char buf[buffer_size];
    explicit LineStream() : pos(0) {}

    bool readline() { 
        pos = 0;
        buf[0] = '\0';
        cin.getline(buf, sizeof(buf));
        return !cin.eof();
    }
    int  operator *  () const { return buf[pos]; }
    void operator ++ ()       { pos++; }
    int  position    () const { return pos; }
};    

int64_t parseLong(LineStream& in) {
    int64_t val = 0;
    bool    neg = false;
    Glucose::skipWhitespace(in);
    if      (*in == '-') neg = true, ++in;
    else if (*in == '+') ++in;
    if (*in < '0' || *in > '9') fprintf(stderr, "PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    while (*in >= '0' && *in <= '9')
        val = val*10 + (*in - '0'),
        ++in;
    return neg ? -val : val; 
}
    
void solveTask(void* solver_, void* task_) {
    PMaxSatSolver* solver = static_cast<PMaxSatSolver*>(solver_);
    solver->cancelUntil(0);
    
    PMaxSatSolver::Msg msg(*static_cast<PMaxSatSolver::Task*>(task_));
    msg.assumptions.copyTo(solver->assumptions);
    
    msg.status = solver->solveTask();
    if(msg.status == l_True) {
        solver->updateUpperBound();
        msg.upperbound = solver->upperbound;
    }
    else if(msg.status == l_False) {
        solver->cancelUntil(0);
        solver->trim();
        solver->conflict.moveTo(msg.core);
    }
    
    static_cast<PMaxSatSolver::Task*>(task_)->done_ = true;
    solver->add(msg);
}
    
PMaxSatSolver::PMaxSatSolver(PMaxSatSolver* master, int id_) : id(id_), taskLevel(-1), taskId(-1), myupperbound(INT64_MAX), upperbound(INT64_MAX), lowerbound(0), tSolve(NULL), lastConflict(0), lastPropagation(0), lastCallCpuTime(0.0) /*, timeBudget(0.0)*/ {
    setIncrementalMode();
    if(master != NULL) solvers.push(master);
    else {
        solvers.push(this);
        for(int i = 1; i < option_threads; i++) solvers.push(new PMaxSatSolver(this, i));
    }
}

PMaxSatSolver::~PMaxSatSolver() {
    for(int i = 1; i < solvers.size(); i++) delete solvers[i];
}

void PMaxSatSolver::interrupt() { 
    this->exit(1);
}

void PMaxSatSolver::sameSoftVar(Lit soft, int64_t weight) {
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

void PMaxSatSolver::addWeightedClause(vec<Lit>& lits, int64_t weight) {
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
    
void PMaxSatSolver::parse(gzFile in_) {
    Glucose::StreamBuffer in(in_);

    upperbound = INT64_MAX;

    bool weighted = false;
    int64_t top = -1;
    int64_t weight = 1;
    
    vec<Lit> lits;
    int vars = 0;
    int count = 0;
    for(;;) {
        Glucose::skipWhitespace(in);
        if(*in == EOF) break;
        if(*in == 'p') {
            ++in;
            if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            ++in;
            if(*in == 'w') { weighted = true; ++in; }
            
            if(Glucose::eagerMatch(in, "cnf")) {
                vars = Glucose::parseInt(in);
                inClauses = Glucose::parseInt(in);
                if(weighted && *in != '\n') top = parseLong(in);
                
                nInVars(vars);
                while(weights.size() < nInVars()) weights.push(0);
                while(nVars() < nInVars()) newVar();
            }
            else {
                cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            }
        }
        else if(*in == 'c') continue;
        else {
            count++;
            if(weighted) weight = parseLong(in);
            readClause(in, *this, lits);
            if(weight == top) addClause_(lits);
            else addWeightedClause(lits, weight);
        }
    }
    
//    if(count != inClauses) cerr << "WARNING! DIMACS header mismatch: wrong number of clauses." << endl, exit(3);
    
    inClauses = clauses.size();
}

void PMaxSatSolver::hardening() {
    cancelUntil(0);
    int j = 0;
    for(int i = 0; i < softLiterals.size(); i++) {
        int64_t w = weights[var(softLiterals[i])];
        if(w == 0) continue;
        if(w + lowerbound >= upperbound) {
            addClause(softLiterals[i]);
            trace_(id, 30, "Hardening of " << softLiterals[i] << " of weight " << w);
            weights[var(softLiterals[i])] = 0;
            continue;
        }
        softLiterals[j++] = softLiterals[i];
    }
    softLiterals.shrink_(softLiterals.size()-j);
}

void PMaxSatSolver::setAssumptions(int64_t limit) {
    assumptions.clear();
    cancelUntil(0);
    int j = 0;
    for(int i = 0; i < softLiterals.size(); i++) {
        int64_t w = weights[var(softLiterals[i])];
        if(w == 0) continue;
//        assert(w + lowerbound < upperbound);
        softLiterals[j++] = softLiterals[i];
        if(w >= limit) assumptions.push(softLiterals[i]);
    }
    softLiterals.shrink_(softLiterals.size()-j);
}

int64_t PMaxSatSolver::computeNextLimit(int64_t limit) const {
    int64_t next = limit;
    for(int i = 0; i < softLiterals.size(); i++) {
        int64_t w = weights[var(softLiterals[i])];
        if(w == 0) continue;
        if(w >= limit) continue;
        if(next == limit || w > next) next = w;
    }
    return next;
}

bool PMaxSatSolver::tautology(Var v, CRef c1, CRef c2) {
    Clause& a = ca[c1];
    Clause& b = ca[c2];
    for(int i = 0; i < a.size(); i++) {
        if(var(a[i]) == v) continue;
        for(int j = 0; j < b.size(); j++) {
            if(var(b[j]) == v) continue;
            if(var(a[i]) == var(b[j]) && sign(a[i]) != sign(b[j])) { return true; }
        }
    }
    return false;    
}

void PMaxSatSolver::removeSoftLiteralsAtLevelZero() {
    assert(decisionLevel() == 0);
    int j = 0;
    for(int i = 0; i < softLiterals.size(); i++) {
        if(value(softLiterals[i]) != l_Undef) { 
            if(value(softLiterals[i]) == l_False) {
                trace_(id, 20, "False soft literal at level 0: " << softLiterals[i] << " of weight " << weights[var(softLiterals[i])]);
                updateLowerBound(weights[var(softLiterals[i])]);
            }
            weights[var(softLiterals[i])] = 0;
            continue;
        }
        softLiterals[j++] = softLiterals[i];
    }
    softLiterals.shrink_(softLiterals.size() - j);
}

void PMaxSatSolver::lastCheck() {
    trace_(id, 2, "Bounds: [" << lowerbound << ":" << upperbound << "]");

    if(upperbound == INT64_MAX) { cout << "s UNSATISFIABLE" << endl; return; }
        
    if(lowerbound < upperbound) {
        if(!ok) status = l_False;
        else {
            trace_(id, 5, "Last check");
            setAssumptions(1);
//            cout << assumptions << endl;
//            softLiterals.moveTo(assumptions);
//            for(int i = 0; i < weights.size(); i++) if(weights[i] < 0) cout << i << " "<< weights[i] << endl;
            PseudoBooleanSolver::solve();
            trace_(id, 5, (status == l_True ? "SAT!" : "UNSAT!"));
        }
        if(status == l_True) updateUpperBound();
        else lowerbound = upperbound;
    }
    trace_(id, 2, "Bounds: [" << lowerbound << ":" << upperbound << "]");
    
    assert(lowerbound == upperbound);
    assert(upperbound < INT64_MAX);
    
    cout << "o " << lowerbound << endl;
    cout << "s OPTIMUM FOUND" << endl;
    for(int i = 0; i < solvers.size(); i++) if(solvers[i]->handleModel()) break;
}

void PMaxSatSolver::quickSort(int left, int right) {
    int i = left, j = right;
    Lit tmp;
    assert((left + right) / 2 >= 0);
    assert_msg((left + right) / 2 < softLiterals.size(), "Accessing element " << (left + right) / 2 << " in array of size " << softLiterals.size());
    double pivot = activity[var(softLiterals[(left + right) / 2])];

    /* partition */
    while (i <= j) {
        while (activity[var(softLiterals[i])] < pivot)
              i++;
        while (activity[var(softLiterals[j])] > pivot)
              j--;
        if (i <= j) {
              tmp = softLiterals[i];
              softLiterals[i] = softLiterals[j];
              softLiterals[j] = tmp;
              i++;
              j--;
        }
    };

    /* recursion */
    if (left < j)
        quickSort(left, j);
    if (i < right)
        quickSort(i, right);
}

int64_t PMaxSatSolver::computeConflictWeight() const {
    int64_t min = INT64_MAX;
    for(int i = 0; i < conflict.size(); i++) if(weights[var(conflict[i])] < min) min = weights[var(conflict[i])];
    return min;
}

void PMaxSatSolver::updateLowerBound(int64_t limit) {
    lowerbound += limit;
    cout << "o " << lowerbound << endl;
}

void PMaxSatSolver::updateUpperBound() {
    assert(status == l_True);
    int64_t newupperbound = lowerbound;
    for(int i = 0; i < softLiterals.size(); i++) {
        if(value(softLiterals[i]) == l_False) newupperbound += weights[var(softLiterals[i])];
    }
    if(newupperbound < upperbound) {
        upperbound = newupperbound;
        myupperbound = newupperbound;
        copyModel();
        trace_(id, 200, "Model: " << model);
    }
}

void PMaxSatSolver::trim() {
    assert(decisionLevel() == 0);
    
    if(conflict.size() <= 1) return;

    int counter = 0;

    do{
        counter++;
        assumptions.clear();
        for(int i = 0; i < conflict.size(); i++) assumptions.push(~conflict[i]);
        PseudoBooleanSolver::solve();
        if(status == l_Undef) { 
            assert(conflict.size() == 0); 
            for(int i = 0; i < assumptions.size(); i++) conflict.push(~assumptions[i]);
        }
        trace_(id, 15, "Trim " << assumptions.size() - conflict.size() << " literals from conflict");
        trace_(id, 100, "Conflict: " << conflict);
        cancelUntil(0);
        if(conflict.size() <= 1) return;
    }while(assumptions.size() > conflict.size());
    
    if(counter % 2 == 1) for(int i = 0; i < assumptions.size(); i++) conflict[i] = ~assumptions[i];
    
    assert(conflict.size() > 1);
}

void PMaxSatSolver::processCore(int64_t limit) {
    assert(decisionLevel() == 0);
    trace_(id, 10, "Use algorithm one");
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
    addConstraint(cc);
}

bool PMaxSatSolver::handleModel() const {
    if(!option_print_model) return true;
    if(myupperbound != upperbound) return false;
    printModel();
    return true;
}

void PMaxSatSolver::stop(Task& task) {
    assert(id == 0);
    if(task.tSolve == NULL) return;
    assert(task.solverId != -1);
    trace_(id, 20, "Stop solver " << task.solverId);

    PMaxSatSolver* solver = solvers[task.solverId];
    solver->asynch_interrupt = true;
    
    assert(task.tSolve->joinable());
    task.tSolve->join();
    delete task.tSolve;
    task.tSolve = NULL;
    assert(task.done());

    solver->asynch_interrupt = false;
    freeSolvers.push(task.solverId);
    trace_(id, 1, "FREE " << freeSolvers);
}

void PMaxSatSolver::assignTasks() {
    trace_(id, 1, "FREE SOLVERs: " << freeSolvers);
    for(int j = 0; j < tasks.size(); j++) {
        assert(tasks[j].assumptions.size() != 0);
        if(freeSolvers.size() == 0) break;
        if(tasks[j].assigned()) continue;
        
        PMaxSatSolver* solver = solvers[freeSolvers.last()];
        freeSolvers.pop();
        tasks[j].solverId = solver->id;
        tasks[j].tSolve = new thread(aspino::solveTask, static_cast<void*>(solver), static_cast<void*>(&tasks[j]));
        
        trace_(id, 20, "Assign task " << tasks[j] << " to solver " << solver->id);
    }
}

void PMaxSatSolver::add(const Msg& msg) {
    if(this != solvers[0]) { solvers[0]->add(msg); return; }
    
    trace_(id, 1, "ADD " << msg);
    
    std::unique_lock<std::mutex> locker(lock);
    msgs.push_back(msg);
    locker.unlock();
    noMsgs.notify_one();
}

lbool PMaxSatSolver::solveTask() { 
    time_t t = time(0); 
    setConfBudget(solvers[0]->budget);
    cout << "START ASS: " << assumptions.size() << endl;
    lbool res = PseudoBooleanSolver::solve();
    cout << "TIME: " << difftime(time(0), t) << "; ASS: " << assumptions.size() << endl;
    budgetOff();
    clearInterrupt();
    return res; 
}

lbool PMaxSatSolver::solve() {
    assert_msg(weights.size() == nVars(), weights.size() << " " << nVars());

    PseudoBooleanSolver::solve();
    if(status == l_False) { cout << "s UNSATISFIABLE" << endl; return l_False; }
    if(status == l_True) updateUpperBound();
    cout << "c " << upperbound << " ub" << endl;
    cancelUntil(0);

    removeSoftLiteralsAtLevelZero();
    
    assert(cconstraints.size() == 0);
    for(int i = 1; i < solvers.size(); i++) {
        solvers[i]->nInVars(nInVars());
        solvers[i]->lowerbound = lowerbound;
        solvers[i]->upperbound = upperbound;
        
        while(solvers[i]->nVars() < nVars()) {
            solvers[i]->newVar();
            solvers[i]->weights.push(0);
            if(value(solvers[i]->nVars()-1) != l_Undef) solvers[i]->addClause(mkLit(solvers[i]->nVars()-1, value(solvers[i]->nVars()-1) == l_False));
        }
        
        vec<Lit> lits;
        for(int j = 0; j < clauses.size(); j++) {
            Clause& clause = ca[clauses[j]];
            for(int k = 0; k < clause.size(); k++) lits.push(clause[k]);
            solvers[i]->addClause_(lits);
            lits.clear();
        }
        
        for(int j = 0; j < softLiterals.size(); j++) {
            lits.push(softLiterals[j]);
            solvers[i]->addWeightedClause(lits, weights[var(softLiterals[j])]);
            lits.clear();
        }
    }

    vec<Lit> soft;
    int softSatUpTo = 0;
    
    vec<Lit> core;
    int coreSatUpTo = 0;

    int64_t limit = computeNextLimit(INT64_MAX);
    int64_t nextLimit;
    
    budget = conflicts;
    
    freeSolvers.clear();
    for(int i = 0; i < solvers.size(); i++) freeSolvers.push(i);

    for(;;) {
        assert(tasks.size() == 0);
//                if(softLiterals.size() > 1) quickSort(0, softLiterals.size()-1);
            
        Task::currLevel += 1;
        for(int i = 0; i < solvers.size(); i++) solvers[i]->hardening();
        if(lowerbound == upperbound) { lastCheck(); return l_True; }
        setAssumptions(limit);
        
        lastConflict = conflicts;
        lastPropagation = propagations;
        lastCallCpuTime = Glucose::cpuTime();

        trace_(id, 2, "Solve with " << assumptions.size() << " assumptions. Current bounds: [" << lowerbound << ":" << upperbound << "]. Difference: " << (upperbound - lowerbound) << ". At least " << ((upperbound-lowerbound) / limit) << " cores to terminate.");
        trace_(id, 100, "Assumptions: " << assumptions);
        
        assert(decisionLevel() == 0);
        if(assumptions.size() == 0 && upperbound != INT64_MAX) status = l_Undef;
        else {
            vec<Lit> lits;
            if(core.size() == 0) {
                int j = 0;
                for(int i = 0; i < softSatUpTo && j < assumptions.size(); i++) if(assumptions[j] == soft[i]) j++;
                softSatUpTo = j;
                coreSatUpTo = softSatUpTo;
                assumptions.moveTo(soft);
                cout << "SOFTSATUPTO " << softSatUpTo << "/" << soft.size() << endl;
                soft.copyTo(lits);
            }
            else {
                for(int i = core.size() - 1; i >= 0; i--) lits.push(~core[i]);
                cout << "CORESATUPTO " << coreSatUpTo << "/" << core.size() << endl;
            }
            
            int m = coreSatUpTo - 1;
            int pr = 1;
            while (m + pr < lits.size() - (core.size() == 0 ? 0 : 1)) {
                tasks.push();
                
                int limit = m + pr + 1;
                for(int j = 0; j < limit; j++) tasks.last().assumptions.push(lits[j]);
                
                if(m + 2 * pr >= lits.size() - (core.size() == 0 ? 0 : 1)) {
                    m += pr;
                    pr = 1;
                } else {
                    pr *= 2;
                }
            }
            
            assignTasks();
        }
        
        bool debug = false;
        int NEXT = tasks[0].assumptions.size();
        
        // handle msgs
        int64_t ub = upperbound;
        for(;;) {
            std::unique_lock<std::mutex> locker(lock);
            
            if(msgs.size() == 0 || debug) {
                bool stop_ = msgs.size() == 0;
                for(int i = 0; i < tasks.size(); i++) if(!tasks[i].done()) { stop_ = false; break; }
                if(stop_) {
                    locker.unlock();
                    break;
                }
                if(debug) noMsgs.wait(locker, [this, NEXT](){ assignTasks(); for(list<Msg>::const_iterator it = msgs.begin(); it != msgs.end(); it++) if(it->assumptions.size() == NEXT) return true; return false; });
                else noMsgs.wait(locker, [this](){ return msgs.size() > 0; });
            }
            
            trace_(id, 5, "Pop from " << msgs.size() << " msgs");
            for(list<Msg>::const_iterator it = msgs.begin(); it != msgs.end(); it++) trace_(id, 5, " >>> " << *it);
            Msg* msg_ = NULL;
            if(debug) {
                list<Msg>::iterator it;
                for(it = msgs.begin(); it != msgs.end(); it++) if(it->assumptions.size() == NEXT) {
                    break;
                }
                assert(it != msgs.end());
                msg_ = new Msg(*it);
                msgs.erase(it);
                NEXT++;
            }
            else {
                msg_ = new Msg(msgs.front());
                msgs.pop_front();
            }
            Msg& msg = *msg_;
            for(list<Msg>::const_iterator it = msgs.begin(); it != msgs.end(); it++) trace_(id, 5, "after pop " << *it);
            locker.unlock();
            
            trace_(id, 10, msg);
            
            assert(msg.level == Task::currLevel);
//            if(msg.level < Task::currLevel) { trace_(id, 20, "ignored level"); continue;} 
            
            if(msg.status == l_True) {
                if(ub > msg.upperbound) { ub = msg.upperbound; assert_msg(ub >= lowerbound, msg); }
                if(msg.assumptions.size() < (core.size() == 0 ? softSatUpTo : coreSatUpTo)) { trace_(id, 20, "ignored sat"); continue; }
//                if(msg.id < Task::minIdSat) { trace_(id, 20, "ignored sat"); continue;} 
//                Task::minIdSat = msg.id;
                if(core.size() == 0) {
                    assert_msg(softSatUpTo < msg.assumptions.size(), softSatUpTo << " " << msg.assumptions.size());
                    softSatUpTo = msg.assumptions.size();
                }
                else {
                    assert(coreSatUpTo < msg.assumptions.size());
                    coreSatUpTo = msg.assumptions.size();
                }

                for(int i = 0; i < tasks.size(); i++) {
                    // TODO ottimizzare
                    if(tasks[i].assumptions.size() <= msg.assumptions.size()) { stop(tasks[i]); assert(!solvers[tasks[i].solverId]->asynch_interrupt); }
                }
                
                assignTasks();
                if(freeSolvers.size() > 0 && tasks.size() > solvers.size()) for(int i = 0; i < tasks.size(); i++) stop(tasks[i]);
            }
            else if(msg.status == l_False) {
                if(core.size() != 0 && msg.core.size() >= core.size()) { trace_(id, 10, "ignored unsat"); continue;}
                    
                for(int i = 0; i < tasks.size(); i++) stop(tasks[i]);
                tasks.clear();
                assert(freeSolvers.size() == solvers.size());
                
                if(core.size() == 0) {
                    assert(msg.core.size() > 0); //if(msg.core.size() == 0) return l_True; // TERMINA (potrebbe non terminare se i thread non sono demoni)
                    msg.core.moveTo(core);
                    
                    int j = core.size() - 1;
                    for(int i = 0; i < softSatUpTo && j >= 0; i++) if(~core[j] == soft[i]) j--;
                    coreSatUpTo = core.size() - 1 - j;
                }
                else {
                    assert(msg.core.size() > 0);
                    
                    int j = 0;
                    for(int i = 0; i < coreSatUpTo && j < msg.core.size(); i++) if(msg.core[j] == core[i]) j++;
                    coreSatUpTo = j;
                    
                    msg.core.moveTo(core);
                }
            }
            else {
                for(int i = 0; i < tasks.size(); i++) {
                    if(tasks[i].assumptions.size() == msg.assumptions.size()) { stop(tasks[i]); assert(!solvers[tasks[i].solverId]->asynch_interrupt); break; }
                }
                assignTasks();
                if(freeSolvers.size() > 0 && tasks.size() > solvers.size()) for(int i = 0; i < tasks.size(); i++) stop(tasks[i]);
            }
        } // msgs
        assert(msgs.size() == 0);
        
        for(int i = 0; i < tasks.size(); i++) assert(tasks[i].done());
        tasks.clear();
        assert(freeSolvers.size() == solvers.size());
        
        if(ub < upperbound) {
            for(int i = 0; i < solvers.size(); i++) solvers[i]->upperbound = ub;
            cout << "c " << upperbound << " ub" << endl;
        }
        
        if(core.size() > 0) {
            if(coreSatUpTo < core.size() - 1) { trace_(id, 10, "ignore core " << core); continue; }
            
            trace_(id, 10, "Process core " << core);
            // TODO update softSatUpTo
//            softSatUpTo = 0;
            
            for(int i = 0; i < solvers.size(); i++) {
                core.copyTo(solvers[i]->conflict);
                solvers[i]->cancelUntil(0);
                solvers[i]->lowerbound += limit;
                solvers[i]->processCore(limit);
            }
            core.clear();
            
            cout << "o " << lowerbound << endl;
        }
        else if(softSatUpTo == soft.size()) {
            nextLimit = computeNextLimit(limit);
            if(nextLimit == limit) {
                trace_(id, 4, (status == l_True ? "SAT!" : "Skip!") << " No other limit to try");
                lastCheck();
                return l_True;
            }
            
            trace_(id, 4, (status == l_True ? "SAT!" : "Skip!") << " Decrease limit to " << nextLimit);
            limit = nextLimit;
        }
        else {
            budget *= 2;
            cout << "budget " << budget << endl;
        }
        
    }  // primo for
    
    return l_True;
}

} // namespace aspino

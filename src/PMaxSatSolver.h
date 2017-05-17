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

#ifndef __PMaxSatSolver_h__
#define __PMaxSatSolver_h__

#include "PseudoBooleanSolver.h"

#include <thread>
#include <mutex>
#include <condition_variable>
#include <list>

namespace aspino {

class PMaxSatSolver : public PseudoBooleanSolver {
public:
    PMaxSatSolver(PMaxSatSolver* master = NULL, int id = 0);
    virtual ~PMaxSatSolver();
    
    virtual void interrupt();
    
    virtual void parse(gzFile in);
    virtual bool okay() const { return true; }
    
    void addWeightedClause(vec<Lit>& lits, int64_t weight);
    
    virtual lbool solve();
    virtual lbool solve(int) { return solve(); }
    
//    void onTick();
    
    std::mutex lock;
    std::condition_variable noMsgs;
    
private:
    vec<PMaxSatSolver*> solvers;
    vec<int> freeSolvers;
    
    int id;
    int taskLevel;
    int taskId;
    int64_t myupperbound;
    int64_t upperbound;
    int64_t lowerbound;
    
    vec<Lit> softLiterals;
    vec<int64_t> weights;
    
    struct Task {
        friend ostream& operator<<(ostream& out, const Task& task) { return out << "[level=" << task.level << "; id=" << task.id << "; assumptions=" << task.assumptions << "]"; }
        inline Task() : level(currLevel), id(nextId++), tSolve(NULL), solverId(-1) {}
//        inline Task(int level_, int id_, vec<Lit>& ass) :level(level_), id(id_) { ass.moveTo(assumptions); }
        int level;
        int id;
        vec<Lit> assumptions;
        thread* tSolve;
        int solverId;
        
        inline bool assigned() const { return solverId != -1; }
        inline bool done() const { return assigned() && tSolve == NULL; }
        
        static int currLevel;
        static int nextId;
        static int minIdSat;
    };
    vec<Task> tasks;
    
    void stop(Task& task);
    void assignTasks();
    inline lbool solveTask() { return PseudoBooleanSolver::solve(); }
    
    struct Msg {
        friend ostream& operator<<(ostream& out, const Msg& msg) { return out << "[level=" << msg.level << "; id=" << msg.id << "; assumptions=" << msg.assumptions << "; status=" << msg.status << "; upperbound=" << msg.upperbound << "; core=" << msg.core << "]"; }
        inline Msg() : level(-1), id(-1), status(l_Undef), upperbound(INT64_MAX) {}
        inline Msg(const Task& task) : level(task.level), id(task.id), status(l_Undef), upperbound(INT64_MAX) { task.assumptions.copyTo(assumptions); }
        inline Msg(const Msg& init) : level(init.level), id(init.id), status(init.status), upperbound(init.upperbound) { init.assumptions.copyTo(assumptions); init.core.copyTo(core); }
        int level;
        int id;
        vec<Lit> assumptions;
        lbool status;
        vec<Lit> core;
        int64_t upperbound;
    };
    list<Msg> msgs;
    void add(const Msg& msg);
    
    friend void solveTask(void* solver_, void* task_);
    
    thread* tSolve;
    
    uint64_t lastConflict;
    uint64_t lastPropagation;
    double lastCallCpuTime;
    double timeBudget;
    
    void lastCheck();

    void quickSort(int left, int right);
    
    void sameSoftVar(Lit soft, int64_t weight);
    
    void removeSoftLiteralsAtLevelZero();
    bool tautology(Var var, CRef c1, CRef c2);
    void preprocess();
    void hardening();
    void setAssumptions(int64_t limit);
    int64_t computeNextLimit(int64_t limit) const;
    void solveCurrentLevel();
    void solve_();
    
    int64_t computeConflictWeight() const;
    
    void updateLowerBound(int64_t limit);
    void updateUpperBound();
    
    void trim();
    void shrink(int64_t weight);
    void processCore(int64_t limit);
    
    bool handleModel() const;
    
    inline void printTask() const { if(taskLevel >= 0 || true) cout << "l " << taskLevel << " " << taskId << " "; }
};



} // namespace aspino

#endif
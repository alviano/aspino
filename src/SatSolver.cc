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

#include "SatSolver.h"

#include "utils/algorithm.h"
#include "utils/ParseUtils.h"
#include "core/Dimacs.h"

#include <iostream>
#include <cmath>

using namespace std;

namespace aspino {

void SatSolver::parse(gzFile in_) {
    Glucose::StreamBuffer in(in_);

    vec<Lit> lits;
    int vars = 0;
    int clauses = 0;
    int count = 0;
    for(;;) {
        skipWhitespace(in);
        if(*in == EOF) break;
        if(*in == 'p') {
            if(eagerMatch(in, "p cnf")) {
                vars = parseInt(in);
                clauses = parseInt(in);
            }
            else {
                cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            }
        }
        else if(*in == 'c')
            skipLine(in);
        else {
            count++;
            readClause(in, *this, lits);
            addClause(lits); 
        }
    }
    if(nVars() != vars)
        cerr << "WARNING! DIMACS header mismatch: wrong number of variables." << endl;
    nInVars(min(vars, nVars()));
    if(count != clauses)
        cerr << "WARNING! DIMACS header mismatch: wrong number of clauses." << endl;
    inVars = nVars();
}

lbool SatSolver::solve() {
    conflict.clear();
    cancelUntil(0);

    solves++;

    solve_();

    return status;
}

lbool SatSolver::solve(int n) {
    conflict.clear();
    cancelUntil(0);

    solves++;

    lbool ret = l_False;
    int count = 0;
    for(;;) {
        solve_();
        if(status == l_False) {
            if(conflict.size() == 0) ok = false;
            break;
        }
        assert(status == l_True);
        if(++count == 1) printStatus();
        
        cout << "c Model " << count << endl;
        copyModel();
        printModel();
        ret = l_True;
        if(--n == 0) break;
        if(decisionLevel() == 0) break;
        learnClauseFromModel();
    }
    if(ret == l_False) printStatus();
    return ret;
}

void SatSolver::quickSort(int left, int right) {
    int i = left, j = right;
    Lit tmp;
    assert((left + right) / 2 >= 0);
    assert_msg((left + right) / 2 < assumptions.size(), "Accessing element " << (left + right) / 2 << " in array of size " << assumptions.size());
    double pivot = activity[var(assumptions[(left + right) / 2])];

    /* partition */
    while (i <= j) {
        while (activity[var(assumptions[i])] < pivot)
              i++;
        while (activity[var(assumptions[j])] > pivot)
              j--;
        if (i <= j) {
              tmp = assumptions[i];
              assumptions[i] = assumptions[j];
              assumptions[j] = tmp;
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

void SatSolver::solve_() {
    if(!ok) { status = l_False; return; }

    status = l_Undef;

    // Search:
    int curr_restarts = 0;
    while(status == l_Undef) {
//        cancelUntil(0);
//        shuffle(assumptions);
//        if(sortAssumptions && assumptions.size() > 1 && curr_restarts % 100 == 0) { cancelUntil(0); quickSort(0, assumptions.size()-1); }
//        if(curr_restarts == 1000)break;
        status = search(0);
        if(!withinBudget()) break;
        curr_restarts++;
//        trace(maxsat, 1, "Restart: " << curr_restarts);
    }

    if(status == l_False && conflict.size() == 0)
        ok = false;
}

// This is almost function search from glucose
lbool SatSolver::search(int)
{
    assert(ok);
    int         backtrack_level;
    int         conflictC = 0;
    vec<Lit>    learnt_clause,selectors;
    unsigned int nblevels,szWoutSelectors;
    bool blocked=false;
    starts++;
    for (;;){
        CRef confl;
        do {
            confl = propagate();
            if(confl != CRef_Undef) break;
            confl = morePropagate();
            if(confl != CRef_Undef) break;
        }while(nextToPropagateByUnit() < trail.size());

        if (confl != CRef_Undef){
            // CONFLICT
	  conflicts++; conflictC++;conflictsRestarts++;
	  if(conflicts%5000==0 && var_decay<0.95)
            var_decay += 0.01;

	  if (decisionLevel() == 0) {
	    return l_False;
	    
	  }
	  
	  trailQueue.push(trail.size());
	  // BLOCK RESTART (CP 2012 paper)
	  if( conflictsRestarts>LOWER_BOUND_FOR_BLOCKING_RESTART && lbdQueue.isvalid()  && trail.size()>R*trailQueue.getavg()) {
	    lbdQueue.fastclear();
	    nbstopsrestarts++;
	    if(!blocked) {lastblockatrestart=starts;nbstopsrestartssame++;blocked=true;}
	  }

            learnt_clause.clear();
	    selectors.clear();
            analyze(confl, learnt_clause, selectors,backtrack_level,nblevels,szWoutSelectors);

	    lbdQueue.push(nblevels);
	    sumLBD += nblevels;
 

            cancelUntil(backtrack_level);

//            if (certifiedUNSAT) {
//              for (int i = 0; i < learnt_clause.size(); i++)
//                fprintf(certifiedOutput, "%i " , (var(learnt_clause[i]) + 1) *
//                            (-2 * sign(learnt_clause[i]) + 1) );
//              fprintf(certifiedOutput, "0\n");
//            }

            if (learnt_clause.size() == 1){
	      uncheckedEnqueue(learnt_clause[0]);nbUn++;
            }else{
                CRef cr = ca.alloc(learnt_clause, true);
		ca[cr].setLBD(nblevels); 
		ca[cr].setSizeWithoutSelectors(szWoutSelectors);
		if(nblevels<=2) nbDL2++; // stats
		if(ca[cr].size()==2) nbBin++; // stats
                learnts.push(cr);
                attachClause(cr);

                claBumpActivity(ca[cr]);
                uncheckedEnqueue(learnt_clause[0], cr);
            }
            varDecayActivity();
            claDecayActivity();

           
        }else{
	  // Our dynamic restart, see the SAT09 competition compagnion paper 
	  if ( asynch_interrupt || // ADDED BY MALVI
	      ( lbdQueue.isvalid() && ((lbdQueue.getavg()*K) > (sumLBD / conflictsRestarts)))) {
	    lbdQueue.fastclear();
	    progress_estimate = progressEstimate();
	    int bt = 0;
	    if(incremental) { // DO NOT BACKTRACK UNTIL 0.. USELESS
	      bt = (decisionLevel()<assumptions.size()) ? decisionLevel() : assumptions.size();
	    }
	    cancelUntil(bt);
	    return l_Undef; }


           // Simplify the set of problem clauses:
	  if (decisionLevel() == 0 && !simplify()) {
	    return l_False;
	  }
	    // Perform clause database reduction !
	    if(static_cast<unsigned>(conflicts) >= curRestart* nbclausesbeforereduce) 
	      {
	
		assert(learnts.size()>0);
		curRestart = (conflicts/ nbclausesbeforereduce)+1;
		reduceDB();
		nbclausesbeforereduce += incReduceDB;
	      }
	    
            Lit next = lit_Undef;
            while (decisionLevel() < assumptions.size()){
                // Perform user provided assumption:
                Lit p = assumptions[decisionLevel()];
                if (value(p) == l_True){
                    // Dummy decision level:
                    newDecisionLevel();
                }else if (value(p) == l_False){
                    analyzeFinal(~p, conflict);
                    return l_False;
                }else{
                    next = p;
                    break;
                }
            }

            if (next == lit_Undef){
                // New variable decision:
                decisions++;
                next = pickBranchLit();

                if (next == lit_Undef){
		  // Model found:
		  return l_True;
		}
            }

            // Increase decision level and enqueue 'next'
            newDecisionLevel();
            uncheckedEnqueue(next);
        }
    }
}

void SatSolver::printStatus() const {
    cout << "s " << (status == l_True ? "SATISFIABLE" : status == l_False ? "UNSATISFIABLE" : "INDETERMINATE") << endl;
}

void SatSolver::copyModel() {
    // Extend & copy model:
    model.growTo(nVars());
    for (int i = 0; i < nVars(); i++) model[i] = value(i);
    Glucose::SimpSolver::extendModel();
}

void SatSolver::printModel() const {
    assert(model.size() >= nInVars());
    cout << "v";
    for(int i = 0; i < nInVars(); i++)
        cout << " " << (model[i] == l_False ? "-" : "") << (i+1);
    cout << endl;
}

void SatSolver::learnClauseFromModel() {
    vec<Lit> lits;
    lits.growTo(trail_lim.size());
    for(int i = 0; i < lits.size(); i++) {
        Lit lit = trail[trail_lim[i]];
        assert(reason(var(lit)) == CRef_Undef);
        assert(level(var(lit)) > 0);
        lits[lits.size() - level(var(lit))] = ~lit;
    }
    trace(sat, 2, "Adding clause from model: " << lits);
    cancelUntil(decisionLevel()-1);
    if (lits.size() == 1)
        uncheckedEnqueue(lits[0]);
    else {
        CRef cr = ca.alloc(lits, true);
        clauses.push(cr);
        attachClause(cr);
//        claBumpActivity(ca[cr]);
        uncheckedEnqueue(lits[0], cr);
    }
}

double SatSolver::luby(double y, int x) {
    int size, seq;
    for (size = 1, seq = 0; size < x+1; seq++, size = 2*size+1);

    while (size-1 != x){
        size = (size-1)>>1;
        seq--;
        x = x % size;
    }

    return pow(y, seq);
}

} // namespace aspino
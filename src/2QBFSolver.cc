/*
 *  Copyright (C) 2015  Mario Alviano (mario@alviano.net)
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

#include "2QBFSolver.h"

#include <core/Dimacs.h>

Glucose::BoolOption option_qbf_printmodel("QBF", "qbf-print-model", "Print (counter-)model if found.", true);

namespace aspino {

QBFSolver::QBFSolver() : PseudoBooleanSolver() {
    
}

QBFSolver::~QBFSolver() {
}

void QBFSolver::parse(gzFile in_) {
    Glucose::StreamBuffer in(in_);

    vec<Lit> lits;
    vec<Lit> forall;
    vec<Lit> exists;
    int vars = 0;
    int count = 0;
    inClauses = 0;
    for(;;) {
        if(inClauses > 0 && count == inClauses) break;
        skipWhitespace(in);
        if(*in == EOF) break;
        if(*in == 'p') {
            ++in;
            if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            ++in;
            
            if(eagerMatch(in, "cnf")) {
                vars = parseInt(in);
                inClauses = parseInt(in);
                
                if(vars <= 0) cerr << "PARSE ERROR! The number of variables must be positive." << endl, exit(3);
                if(inClauses <= 0) cerr << "PARSE ERROR! The number of clauses must be positive." << endl, exit(3);
                
                nInVars(vars);
                while(nVars() < nInVars()) { newVar(); g.add(nVars()-1); }
                outSolver.nInVars(nInVars());
                while(outSolver.nVars() < nVars()) outSolver.newVar();
                inSolver.nInVars(nInVars());
                while(inSolver.nVars() < nVars()) inSolver.newVar();
            }
            else {
                cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            }
        }
        else if(*in == 'c')
            skipLine(in);
        else if(eagerMatch(in, "a")) {
            if(forall.size() > 0) cerr << "PARSE ERROR! Input formula is not 2-QBF." << endl, exit(3);
            skipWhitespace(in);
            readClause(in, *this, forall);
            universal.growTo(nVars(), false);
            for(int i = 0; i < forall.size(); i++)
                universal[var(forall[i])] = true;
        }
        else if(eagerMatch(in, "e")) {
            if(exists.size() > 0) cerr << "PARSE ERROR! Input formula is not 2-QBF." << endl, exit(3);
            skipWhitespace(in);
            readClause(in, *this, exists);
        }
        else {
            count++;
            readClause(in, *this, lits);
            
            bool hasUniversal = false;
            for(int j = 0; j < lits.size(); j++) {
                Lit lj = lits[j];
                if(universal[var(lj)]) { hasUniversal = true; continue; }
                for(int k = j+1; k < lits.size(); k++) {
                    Lit lk = lits[k];
                    if(universal[var(lk)]) continue;
                    g.add(var(lj), var(lk));
                }
            }
    
            if(hasUniversal)
                addClause_(lits);
            else
                inSolver.addClause(lits);
        }
    }
    if(count != inClauses)
        cerr << "WARNING! DIMACS header mismatch: wrong number of clauses." << endl, exit(3);
}

void QBFSolver::interrupt() { 
    cout << "s UNKNOWN" << endl;
    //if(lowerbound != -1 && option_fairsat_printmodel) printModel();
    this->exit(1);
}

lbool QBFSolver::solve() {
    outSolver.newVar();
    outSolver.setFrozen(outSolver.nVars()-1, true);
    inSolver.newVar();
    inSolver.setFrozen(inSolver.nVars()-1, true);
    
    for(int i = 0; i < clauses.size(); i++) {
        Clause& clause = ca[clauses[i]];
        outSolver.newVar();
        outSolver.setFrozen(outSolver.nVars()-1, true);
        inSolver.newVar();
        inSolver.setFrozen(inSolver.nVars()-1, true);
        
        Lit c = mkLit(getClauseVar(i));
        
        vec<Lit> inClause;
        inClause.push(~c);
        
        for(int j = 0; j < clause.size(); j++) {
            Lit l = clause[j];
            if(universal[var(l)]) {
                outSolver.addClause(~c, ~l);
                //cout << ~c << " | " << ~l << endl;
            }
            else
                inClause.push(l);
        }
        
        inSolver.addClause(inClause);
    }
    
    
    vec<int> in(clauses.size(), 0);
    int inT = 0;
    vec<Lit> assumptions;
    vec<int> incompatible[2];
    
    
    while(true) {
        assumptions.clear();
        inT++;
        
        for(int i = 0; i < incompatible[0].size(); i++) {
            if(in[incompatible[0][i]] == inT) continue;
            assumptions.push(~mkLit(getClauseVar(incompatible[0][i])));
            in[incompatible[0][i]] = inT;
        }
        for(int i = 0; i <  clauses.size(); i++) {
            if(in[i] == inT) continue;
            assumptions.push(mkLit(getClauseVar(i)));
        }
        inSolver.setAssumptions(assumptions);
        if(inSolver.solve() == l_False) {
            vec<Lit> conflict;
            inSolver.conflict.moveTo(conflict);
            
//            cout << "IN CONFLICT: " << conflict << endl;
            
            assumptions.clear();
            for(int i = 0; i < conflict.size(); i++) {
                assumptions.push(~conflict[i]);
                vec<Lit> tmp;
                assumptions.copyTo(tmp);
                inSolver.setAssumptions(tmp);
                if(inSolver.solve() == l_False) {
                    inSolver.conflict.moveTo(conflict);
                    break;
                }
            }

//            cout << "IN CONFLICT (after): " << conflict << endl;
            
            assumptions.clear();
            for(int i = 0; i < conflict.size(); i++) assumptions.push(~conflict[i]);
            outSolver.setAssumptions(assumptions);
            if(outSolver.solve() == l_True) {
                cout << "INVALID" << endl;
                outSolver.copyModel();
                outSolver.printModel();
                return l_False;
            }
            
//            cout << "OUT CONFLICT: " << outSolver.conflict << endl;
            assert(outSolver.conflict.size() == 2);
            incompatible[0].push(var(outSolver.conflict[0]) - nVars() - 1);
            incompatible[1].push(var(outSolver.conflict[1]) - nVars() - 1);
        }
        else {
            while(true) {
                if(incompatible[0].size() == 0) {
                    cout << "VALID" << endl;
                    return l_True;
                }
                incompatible[0].last() = incompatible[1].last();
                incompatible[1].last() = -1;
                if(incompatible[0].last() != -1) break;
                incompatible[0].pop();
                incompatible[1].pop();
            }
            cout << "INCOMPATIBILITIES:\n\t" << incompatible[0] << "\n\t" << incompatible[1] << endl;
        }
    }

    //if(option_qbf_printmodel) printModel();
    
    return status;
}

} // namespace aspino

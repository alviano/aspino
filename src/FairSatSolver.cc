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

#include "FairSatSolver.h"

#include <core/Dimacs.h>

Glucose::EnumOption option_fairsat_alg("FAIRSAT", "fairsat-alg", "Set search algorithm.", "bb|linear-unsat-sat|binary|progression");
Glucose::BoolOption option_fairsat_printmodel("FAIRSAT", "fairsat-print-model", "Print optimal model if found.", true);

namespace aspino {

template<class B, class Solver>
static void readGE(B& in, Solver& S, WeightConstraint& wc) {
    int parsed_lit, var;
    wc.clear();
    int size = parseInt(in);
    for(int i = 0; i < size; i++) {
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
        var = abs(parsed_lit)-1;
        while (var >= S.nVars()) S.newVar();
        wc.lits.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
    }
    for(int i = 0; i < size; i++) {
        wc.coeffs.push(parseLong(in));
        if(wc.coeffs.last() <= 0) cerr << "PARSE ERROR! The weights must be nonnegative." << endl, exit(3);
    }
    wc.bound = parseLong(in);
    parsed_lit = parseInt(in);
    if (parsed_lit != 0) cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);    
}

template<class B, class Solver>
static void readObjFunc(B& in, Solver& S, vec<Lit>& lits, vec<int64_t>& coeffs) {
    int parsed_lit, var;
    lits.clear();
    coeffs.clear();
    int size = parseInt(in);
    for(int i = 0; i < size; i++) {
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
        var = abs(parsed_lit)-1;
        while (var >= S.nVars()) S.newVar();
        lits.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
    }
    for(int i = 0; i < size; i++) {
        coeffs.push(parseLong(in));
        if(coeffs.last() <= 0) cerr << "PARSE ERROR! The weights must be nonnegative." << endl, exit(3);
    }
    parsed_lit = parseInt(in);
    if (parsed_lit != 0) cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);    
}

void ObjectFunction::init(vec<Lit>& lits_, vec<int64_t>& coeffs_) {
    assert(lits_.size() == coeffs_.size());
    
    lits_.moveTo(lits);
    coeffs_.moveTo(coeffs);
}

FairSatSolver::FairSatSolver() : PseudoBooleanSolver(), lowerbound(-1) {
    if(strcmp(option_fairsat_alg, "bb") == 0) search_alg = &FairSatSolver::search_alg_bb;
    else if(strcmp(option_fairsat_alg, "linear-unsat-sat") == 0) search_alg = &FairSatSolver::search_alg_linear_unsat_sat;
    else if(strcmp(option_fairsat_alg, "binary") == 0) search_alg = &FairSatSolver::search_alg_binary;
    else if(strcmp(option_fairsat_alg, "progression") == 0) search_alg = &FairSatSolver::search_alg_progression;
    else assert(0);
}

FairSatSolver::~FairSatSolver() {
    for(int i = 0; i < objectFunctions.size(); i++) delete objectFunctions[i];
    objectFunctions.clear();
}

void FairSatSolver::parse(gzFile in_) {
    Glucose::StreamBuffer in(in_);

    int64_t objFuncs = -1;
    
    vec<Lit> lits;
    vec<int64_t> coeffs;
    WeightConstraint wc;
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
            
            if(eagerMatch(in, "fcnf")) {
                vars = parseInt(in);
                inClauses = parseInt(in);
                objFuncs = parseLong(in);
                
                if(vars <= 0) cerr << "PARSE ERROR! The number of variables must be positive." << endl, exit(3);
                if(inClauses <= 0) cerr << "PARSE ERROR! The number of clauses must be positive." << endl, exit(3);
                if(objFuncs <= 0) cerr << "PARSE ERROR! The number of object functions must be positive." << endl, exit(3);
                
                nInVars(vars);
                while(nVars() < nInVars()) newVar();
            }
            else {
                cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            }
        }
        else if(*in == 'c')
            skipLine(in);
        else if(*in == 'o') {
            ++in;
            if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            ++in;
            
            lowerbound = parseLong(in) - 1;
        }
        else if(eagerMatch(in, "ge")) {
            count++;
            skipWhitespace(in);
            readGE(in, *this, wc);
            addConstraint(wc);
            morePropagate();
        }
        else {
            count++;
            readClause(in, *this, lits);
            addClause_(lits);
        }
    }
    if(count != inClauses)
        cerr << "WARNING! DIMACS header mismatch: wrong number of clauses." << endl, exit(3);

    count = 0;
    for(;;) {
        skipWhitespace(in);
        if(*in == EOF) break;
        if(*in == 'c')
            skipLine(in);
        else {
            count++;
            readObjFunc(in, *this, lits, coeffs);
            addObjectFunction(lits, coeffs);
        }
    }

    if(count != objFuncs)
        cerr << "WARNING! DIMACS header mismatch: wrong number of object functions." << endl, exit(3);

    upperbound = processObjectFunctions();
    assert(upperbound >= 0);
    freeze();
}

void FairSatSolver::interrupt() { 
    cout << "s UNKNOWN" << endl;
    if(lowerbound != -1 && option_fairsat_printmodel) printModel();
    this->exit(1);
}

void FairSatSolver::addObjectFunction(vec<Lit>& lits, vec<int64_t>& coeffs) {
    ObjectFunction* objF = new ObjectFunction;
    objF->init(lits, coeffs);
    objectFunctions.push(objF);
}

int64_t FairSatSolver::processObjectFunctions() {
    int64_t min = INT64_MAX;
    for(int idx = 0; idx < objectFunctions.size(); idx++) {
        ObjectFunction& objF = *objectFunctions[idx];
        
        WeightConstraint wc;
        objF.lits.copyTo(wc.lits);
        objF.coeffs.copyTo(wc.coeffs);
        for(int i = 0; i < objF.coeffs.size(); i++) wc.bound += objF.coeffs[i];
        objF.sumOfInCoeffs = wc.bound;
        if(wc.bound < min) min = wc.bound;
        
        int64_t n = floor(log2(wc.bound));
        for(int64_t i = 0, w = 1; i <= n; i++, w *= 2) {
            newVar(); setFrozen(nVars()-1, true);
            wc.lits.push(~mkLit(nVars()-1));
            wc.coeffs.push(w);
            objF.selectorVars.push(nVars()-1);
        }
        addConstraint(wc);
    }
    return min;
}

void FairSatSolver::setMinObjectFunction(int64_t min) {
    assumptions.clear();
    for(int i = 0; i < objectFunctions.size(); i++) {
        int64_t mask = objectFunctions[i]->sumOfInCoeffs - min;
        assert(mask >= 0);
        for(int j = objectFunctions[i]->selectorVars.size() - 1; j >= 0; j--) {
//        for(int j = 0; j < objectFunctions[i]->selectorVars.size(); j++) {
            assumptions.push(mkLit(objectFunctions[i]->selectorVars[j], (1ll << j) & mask));
        }
    }
}

int64_t FairSatSolver::analyzeConflict() {
//    cout << assumptions << " " << conflict << endl;
    int64_t max = INT64_MIN;
    int k = conflict.size() - 1;
    for(int i = 0; i < objectFunctions.size(); i++) {
        int64_t sum = objectFunctions[i]->sumOfInCoeffs;
        for(int j = objectFunctions[i]->selectorVars.size() - 1; j >= 0; j--) {
//        for(int j = 0; j < objectFunctions[i]->selectorVars.size(); j++) {
            if(k >= 0 && var(conflict[k]) == objectFunctions[i]->selectorVars[j]) k--;
            else sum -= 1ll << j;;
        }
        if(sum > max) max = sum;
//        cout << k << " " << sum << endl;
    }
//    cout << max << endl;
//    static int count = 0; count++;
//    if(count > 100) exit(1);
    return max;
}

void FairSatSolver::updateLowerBound() {
    copyModel();
    int64_t min = INT64_MAX;
    for(int i = 0; i < objectFunctions.size(); i++) {
        ObjectFunction& objF = *objectFunctions[i];
        objF.modelValue = 0;
        for(int j = 0; j < objF.lits.size(); j++) {
            if(value(objF.lits[j]) != l_False) objF.modelValue += objF.coeffs[j];
        }
        if(objF.modelValue < min) min = objF.modelValue;
    }
    assert_msg(min > lowerbound, "min = " << min << "; lowerbound = " << lowerbound);
    lowerbound = min;
}

lbool FairSatSolver::solve() {
    cout << "c searching in [" << lowerbound << ".." << upperbound << "]" << endl;
    (this->*search_alg)();
    if(lowerbound == -1)
        cout << "s UNSATISFIABLE" << endl;
    else {
        cout << "s OPTIMUM FOUND" << endl;
        if(option_fairsat_printmodel) printModel();
    }
    
    return status;
}

void FairSatSolver::search_alg_bb() {
    for(;;) {
        setMinObjectFunction(lowerbound + 1);
        PseudoBooleanSolver::solve();
        if(status == l_False) break;
        if(status == l_True) {
            updateLowerBound();

            cout << "o " << lowerbound << endl;
//            cout << "c object function values:";
//            for(int i = 0; i < objectFunctions.size(); i++) cout << " " << objectFunctions[i]->modelValue;
//            cout << endl;
//            cout << "c unsatisfied values:";
//            for(int i = 0; i < objectFunctions.size(); i++) cout << " " << (objectFunctions[i]->sumOfInCoeffs - objectFunctions[i]->modelValue);
//            cout << endl;
        }
        if(lowerbound == upperbound) break;
    }
}

void FairSatSolver::search_alg_linear_unsat_sat() {
    for(;;) {
        setMinObjectFunction(upperbound);
        PseudoBooleanSolver::solve();
        if(status == l_False) {
            cout << (upperbound - analyzeConflict() + 1) << endl;
            upperbound = analyzeConflict() - 1;
            cout << "c ub " << upperbound << endl;
        }
        else if(status == l_True) {
            updateLowerBound();

            cout << "o " << lowerbound << endl;
            break;
        }
        if(lowerbound == upperbound) break;
    }
}

void FairSatSolver::search_alg_binary() {
    for(;;) {
        assert(lowerbound < upperbound);
        int64_t mid = (lowerbound + upperbound) / 2;
        if(mid == lowerbound) mid++;
        setMinObjectFunction(mid);
        PseudoBooleanSolver::solve();
        if(status == l_False) {
            upperbound = mid-1;
            cout << "c ub " << upperbound << endl;
        }
        else if(status == l_True) {
            updateLowerBound();

            cout << "o " << lowerbound << endl;
//            cout << "c object function values:";
//            for(int i = 0; i < objectFunctions.size(); i++) cout << " " << objectFunctions[i]->modelValue;
//            cout << endl;
//            cout << "c unsatisfied values:";
//            for(int i = 0; i < objectFunctions.size(); i++) cout << " " << (objectFunctions[i]->sumOfInCoeffs - objectFunctions[i]->modelValue);
//            cout << endl;
        }
        if(lowerbound == upperbound) break;
    }    
}

void FairSatSolver::search_alg_progression() {
    int64_t progression = 1;
    for(;;) {
        assert(lowerbound < upperbound);
        if(lowerbound + progression > upperbound) progression = 1;
        setMinObjectFunction(lowerbound + progression);
//        cout << progression << endl;
//        setConfBudget(100);
        PseudoBooleanSolver::solve();
//        budgetOff();
        if(status == l_False) {
            upperbound = lowerbound + progression - 1;
            cout << "c ub " << upperbound << endl;
        }
        else if(status == l_True) {
            updateLowerBound();

            cout << "o " << lowerbound << endl;
//            cout << "c object function values:";
//            for(int i = 0; i < objectFunctions.size(); i++) cout << " " << objectFunctions[i]->modelValue;
//            cout << endl;
//            cout << "c unsatisfied values:";
//            for(int i = 0; i < objectFunctions.size(); i++) cout << " " << (objectFunctions[i]->sumOfInCoeffs - objectFunctions[i]->modelValue);
//            cout << endl;
        }
//        else { progression = 1; continue; }
        if(lowerbound == upperbound) break;
        progression *= 2;
    }    
}

} // namespace aspino

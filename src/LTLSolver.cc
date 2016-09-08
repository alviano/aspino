/*
 *  Copyright (C) 2016 Mario Alviano (mario@alviano.net)
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

#include "LTLSolver.h"

#include "utils/ParseUtils.h"
#include "core/Dimacs.h"

namespace aspino {

void LTLFormula::init(LTLSolver& solver) {
    vars.clear();
    solver.newVar();
    vars.push(solver.nVars()-1);
}

void LTLFormula::maxLength(LTLSolver& solver, int value) {
    while(value > vars.size()) { 
        solver.newVar();
        vars.push(solver.nVars()-1);
    }
}

void LTLConstant::addClauses(LTLSolver& solver, int step) {
    if(step == 0) solver.addClause(mkLit(vars[0], !value));
}

void LTLFormulaX::addClauses(LTLSolver& solver, int i) {
    solver.addClause(~mkLit(var(i)), mkLit(solver.getFormula(child).var(i+1)));
    solver.addClause(~mkLit(var(i)), ~mkLit(solver.getLengthVar(i)));
    solver.addClause(mkLit(var(i)), ~mkLit(solver.getFormula(child).var(i+1)), mkLit(solver.getLengthVar(i)));
}

void LTLFormulaW::addClauses(LTLSolver& solver, int i) {
    solver.addClause(~mkLit(var(i)), mkLit(solver.getFormula(child).var(i+1)), mkLit(solver.getLengthVar(i)));
    solver.addClause(mkLit(var(i)), ~mkLit(solver.getFormula(child).var(i+1)));
    solver.addClause(mkLit(var(i)), ~mkLit(solver.getLengthVar(i)));
}

void LTLFormulaF::addClauses(LTLSolver& solver, int i) {
    solver.addClause(~mkLit(var(i)), mkLit(solver.getFormula(child).var(i)), mkLit(var(i+1)));
    solver.addClause(~mkLit(var(i)), mkLit(solver.getFormula(child).var(i)), ~mkLit(solver.getLengthVar(i)));
    solver.addClause(mkLit(var(i)), ~mkLit(solver.getFormula(child).var(i)));
    solver.addClause(mkLit(var(i)), ~mkLit(var(i+1)), mkLit(solver.getLengthVar(i)));
}

void LTLFormulaG::addClauses(LTLSolver& solver, int i) {
    solver.addClause(~mkLit(var(i)), mkLit(solver.getFormula(child).var(i)));
    solver.addClause(~mkLit(var(i)), mkLit(var(i+1)), mkLit(solver.getLengthVar(i)));
    solver.addClause(mkLit(var(i)), ~mkLit(solver.getFormula(child).var(i)), ~mkLit(var(i+1)));
    solver.addClause(mkLit(var(i)), ~mkLit(solver.getFormula(child).var(i)), ~mkLit(solver.getLengthVar(i)));
}

void LTLFormulaU::addClauses(LTLSolver& solver, int i) {
    solver.addClause(~mkLit(var(i)), mkLit(solver.getFormula(right).var(i)), mkLit(solver.getFormula(left).var(i)));
    solver.addClause(~mkLit(var(i)), mkLit(solver.getFormula(right).var(i)), mkLit(var(i+1)));
    solver.addClause(~mkLit(var(i)), mkLit(solver.getFormula(right).var(i)), ~mkLit(solver.getLengthVar(i)));

    solver.addClause(mkLit(var(i)), ~mkLit(solver.getFormula(right).var(i)));
    vec<Lit> lits;
    lits.push(mkLit(var(i)));
    lits.push(~mkLit(solver.getFormula(left).var(i)));
    lits.push(~mkLit(var(i+1)));
    lits.push(mkLit(solver.getLengthVar(i)));
    solver.addClause(lits);
}

void LTLFormulaNot::addClauses(LTLSolver& solver, int i) {
    solver.addClause(~mkLit(var(i)), ~mkLit(solver.getFormula(child).var(i)));
    solver.addClause(mkLit(var(i)), mkLit(solver.getFormula(child).var(i)));
}

void LTLFormulaAnd::addClauses(LTLSolver& solver, int i) {
    solver.addClause(~mkLit(var(i)), mkLit(solver.getFormula(left).var(i)));
    solver.addClause(~mkLit(var(i)), mkLit(solver.getFormula(right).var(i)));
    solver.addClause(mkLit(var(i)), ~mkLit(solver.getFormula(left).var(i)), ~mkLit(solver.getFormula(right).var(i)));
}

void LTLFormulaOr::addClauses(LTLSolver& solver, int i) {
    solver.addClause(~mkLit(var(i)), mkLit(solver.getFormula(left).var(i)), mkLit(solver.getFormula(right).var(i)));
    solver.addClause(mkLit(var(i)), ~mkLit(solver.getFormula(left).var(i)));
    solver.addClause(mkLit(var(i)), ~mkLit(solver.getFormula(right).var(i)));
}

void LTLFormulaImpl::addClauses(LTLSolver& solver, int i) {
    solver.addClause(~mkLit(var(i)), ~mkLit(solver.getFormula(left).var(i)), mkLit(solver.getFormula(right).var(i)));
    solver.addClause(mkLit(var(i)), mkLit(solver.getFormula(left).var(i)));
    solver.addClause(mkLit(var(i)), ~mkLit(solver.getFormula(right).var(i)));
}

void LTLFormulaEq::addClauses(LTLSolver& solver, int i) {
    solver.addClause(~mkLit(var(i)), ~mkLit(solver.getFormula(left).var(i)), mkLit(solver.getFormula(right).var(i)));
    solver.addClause(~mkLit(var(i)), mkLit(solver.getFormula(left).var(i)), ~mkLit(solver.getFormula(right).var(i)));
    solver.addClause(mkLit(var(i)), ~mkLit(solver.getFormula(left).var(i)), ~mkLit(solver.getFormula(right).var(i)));
    solver.addClause(mkLit(var(i)), mkLit(solver.getFormula(left).var(i)), mkLit(solver.getFormula(right).var(i)));
}


void LTLVariable::printModel(const vec<lbool>& model, int step) const {
    if(model[vars[step]] == l_True) cout << " " << name;
    else if(model[vars[step]] == l_False) cout << " -" << name;
}


void LTLSolver::parse(gzFile in_) {
    Glucose::StreamBuffer in(in_);
    int id;
    char buff[1024];

    for(;;) {
        skipWhitespace(in);
        if(*in == EOF) break;
        if(*in == 'p') {
            ++in;
            if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            ++in;
            
            if(eagerMatch(in, "ltlcnf")) {
                maxsize = parseInt(in);
            }
            else {
                cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            }
        }
        else if(*in == 'c')
            skipLine(in);
        else if(*in == 'f') {
            ++in;
            if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            ++in;
            
            id = parseInt(in);
            roots.push(id);
            while(id >= formulas.size()) formulas.push(NULL);
        }
        else {
            id = parseInt(in);
            ++in;
            while(id >= formulas.size()) formulas.push(NULL);
            assert(formulas[id] == NULL);
            
            if(*in == 'X') {
                ++in;
                if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
                formulas[id] = new LTLFormulaX(id, parseInt(in));
            }
            else if(*in == 'W') {
                ++in;
                if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
                formulas[id] = new LTLFormulaW(id, parseInt(in));
            }
            else if(*in == 'F') {
                if(eagerMatch(in, "FALSE")) formulas[id] = new LTLConstant(id, false);
                else {
                    if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
                    formulas[id] = new LTLFormulaF(id, parseInt(in));
                }
            }
            else if(eagerMatch(in, "TRUE")) {
                formulas[id] = new LTLConstant(id, true);
            }
            else if(*in == 'G') {
                ++in;
                if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
                formulas[id] = new LTLFormulaG(id, parseInt(in));
            }
            else if(*in == 'U') {
                ++in;
                if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
                formulas[id] = new LTLFormulaU(id, parseInt(in), parseInt(in));
            }
            else if(*in == '!') {
                ++in;
                if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
                formulas[id] = new LTLFormulaNot(id, parseInt(in));
            }
            else if(*in == '&') {
                ++in;
                if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
                formulas[id] = new LTLFormulaAnd(id, parseInt(in), parseInt(in));
            }
            else if(*in == '|') {
                ++in;
                if(*in != ' ') cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
                formulas[id] = new LTLFormulaOr(id, parseInt(in), parseInt(in));
            }
            else if(*in == '-') {
                if(!eagerMatch(in, "->")) cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
                formulas[id] = new LTLFormulaImpl(id, parseInt(in), parseInt(in));
                
            }
            else if(*in == '<') {
                if(!eagerMatch(in, "<->")) cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
                formulas[id] = new LTLFormulaEq(id, parseInt(in), parseInt(in));
            }
            else {
                int i = 0;
                while(*in != '\n') { buff[i++] = *in; ++in; assert(i < 1024); }
                buff[i] = '\0';
                formulas[id] = new LTLVariable(id, buff);
            }
        }
    }
}

LTLSolver::~LTLSolver() {
    for(int i = 0; i < formulas.size(); i++) delete formulas[i];
}


lbool LTLSolver::solve() {
    lengthVars.clear();
    newVar();
    lengthVars.push(nVars()-1);

    for(int i = 1; i < formulas.size(); i++) {
        assert(formulas[i] != NULL);
        formulas[i]->init(*this);
    }

    // roots must be true at step 0
    for(int i = 0; i < roots.size(); i++) {
        LTLFormula* formula = formulas[roots[i]];
        assert(formula != NULL);
        addClause(mkLit(formula->var(0)));
    }
    
    for(int iprogression = 0, prec = 0, progression = 1; iprogression < maxsize; iprogression++, prec = progression, progression *= 2) {
        // l[i] -> l[i+1]
        while(progression > lengthVars.size()) {
            newVar();
            addClause(~mkLit(lengthVars.last()), mkLit(nVars()-1));
            lengthVars.push(nVars()-1);
        }

        for(int i = 1; i < formulas.size(); i++) formulas[i]->maxLength(*this, progression+1);
        for(int i = 1; i < formulas.size(); i++) {
            assert(formulas[i] != NULL);
            for(int j = prec; j < progression; j++) {
                formulas[i]->addClauses(*this, j);
            }
        }
        

//for(int i = 0; i < clauses.size();i++) cout << ca[clauses[i]] << endl;
        
        // assume the last l[i]
        assumptions.clear();
        assumptions.push(mkLit(lengthVars.last()));
        
        cout << "c solving..." << endl;
        SatSolver::solve();
        cout << "c done" << endl;
        if(status == l_True) {
            cout << "s SATISFIABLE" << endl;
            copyModel();
            printModel();
            return l_True;
        }
        
        if(conflict.size() == 0) { cout << "o infinity" << endl; break; }
        assert(conflict.size() == 1);
        cancelUntil(0);
        addClause(conflict);
        cout << "o " << progression << endl;
        
        for(int i = 1; i < formulas.size(); i++) {
            //cout << value(formulas[i]->var(0)) << endl;
        }
    }
    
    cout << "s UNSATISFIABLE" << endl;    
    return l_False;
}

void LTLSolver::printModel() const {
    for(int i = 0; i < lengthVars.size(); i++) {
        cout << "c Step " << i+1 << endl;
        cout << "v";
        for(int j = 1; j < formulas.size(); j++) formulas[j]->printModel(model, i);
        cout << endl;
        if(model[lengthVars[i]] == l_True) break;
    }
}

} // namespace aspino

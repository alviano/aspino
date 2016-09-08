/*
 *  Copyright (C) 2016  Mario Alviano (mario@alviano.net)
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

#ifndef __LTLSolver_h__
#define __LTLSolver_h__

#include "SatSolver.h"

namespace aspino {

class LTLSolver;
    
class LTLFormula {
public:
    LTLFormula(int id_) : id(id_) {}
    virtual ~LTLFormula() {}
    
    void init(LTLSolver& solver);
    virtual void maxLength(LTLSolver& solver, int value);
    virtual Var var(int step) const { assert(step < vars.size()); return vars[step]; }
    
    virtual void addClauses(LTLSolver& solver, int i) = 0;
    
    virtual void printModel(const vec<lbool>&, int) const {}
    
protected:
    int id;
    vec<Var> vars;
};

class LTLConstant : public LTLFormula {
public:
    LTLConstant(int id, bool value_) : LTLFormula(id), value(value_) {}
    
    virtual void addClauses(LTLSolver& solver, int step);
    virtual void maxLength(LTLSolver&, int) {}
    virtual Var var(int) const { assert(vars.size() == 1); return vars[0]; }
    
protected:    
    bool value;
};

class LTLVariable : public LTLFormula {
public:
    LTLVariable(int id, const char* name_) : LTLFormula(id), name(name_) {}
    
    virtual void addClauses(LTLSolver&, int) {};
    
    virtual void printModel(const vec<lbool>& model, int step) const;

protected:    
    string name;
};

class LTLFormulaUnary : public LTLFormula {
public:
    LTLFormulaUnary(int id_, int child_) : LTLFormula(id_), child(child_) {}
    
    int child;
};

class LTLFormulaBinary : public LTLFormula {
public:
    LTLFormulaBinary(int id_, int left_, int right_) : LTLFormula(id_), left(left_), right(right_) {}
    
    int left;
    int right;
};

class LTLFormulaX : public LTLFormulaUnary {
public:
    LTLFormulaX(int id, int child) : LTLFormulaUnary(id, child) {}
    
    virtual void addClauses(LTLSolver& solver, int i);
};

class LTLFormulaW : public LTLFormulaUnary {
public:
    LTLFormulaW(int id, int child) : LTLFormulaUnary(id, child) {}
    
    virtual void addClauses(LTLSolver& solver, int i);
};

class LTLFormulaF : public LTLFormulaUnary {
public:
    LTLFormulaF(int id, int child) : LTLFormulaUnary(id, child) {}
    
    virtual void addClauses(LTLSolver& solver, int i);
};

class LTLFormulaG : public LTLFormulaUnary {
public:
    LTLFormulaG(int id, int child) : LTLFormulaUnary(id, child) {}
    
    virtual void addClauses(LTLSolver& solver, int i);
};

class LTLFormulaU : public LTLFormulaBinary {
public:
    LTLFormulaU(int id, int left, int right) : LTLFormulaBinary(id, left, right) {}
    
    virtual void addClauses(LTLSolver& solver, int i);
};

class LTLFormulaNot : public LTLFormulaUnary {
public:
    LTLFormulaNot(int id, int child) : LTLFormulaUnary(id, child) {}
    
    virtual void addClauses(LTLSolver& solver, int i);
};

class LTLFormulaAnd : public LTLFormulaBinary {
public:
    LTLFormulaAnd(int id, int left, int right) : LTLFormulaBinary(id, left, right) {}
    
    virtual void addClauses(LTLSolver& solver, int i);
};

class LTLFormulaOr : public LTLFormulaBinary {
public:
    LTLFormulaOr(int id, int left, int right) : LTLFormulaBinary(id, left, right) {}
    
    virtual void addClauses(LTLSolver& solver, int i);
};

class LTLFormulaImpl : public LTLFormulaBinary {
public:
    LTLFormulaImpl(int id, int left, int right) : LTLFormulaBinary(id, left, right) {}
    
    virtual void addClauses(LTLSolver& solver, int i);
};

class LTLFormulaEq : public LTLFormulaBinary {
public:
    LTLFormulaEq(int id, int left, int right) : LTLFormulaBinary(id, left, right) {}
    
    virtual void addClauses(LTLSolver& solver, int i);
};


    
class LTLSolver : public SatSolver {
public:
    virtual ~LTLSolver();
    virtual void interrupt() { this->exit(-1); }
    
    virtual void parse(gzFile in);

    virtual void newVar() { SatSolver::newVar(); setFrozen(nVars()-1, false); }

    virtual lbool solve();
    virtual lbool solve(int) { return solve(); }
    
    const LTLFormula& getFormula(int idx) const { return *formulas[idx]; }
    Var getLengthVar(int idx) const { assert(idx < lengthVars.size()); return lengthVars[idx]; }
    
    virtual void printModel() const;

private:
    int maxsize;
    vec<int> roots;
    vec<LTLFormula*> formulas;
    vec<Var> lengthVars;
};

} // namespace aspino

#endif
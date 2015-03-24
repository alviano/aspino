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

#ifndef __AspSolver_h__
#define __AspSolver_h__

#include "MaxSatSolver.h"
#include "glucose-syrup/mtl/Map.h"

using Glucose::Map;

namespace aspino {

class Literal {
public:
    Literal(int id=0, int type=0);
    inline Literal(const Literal& init) : id(init.id), type(init.type) {}
    
    inline Lit toLit() const { return mkLit(id, type == POS || type == DNEG); }
    
    inline bool operator==(const Literal& right) const { return id == right.id && type == right.type; }
    inline bool operator!=(const Literal& right) const { return !(*this == right); }
    
    inline bool head() const { return SHEAD <= type && type <= UHEAD; }
    inline bool body() const { return POS <= type && type <= DNEG; }
    
    unsigned id : 29;
    unsigned type : 3;
    
    static const int SHEAD;
    static const int UHEAD;
    static const int POS;
    static const int NEG;
    static const int DNEG;

private:
};

class AspSolver : public MaxSatSolver {
public:
    AspSolver();
    virtual ~AspSolver();
    
    virtual void interrupt() { cout << "s UNKNOWN" << endl; this->exit(1); }
    
    virtual void parse(gzFile in);
    
    virtual void newVar();
    
    virtual bool eliminate(bool turn_off_elim);
    virtual lbool solve(int n);
    
    virtual void printModel() const;
    
private:
    vec<vec<Literal>*> program;
    
    vec<int> idmap;
    Map<int, string> namemap;
    vec<vec<vec<Literal>*> > occ[5];
    vec<unsigned> tag;
    unsigned tagCalls;
    
    vec<bool> supported;
    vec<Var> supportInference;
    
    unsigned propagated;
    
    int nModels;
    
    int getId(int input_id);
    string getName(int atom) const;
    
    void parseNormalRule(Glucose::StreamBuffer& in);
    void parseChoiceRule(Glucose::StreamBuffer& in);

    void propagate();
    void propagateTrue(Var v);
    void propagateFalse(Var v);
    void onTrueHead(vec<Literal>& rule, Literal headAtom);
    void onFalseHead(vec<Literal>& rule, Literal headAtom);
    void onTrueBody(vec<Literal>& rule, Literal bodyAtom);
    void onFalseBody(vec<Literal>& rule);
    
    void finalPropagation();
    void processComponents();
    void completion();
    void clearParsingStructures();
};
    
} // namespace aspino

#endif
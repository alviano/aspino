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
    friend ostream& operator<<(ostream& out, const Literal& lit);
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

class Sum {
    friend ostream& operator<<(ostream& out, const Sum& sum);
public:
    inline Sum() : bound(0), posTrueWeight(0) {}
    
    Literal head;
    vec<Literal> body;
    vec<long> weights;
    long bound;
    long posTrueWeight;
    
    inline int size() const { return body.size(); }
    inline void shrink(int size) { body.shrink(size); weights.shrink(size); }
    
    static inline int type2idx(int type) { assert(type == Literal::POS || type == Literal::NEG); return type == Literal::POS ? 0 : 1; }
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
    
protected:
    virtual CRef morePropagate();
    virtual bool moreReason(Lit lit, vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC);
    void _moreReasonWF(Lit lit, vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC);
    virtual bool moreReason(Lit lit);
    void _moreReasonWF(Lit lit);
    virtual bool moreConflict(vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC);
    void _moreConflictWF(vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC);
    
private:
    vec<vec<Literal>*> program;
    vec<Sum*> sums;
    
    vec<int> idmap;
    Map<int, string> namemap;
    vec<vec<vec<Literal>*> > occ[5];
    vec<Sum*> sumId;
    vec<vec<Sum*> > sumOcc[2];
    vec<unsigned> tag;
    unsigned tagCalls;
    
    vec<bool> supported;
    vec<Var> supportInference;
    
    int nextToPropagate;
    
    int nModels;
    
    vec<vec<Literal>*> recursiveRules;
    vec<Lit> body;
    vec<vec<int> > recBody;
    vec<vec<int> > isBodyOf[2];
    vec<int> withoutSourcePointer;
    vec<int> possibleSourcePointerOf;
    vec<int> sourcePointer;
    vec<vec<int> > possibleSourcePointers;
    vec<vec<int> > inBody[2];
    vec<vec<int> > inRecBody;
    vec<int> unfounded;
    
    vec<int> moreReasonWF;
    vec<Var> moreReasonVars;
    vec<vec<Lit> > moreReasonWFVec;
//    WeightConstraint* moreConflictWC;
//    Lit moreConflictLit;
    int moreConflictWF;

    
    CRef morePropagate(Lit lit);
    bool findSourcePointers();
    
    int getId(int input_id);
    string getName(int atom) const;
    
    void parseNormalRule(Glucose::StreamBuffer& in);
    void parseChoiceRule(Glucose::StreamBuffer& in);
    void parseCountRule(Glucose::StreamBuffer& in);
    void parseSumRule(Glucose::StreamBuffer& in);
    
    void add(Sum* sum);

    void propagate();
    void propagateTrue(Var v);
    void propagateFalse(Var v);
    void onTrueHead(vec<Literal>& rule, Literal headAtom);
    void onFalseHead(vec<Literal>& rule, Literal headAtom);
    void onTrueBody(vec<Literal>& rule, Literal bodyAtom);
    void onFalseBody(vec<Literal>& rule);
    void onTrueHead(Sum& rule);
    void onFalseHead(Sum& rule);
    void onTrueBody(Sum& rule, Literal lit);
    void onFalseBody(Sum& rule, Literal lit);
    
    void finalPropagation();
    void processComponents();
    void completion(Var i, vec<Lit>& supp);
    void clearParsingStructures();
    
    virtual void onCancel();
    
    bool hasSupportInference(Var atom) const;
};
    
} // namespace aspino

#endif
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

#include "AspSolver.h"

#include <sstream>
#include "DepGraph.h"

extern Glucose::IntOption option_n;

namespace aspino {

template<class B>
static int parseLine(B& in, char* out, int maxLenght = 1024) {
    skipWhitespace(in);
    int i = 0;
    for(;;) {
        if(i >= maxLenght) return -1;
        if(*in == '\n' || *in == EOF) break;
        out[i++] = *in;
        ++in;
    }
    out[i] = '\0';
    return i;
}

const int Literal::SHEAD = 0;
const int Literal::UHEAD = 1;
const int Literal::POS = 2;
const int Literal::NEG = 3;
const int Literal::DNEG = 4;

Literal::Literal(int id_, int type_)
: id(id_), type(type_) {
    assert(SHEAD <= type && type <= DNEG);
}


AspSolver::AspSolver() {
    tagCalls = 0;
    nextToPropagate = 0;
    nModels = 0;
}

AspSolver::~AspSolver() {
    
}

void AspSolver::parseNormalRule(Glucose::StreamBuffer& in) {
    trace(asp_pre, 5, "Parsing normal rule...");
    int head = parseInt(in);
    vec<Literal>* rule = new vec<Literal>();
    
    int size = parseInt(in);
    int neg = parseInt(in);
    

    tagCalls += 3;
    if(head != 1) {
        rule->push(Literal(getId(head), Literal::SHEAD));
        if(supported[getId(head)]) {
            trace(asp_pre, 8, "Already supported head. Skip this rule.");
            delete rule;
            while(size-- > 0) parseInt(in);
            return;
        }
        tag[getId(head)] = tagCalls - 2;
    }
    
    while(size-- > 0) {
        int atom = getId(parseInt(in));
        if(neg-- > 0) {
            if(tag[atom] == tagCalls - 2) {
                trace(asp_pre, 8, "Head atom in negative body. Replace by constraint.");
                rule->operator[](0) = rule->operator[](rule->size()-1);
                rule->shrink(1);
                head = 1;
                continue;
            }
            if(tag[atom] == tagCalls - 1) continue;
            if(value(atom) == l_False) continue;
            if(value(atom) == l_True) {
                trace(asp_pre, 8, "Negative body is false. Skip rule.");
                delete rule;
                while(size-- > 0) parseInt(in);
                return;
            }

            rule->push(Literal(atom, Literal::NEG));
            tag[atom] = tagCalls - 1;
        }
        else {
            if(tag[atom] == tagCalls - 2 || tag[atom] == tagCalls - 1 || value(atom) == l_False) {
                trace(asp_pre, 8, "Head atom in positive body, or inconsistent body, or false body. Skip rule.");
                delete rule;
                while(size-- > 0) parseInt(in);
                return;
            }
            if(tag[atom] == tagCalls) continue;
                
            rule->push(Literal(atom, Literal::POS));
            tag[atom] = tagCalls;
        }
    }
    
    if(head != 1) {
        if(rule->size() == 1) {
            trace(asp_pre, 8, "It is a fact!");
            supported[rule->operator[](0).id] = true;
            addClause(rule->operator[](0).toLit());
            delete rule;
            return;
        }
        program.push(rule);
        for(int i = 0; i < rule->size(); i++) {
            Literal lit = rule->operator[](i);
            occ[lit.type][lit.id].push(rule);
        }
    }
    
    vec<Lit> lits;
    for(int i = 0; i < rule->size(); i++) {
        Literal lit = rule->operator[](i);
        lits.push(lit.toLit());
    }
    addClause(lits);
}

void AspSolver::parseChoiceRule(Glucose::StreamBuffer& in) {
    trace(asp_pre, 5, "Parsing choice rule...");
    
    int headSize = parseInt(in);
    vec<int> head;
    while(headSize-- > 0) {
        int atom = getId(parseInt(in));
        if(value(atom) == l_True && supported[atom]) continue;
        head.push(atom);
    }
    
    vec<Literal>* rule = new vec<Literal>();
    
    int size = parseInt(in);
    int neg = parseInt(in);

    tagCalls += 2;
    
    while(size-- > 0) {
        int atom = getId(parseInt(in));
        if(neg-- > 0) {
            if(tag[atom] == tagCalls - 1) continue;
            if(value(atom) == l_False) continue;
            if(value(atom) == l_True) {
                trace(asp_pre, 8, "Negative body is false. Skip rule.");
                delete rule;
                while(size-- > 0) parseInt(in);
                return;
            }

            rule->push(Literal(atom, Literal::NEG));
            tag[atom] = tagCalls - 1;
        }
        else {
            if(tag[atom] == tagCalls - 1 || value(atom) == l_False) {
                trace(asp_pre, 8, "Inconsistent body, or false body. Skip rule.");
                delete rule;
                while(size-- > 0) parseInt(in);
                return;
            }
            if(tag[atom] == tagCalls) continue;
                
            rule->push(Literal(atom, Literal::POS));
            tag[atom] = tagCalls;
        }
    }
    
    int j = 0;
    for(int i = 0; i < head.size(); i++) {
        if(tag[head[i]] == tagCalls - 1) continue;
        if(tag[head[i]] == tagCalls) continue;
        head[j++] = head[i];
    }
    head.shrink(head.size() - j);
    if(head.size() == 0) {
        delete rule;
        return;
    }
    
    for(int i = 0; i < head.size(); i++) {
        vec<Literal>* r;
        if(i == head.size()-1) r = rule;
        else {
            r = new vec<Literal>();
            rule->copyTo(*r);
        }
        
        r->push(Literal(head[i], Literal::SHEAD));
        r->push(Literal(head[i], Literal::DNEG));
        program.push(r);

        for(int j = 0; j < r->size(); j++) {
            Literal lit = r->operator[](j);
            occ[lit.type][lit.id].push(r);
        }
    }
}

void AspSolver::parseCountRule(Glucose::StreamBuffer& in) {
    trace(asp_pre, 5, "Parsing count rule...");
    int head = parseInt(in);
    assert(head != 1);
    vec<Literal>* rule = new vec<Literal>();
    
    int size = parseInt(in);
    int neg = parseInt(in);
    int bound = parseInt(in);    

    while(size-- > 0) {
        int atom = getId(parseInt(in));
        if(neg-- > 0) rule->push(Literal(atom, Literal::NEG));
        else rule->push(Literal(atom, Literal::POS));
    }
    delete rule;
}

void AspSolver::parseSumRule(Glucose::StreamBuffer& in) {
    trace(asp_pre, 5, "Parsing sum rule...");
    int head = parseInt(in);
    assert(head != 1);
    vec<Literal>* rule = new vec<Literal>();
    
    int bound = parseInt(in);    
    int size = parseInt(in);
    int neg = parseInt(in);

    while(size-- > 0) {
        int atom = getId(parseInt(in));
        if(neg-- > 0) rule->push(Literal(atom, Literal::NEG));
        else rule->push(Literal(atom, Literal::POS));
    }
    
    vec<int> weights;
    while(weights.size() < rule->size()) weights.push(parseInt(in));
    delete rule;
}

void AspSolver::parse(gzFile in_) {
    trace(asp_pre, 1, "Start parsing...");
    Glucose::StreamBuffer in(in_);
    
    for(;;) {
        skipWhitespace(in);
        if(*in == EOF) cerr << "PARSE ERROR! Unexpected EOF." << endl, exit(3);
        int type = parseInt(in);
        if(type == 0) break;
        else if(type == 1) parseNormalRule(in);
        else if(type == 3) parseChoiceRule(in);
        else if(type == 2) parseCountRule(in);
        else if(type == 5) parseSumRule(in);
        else cerr << "PARSE ERROR! Unexpected rule type: " << type << endl, exit(3);
        propagate();
    }
    
    trace(asp_pre, 1, "Parsing names...");
    int id;
    char buff[1024];
    for(;;) {
        skipWhitespace(in);
        id = parseInt(in);
        if(id == 0) break;
        if(parseLine(in, buff, 1024) == -1) cerr << "PARSE ERROR! Atom name too long." << endl, exit(3);
        assert(!namemap.has(idmap[id]));
        namemap.insert(idmap[id], buff);
        trace(asp_pre, 15, "Atom " << idmap[id] << ": id=" << id << " name=" << getName(idmap[id]));
    }
    
    finalPropagation();
    processComponents();
    clearParsingStructures();
//    for(int i = 0; i < nVars(); i++) setFrozen(i, true); // TODO: remove
    
    trace(asp_pre, 1, "All done.");
}

int AspSolver::getId(int input_id) {
    assert(idmap.size() == occ[0].size() && idmap.size() == occ[1].size() && idmap.size() == occ[2].size() && idmap.size() == occ[3].size() && idmap.size() == occ[4].size());
    while(input_id >= idmap.size()) {
        idmap.push(-1);
        for(int i = Literal::SHEAD; i <= Literal::DNEG; i++) occ[i].push();
    }
    assert(idmap.size() == occ[0].size() && idmap.size() == occ[1].size() && idmap.size() == occ[2].size() && idmap.size() == occ[3].size() && idmap.size() == occ[4].size());
    if(idmap[input_id] == -1) { 
        trace(asp_pre, 10, "Atom " << nVars() << ": id=" << input_id);
        idmap[input_id] = nVars(); 
        newVar(); 
    }
    return idmap[input_id];
}

string AspSolver::getName(int atom) const {
    if(namemap.has(atom)) return namemap[atom];
    stringstream ss;
    ss << atom;
    return ss.str();
}

void AspSolver::newVar() {
    MaxSatSolver::newVar();
    tag.push(0);
    supported.push(false);
    sourcePointer.push(-1);
    possibleSourcePointers.push();
    inBody[0].push();
    inBody[1].push();
    inRecBody.push();
    isBodyOf[0].push();
    isBodyOf[1].push();
    moreReasonWF.push(-1);
}

bool AspSolver::eliminate(bool turn_off_elim) {
    if(option_n != 1) {
        for(int i = 0; i < nVars(); i++)
            if(namemap.has(i))
                setFrozen(i, true);
    }
    return MaxSatSolver::eliminate(turn_off_elim);
}

lbool AspSolver::solve(int n) {
    for(;;) {
        PseudoBooleanSolver::solve_();
        if(status == l_False) break;

        nModels++;
        copyModel();
        printModel();

        if(nModels == n) break;
        if(decisionLevel() == 0) break;
        learnClauseFromModel();
    }
    
    if(nModels == 0) cout << "UNSATISFIABLE" << endl;
    
    return nModels == 0 ? l_False : l_True;
}

void AspSolver::printModel() const {
    cout << "Answer: " << nModels << endl;
    for(int i = 0; i < nVars(); i++)
        if(namemap.has(i) && model[i] == l_True)
            cout << namemap[i] << " ";
    cout << endl;
}

void AspSolver::propagate() {
    trace(asp_pre, 5, "Start propagate()...");
    while(nextToPropagate < nextToPropagateByUnit()) {
        Var v = var(trail[nextToPropagate++]);
        if(value(v) == l_True) propagateTrue(v);
        else propagateFalse(v);
    }
}

void AspSolver::propagateTrue(Var v) {
    trace(asp_pre, 10, "Propagate true: " << getName(v));
    
    for(int i = 0; i < occ[Literal::SHEAD][v].size(); i++)
        onTrueHead(*(occ[Literal::SHEAD][v][i]), Literal(v, Literal::SHEAD));
    for(int i = 0; i < occ[Literal::UHEAD][v].size(); i++)
        onTrueHead(*(occ[Literal::UHEAD][v][i]), Literal(v, Literal::UHEAD));
        
    for(int i = 0; i < occ[Literal::POS][v].size(); i++)
        onTrueBody(*(occ[Literal::POS][v][i]), Literal(v, Literal::POS));
    for(int i = 0; i < occ[Literal::DNEG][v].size(); i++)
        onTrueBody(*(occ[Literal::DNEG][v][i]), Literal(v, Literal::DNEG));

    for(int i = 0; i < occ[Literal::NEG][v].size(); i++)
        onFalseBody(*(occ[Literal::NEG][v][i]));
}

void AspSolver::propagateFalse(Var v) {
    trace(asp_pre, 10, "Propagate false: " << getName(v));

    for(int i = 0; i < occ[Literal::SHEAD][v].size(); i++)
        onFalseHead(*(occ[Literal::SHEAD][v][i]), Literal(v, Literal::SHEAD));
    for(int i = 0; i < occ[Literal::UHEAD][v].size(); i++)
        onFalseHead(*(occ[Literal::UHEAD][v][i]), Literal(v, Literal::UHEAD));


    for(int i = 0; i < occ[Literal::NEG][v].size(); i++)
        onTrueBody(*(occ[Literal::NEG][v][i]), Literal(v, Literal::NEG));

    for(int i = 0; i < occ[Literal::POS][v].size(); i++)
        onFalseBody(*(occ[Literal::POS][v][i]));
    for(int i = 0; i < occ[Literal::DNEG][v].size(); i++)
        onFalseBody(*(occ[Literal::DNEG][v][i]));
}

void AspSolver::onTrueHead(vec<Literal>& rule, Literal headAtom) {
    if(rule.size() == 0) return;
    
    for(int i = 0; i < rule.size(); i++) {
        if(rule[i] == headAtom) continue;
        if(rule[i].type == Literal::SHEAD) {
            rule[i].type = Literal::UHEAD;
            vec<vec<Literal>*>& o = occ[Literal::SHEAD][headAtom.id];
            int k = 0;
            for(int j = 0; j < o.size(); j++) {
                if(o[j] == &rule) continue;
                o[k++] = o[j];
            }
            o.shrink(k-o.size());
            if(o.size() == 0 || (o.size() == 1 && value(rule[i].id) == l_True && !supported[i])) supportInference.push(rule[i].id);
            occ[Literal::UHEAD][headAtom.id].push(&rule);
        }
    }
}

void AspSolver::onFalseHead(vec<Literal>& rule, Literal headAtom) {
    onTrueBody(rule, headAtom);
}

void AspSolver::onTrueBody(vec<Literal>& rule, Literal bodyAtom) {
    if(rule.size() == 0) return;

    bool trueBody = true;
    bool hasMustBeTrueBody = false;
    int headAtoms = 0;
    for(int i = 0; i < rule.size(); i++) {
        if(rule[i] == bodyAtom) {
            if(bodyAtom.type != Literal::POS || supported[bodyAtom.id]) {
                rule[i] = rule[rule.size()-1];
                rule.shrink(1);
            }
            else hasMustBeTrueBody = true;
            continue;
        }
        
        if(rule[i].head()) headAtoms++;
        else if(value(rule[i].id) == l_Undef) trueBody = false;
        else if(value(rule[i].id) == l_True) {
            if(rule[i].type == Literal::POS) hasMustBeTrueBody = true;
            else if(rule[i].type != Literal::DNEG) trueBody = false;
        }
        else if(rule[i].type == Literal::POS || rule[i].type == Literal::DNEG) trueBody = false;
    }
    
    if(headAtoms >= 2) return;
    if(headAtoms == 0) {
        vec<Lit> lits;
        for(int i = 0; i < rule.size(); i++) lits.push(rule[i].toLit());
        addClause(lits);
        rule.clear();
        return;
    }
    if(!trueBody) {
        for(int i = 0; i < rule.size(); i++) {
            if(rule[i].head()) {
                if(value(rule[i].id) != l_True) addClause(mkLit(rule[i].id));
                if(!hasMustBeTrueBody) supported[rule[i].id] = true;
                break;
            }
        }
    }
}

void AspSolver::onFalseBody(vec<Literal>& rule) {
    for(int i = 0; i < rule.size(); i++) {
        if(rule[i].type != Literal::SHEAD) continue;
        trace(asp_pre, 15, "Remove supporting rule for atom " << getName(rule[i].id));
        vec<vec<Literal>*>& heads = occ[Literal::SHEAD][rule[i].id];
        for(int j = 0; j < heads.size(); j++) {
            if(heads[j] != &rule) continue;
            heads[j] = heads[heads.size()-1];
        }
        heads.shrink(1);
        if(heads.size() == 0 || (heads.size() == 1 && value(rule[i].id) == l_True && !supported[rule[i].id])) supportInference.push(rule[i].id);
    }
    rule.clear();
}

void AspSolver::finalPropagation() {
    trace(asp_pre, 1, "Start finalPropagation()...");
    assert(nextToPropagate == nextToPropagateByUnit());
    
    supportInference.clear();
    for(int i = 0; i < nVars(); i++) {
        if(supported[i]) continue;
        if(value(i) == l_False) continue;
        vec<vec<Literal>*>& heads = occ[Literal::SHEAD][i];
        if(heads.size() == 0 || (heads.size() == 1 && value(i) == l_True)) supportInference.push(i);
    }
    
    while(supportInference.size() > 0) {
        Var v = supportInference[supportInference.size()-1];
        supportInference.shrink(1);
        
        vec<vec<Literal>*>& heads = occ[Literal::SHEAD][v];
        if(heads.size() == 0) {
            trace(asp_pre, 8, "Atom " << getName(v) << " has no support.");
            addClause(mkLit(v, true));
            continue;
        }
        
        assert(heads.size() == 1 && value(v) == l_True && !supported[v]);
        trace(asp_pre, 8, "Atom " << getName(v) << " is true and has a last supporting rule.");
        supported[v] = true;
        vec<Literal>& rule = *heads[0];
        for(int i = 0; i < rule.size(); i++) {
            if(rule[i] == Literal(v, Literal::SHEAD)) continue;
            addClause(mkLit(rule[i].id, rule[i].type != Literal::POS && rule[i].type != Literal::DNEG));
            if(rule[i].type == Literal::POS) supported[v] = false;
        }
        
        propagate();
    }
}

void AspSolver::processComponents() {
    if(!ok) return;
    
    trace(asp_pre, 1, "Start processComponents()...");
    DepGraph dg;
    for(int i = 0; i < nVars(); i++) {
        dg.add(i);
        if(value(i) == l_False) continue;
        if(supported[i]) continue;

        for(int o = Literal::SHEAD; o <= Literal::UHEAD; o++) {
            vec<vec<Literal>*>& heads = occ[o][i];
            for(int j = 0; j < heads.size(); j++) {
                vec<Literal>& rule = *heads[j];
                for(int k = 0; k < rule.size(); k++) {
                    if(rule[k].type == Literal::POS) {
                        assert(i != rule[k].id);
                        dg.add(i, rule[k].id);
                    }
                }
            }
        }
    }
    
    vec<int> atom2comp(nVars());
    vec<vec<int> > components;
    bool tight;
    dg.sccs(atom2comp, components, tight);
    assert(components.size() > 0);
    vec<Lit> supp;
    for(int i = 0; i < components.size(); i++) {
        vec<int>& c = components[i];
        if(c.size() == 1) { 
            if(value(i) == l_False) continue;
            if(supported[i]) continue;
            completion(c[0], supp);
            continue;
        }
        
        trace(asp_pre, 2, "Processing recursive component " << i);
        for(int j = 0; j < c.size(); j++) {
            completion(c[j], supp);
            setFrozen(c[j], true);
            
            unfounded.push(c[j]);
            vec<vec<Literal>*>& heads = occ[Literal::SHEAD][c[j]];
            for(int k = 0; k < heads.size(); k++) {
                vec<Literal>* rule = new vec<Literal>();
                heads[k]->copyTo(*rule);
                recursiveRules.push(rule);
                body.push(supp.size() == 0 ? mkLit(c[j]) : supp[k]);
                recBody.push();
                isBodyOf[sign(body[body.size()-1])][var(body[body.size()-1])].push(recursiveRules.size()-1);
                withoutSourcePointer.push(0);
                possibleSourcePointerOf.push(c[j]);
                possibleSourcePointers[c[j]].push(recursiveRules.size()-1);
                for(int idx = 0; idx < rule->size(); idx++) {
                    Literal lit = rule->operator[](idx);
                    if(lit == Literal(c[j], Literal::SHEAD)) continue;
                    if(lit == Literal(c[j], Literal::UHEAD)) continue;
                    inBody[lit.type == Literal::POS || lit.type ==  Literal::DNEG ? 1 : 0][lit.id].push(recursiveRules.size()-1);
                    if(lit.type == Literal::POS && atom2comp[lit.id] == i) {
                        inRecBody[lit.id].push(recursiveRules.size()-1); 
                        withoutSourcePointer[withoutSourcePointer.size()-1]++;
                        recBody[recBody.size()-1].push(lit.id);
                    }
                }
            }
        }
        if(!findSourcePointers()) ok = false;
    }
}

void AspSolver::completion(Var i, vec<Lit>& supp) {
    trace(asp_pre, 3, "Completion of atom " << i << "/" << getName(i));
    supp.clear();
    assert(value(i) != l_False);
    assert(!supported[i]);
    vec<vec<Literal>*>& heads = occ[Literal::SHEAD][i];
    assert_msg(heads.size() > 0, "Atom " << getName(i) << " has " << heads.size() << " definitions.");
    trace(asp_pre, 5, "Processing atom " << i << " with " << heads.size() << " possibly supporting rules.");
    if(heads.size() == 1) {
        if(value(i) == l_True) return;
        assert(!supported[i]);
        assert(value(i) == l_Undef);
        vec<Literal>& rule = *heads[0];
        vec<Lit> lits;
        for(int j = 0; j < rule.size(); j++) {
            lits.push(rule[j].toLit());
            if(rule[j] == Literal(i, Literal::SHEAD)) continue;
            addClause(~mkLit(i), ~rule[j].toLit());
        }
        addClause(lits);
        return;
    }
    
    for(int j = 0; j < heads.size(); j++) {
        vec<Literal>& rule = *heads[j];
        assert(rule.size() >= 2);
        Lit aux;
        if(rule.size() == 2)
            aux = rule[0] == Literal(i, Literal::SHEAD) ? ~rule[1].toLit() : ~rule[0].toLit();
        else {
            aux = mkLit(nVars());
            newVar();
            addClause(~aux, mkLit(i));
            vec<Lit> lits;
            lits.push(aux);
            for(int k = 0; k < rule.size(); k++) {
                if(rule[k] == Literal(i, Literal::SHEAD)) continue;
                lits.push(rule[k].toLit());
                addClause(~aux, ~rule[k].toLit());
            }
            addClause(lits);
            assert(nextToPropagate == nextToPropagateByUnit());
        }
        supp.push(aux);
    }
    supp.push(~mkLit(i));
    addClause(supp);
    supp.pop();
}

CRef AspSolver::morePropagate() {
    CRef ret = MaxSatSolver::morePropagate();
    if(ret != CRef_Undef) return ret;
    while(nextToPropagate < nextToPropagateByUnit()) {
        ret = morePropagate(mainTrail(nextToPropagate++));
        if(ret != CRef_Undef) return ret;
    }
    return ret;
}

CRef AspSolver::morePropagate(Lit lit) {
    trace(asp, 10, "Propagating " << lit << "@" << level(var(lit)) << " (atom " << var(lit) << "/" << getName(var(lit)) << ")");

    tagCalls++;
    assert(unfounded.size() == 0);
    
    vec<int>& in = isBodyOf[1-sign(lit)][var(lit)];
    for(int i = 0; i < in.size(); i++) {
        int ruleIdx = in[i];
        int head = possibleSourcePointerOf[ruleIdx];
        if(value(mkLit(head)) == l_False) continue;
        if(sourcePointer[head] == ruleIdx) {
            int last = unfounded.size();
            tag[head] = tagCalls;
            unfounded.push(head);
            do{
                int a = unfounded[last++];
                trace(asp, 15, "Atom " << a << "/" << getName(a) << " unfounded");
                for(int k = 0; k < inRecBody[a].size(); k++) {
                    if(value(mkLit(possibleSourcePointerOf[inRecBody[a][k]])) == l_False) continue;
                    withoutSourcePointer[inRecBody[a][k]]++;
                    if(tag[possibleSourcePointerOf[inRecBody[a][k]]] != tagCalls) {
                        tag[possibleSourcePointerOf[inRecBody[a][k]]] = tagCalls;
                        unfounded.push(possibleSourcePointerOf[inRecBody[a][k]]);
                    }
                }
            }while(last < unfounded.size());
        }
    }
    
//    vec<int>& in = inBody[1-sign(lit)][var(lit)];
//    for(int i = 0; i < in.size(); i++) {
//        cout << i << endl;
//        int ruleIdx = in[i];
//        int head = possibleSourcePointerOf[ruleIdx];
//        if(sourcePointer[head] == ruleIdx) {
//            int last = unfounded.size();
//            unfounded.push(head);
//            trace(asp, 15, "Atom " << head << "/" << getName(head) << " unfounded");
//            do{
//                int a = unfounded[last];
//                tag[a] = tagCalls;
//                for(int k = 0; k < inRecBody[a].size(); k++) {
//                    withoutSourcePointer[inRecBody[a][k]]++;
//                    if(tag[possibleSourcePointerOf[inRecBody[a][k]]] != tagCalls)
//                        ;//unfounded.push(tag[possibleSourcePointerOf[inRecBody[a][k]]]);
//                }
//            }while(++last < unfounded.size());
//        }
//    }
    if(!findSourcePointers()) return CRef_MoreConflict;

    return CRef_Undef;
}

bool AspSolver::findSourcePointers() {
    trace(asp, 10, "Find source pointers for " << unfounded);
    tagCalls++;
    for(int i = 0; i < unfounded.size(); i++) {
        int atom = unfounded[i];
        assert(value(atom) != l_False);
        if(tag[atom] == tagCalls) continue;
        vec<int>& poss = possibleSourcePointers[atom];
        for(int j = 0; j < poss.size(); j++) {
            if(value(body[poss[j]]) == l_False) continue;
            if(withoutSourcePointer[poss[j]] != 0) continue;
            
            vec<int> stackA;
            vec<int> stackR;
            tag[atom] = tagCalls;
            stackA.push(atom);
            stackR.push(poss[j]);
            
            do{
                int a = stackA[stackA.size()-1];
                sourcePointer[a] = stackR[stackR.size()-1];
                stackA.shrink(1);
                stackR.shrink(1);
                trace(asp, 15, "New source pointer for atom " << a << "/" << getName(a) << ": " << sourcePointer[a]);
                for(int k = 0; k < inRecBody[a].size(); k++) {
                    if(value(possibleSourcePointerOf[inRecBody[a][k]]) == l_False) continue;
                    if(--withoutSourcePointer[inRecBody[a][k]] == 0 && tag[possibleSourcePointerOf[inRecBody[a][k]]] != tagCalls) {
                        tag[possibleSourcePointerOf[inRecBody[a][k]]] = tagCalls;
                        stackA.push(possibleSourcePointerOf[inRecBody[a][k]]);
                        stackR.push(inRecBody[a][k]);
                    }
                }
            }while(stackA.size() > 0);

            break;
        }
    }
    
    // remove founded atoms
    int j = 0;
    for(int i = 0; i < unfounded.size(); i++) {
        int atom = unfounded[i];
        if(tag[atom] == tagCalls) continue;
        unfounded[j++] = unfounded[i];
        withoutSourcePointer[sourcePointer[atom]] = 0;
    }
    unfounded.shrink_(unfounded.size()-j);
    
    // look for inconsistency
    for(int i = 0; i < unfounded.size(); i++) {
        if(value(unfounded[i]) == l_Undef) continue;
        assert(value(unfounded[i]) == l_True);
        trace(asp, 15, "True atom is unfounded: " << unfounded[i] << "/" << getName(unfounded[i]));
        int tmp = unfounded[i];
        unfounded[i] = unfounded[0];
        unfounded[0] = tmp;
        break;
    }
    
    // identify small unfounded sets and reasons
    tagCalls++;
    for(int i = 0; i < unfounded.size(); i++) {
        if(tag[unfounded[i]] == tagCalls) continue;
        tag[unfounded[i]] = tagCalls;
        
        moreReasonWFVec.push();
        vec<Lit>& reasons = moreReasonWFVec.last();
        vec<int> uSet;
        uSet.push(unfounded[i]);
        int processed = 0;
        while(processed < uSet.size()) {
            int atom = uSet[processed++];
            vec<int>& poss = possibleSourcePointers[atom];
            for(int j = 0; j < poss.size(); j++) {
                if(value(body[poss[j]]) == l_False) {
                    reasons.push(body[poss[j]]);
                    continue;
                }
                vec<int>& rec = recBody[poss[j]];
                assert(rec.size() > 0);
                int toBeAdd = -1;
                for(int i = 0; i < rec.size(); i++) {
                    if(tag[rec[i]] == tagCalls) {
                        if(value(rec[i]) != l_False) continue;
                        reasons.push(mkLit(rec[i]));
                        toBeAdd = -1;
                        break;
                    }
                    toBeAdd = rec[i];
                }
                if(toBeAdd != -1) {
                    uSet.push(toBeAdd);
                    tag[toBeAdd] = tagCalls;
                }
            }
        }
        trace(asp, 15, "Unfounded set (atoms) " << uSet << " with reasons (literals) " << reasons);
        for(int j = 0; j < uSet.size(); j++) {
            int v = uSet[j];
            trace(asp, 20, "Inferring atom " << v << "/" << getName(v) << "@" << decisionLevel() << " from unfounded set");
            assert(moreReasonWF[v] == -1);
            moreReasonWF[v] = moreReasonWFVec.size()-1;
            moreReasonVars.push(v);
            if(value(v) == l_True) {
                moreConflictLit = ~mkLit(v);
                moreConflictWF = moreReasonWFVec.size() - 1;
                unfounded.clear();
                return false;
            }
            else uncheckedEnqueue(~mkLit(v));
        }
    }
    
    unfounded.clear();
    return true;
}

bool AspSolver::moreReason(Lit lit, vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC) {
    if(PseudoBooleanSolver::moreReason(lit, out_learnt, selectors, pathC)) return true;
    if(moreReasonWF[var(lit)] != -1) { _moreReasonWF(lit, out_learnt, selectors, pathC); return true; }
    return false;
}

void AspSolver::_moreReasonWF(Lit lit, vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC) {
    assert(decisionLevel() != 0);
    assert(reason(var(lit)) == CRef_Undef);
    vec<Lit>& reasons = moreReasonWFVec[moreReasonWF[var(lit)]];
    for(int i = 0; i < reasons.size(); i++) {
        Lit q = reasons[i];
        assert(value(q) == l_False);
        assert(level(var(q)) <= level(var(lit)));
        
        if(seen[var(q)]) continue;
        if(level(var(q)) == 0) continue;
        
        if(!isSelector(var(q)))
            varBumpActivity(var(q));
        
        seen[var(q)] = 1;
        
        if(level(var(q)) >= decisionLevel()) {
            pathC++;
            // UPDATEVARACTIVITY trick (see competition'09 companion paper)
            if(!isSelector(var(q)) && (reason(var(q)) != CRef_Undef) && ca[reason(var(q))].learnt())
                lastDecisionLevel.push(q);
        }
        else {
            if(isSelector(var(q))) {
                assert(value(q) == l_False);
                selectors.push(q);
            }
            else 
                out_learnt.push(q);
        }
    }        
}

bool AspSolver::moreReason(Lit lit) {
    if(moreReasonWF[var(lit)] != -1) { _moreReasonWF(lit); return true; }
    return false;
}

void AspSolver::_moreReasonWF(Lit lit) {
    assert(decisionLevel() != 0);
    assert(reason(var(lit)) == CRef_Undef);
    vec<Lit>& reasons = moreReasonWFVec[moreReasonWF[var(lit)]];
    for(int i = 0; i < reasons.size(); i++) {
        Lit l = reasons[i];
        assert(value(l) == l_False);
        assert(level(var(l)) <= level(var(lit)));
        if(level(var(l)) == 0) continue;
        seen[var(l)] = 1;
    }
}

bool AspSolver::moreConflict(vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC) {
    if(moreConflictWF != -1) { _moreConflictWF(out_learnt, selectors, pathC); return true; }
    return false;
}
    
void AspSolver::_moreConflictWF(vec<Lit>& out_learnt, vec<Lit>& selectors, int& pathC) {
    assert(decisionLevel() != 0);
    
    if(!seen[var(moreConflictLit)] && level(var(moreConflictLit)) > 0) {
        if(!isSelector(var(moreConflictLit))) varBumpActivity(var(moreConflictLit));
        seen[var(moreConflictLit)] = 1;
        assert(level(var(moreConflictLit)) == decisionLevel());
        pathC++;
        // UPDATEVARACTIVITY trick (see competition'09 companion paper)
        if(!isSelector(var(moreConflictLit)) && (reason(var(moreConflictLit)) != CRef_Undef) && ca[reason(var(moreConflictLit))].learnt())
            lastDecisionLevel.push(moreConflictLit);
    }
    
    vec<Lit>& reasons = moreReasonWFVec[moreConflictWF];
    moreConflictWF = -1;
    for(int i = 0; i < reasons.size(); i++) {
        Lit q = reasons[i];
        assert(value(q) == l_False);
        
        if(seen[var(q)]) continue;
        if(level(var(q)) == 0) continue;
        
        if(!isSelector(var(q)))
            varBumpActivity(var(q));
        
        seen[var(q)] = 1;
        
        if(level(var(q)) >= decisionLevel()) {
            pathC++;
            // UPDATEVARACTIVITY trick (see competition'09 companion paper)
            if(!isSelector(var(q)) && (reason(var(q)) != CRef_Undef) && ca[reason(var(q))].learnt())
                lastDecisionLevel.push(q);
        }
        else {
            if(isSelector(var(q))) {
                assert(value(q) == l_False);
                selectors.push(q);
            }
            else 
                out_learnt.push(q);
        }
    }
}

void AspSolver::clearParsingStructures() {
    for(int i = 0; i < program.size(); i++) delete program[i];
    { vec<vec<Literal>*> tmp; program.moveTo(tmp); }
    
    for(int i = 0; i < 5; i++) { vec<vec<vec<Literal>*> > tmp; occ[i].moveTo(tmp); }
    { vec<bool> tmp; supported.moveTo(tmp); }
    
    { vec<Var> tmp; supportInference.moveTo(tmp); }
}

void AspSolver::onCancel() {
    MaxSatSolver::onCancel();
    
    trace(asp, 2, "Cancel until level " << decisionLevel());
    while(nextToPropagate > nextToPropagateByUnit()) { 
        int v = var(mainTrail(--nextToPropagate));
        propagated[v] = false;
    }

    while(moreReasonVars.size() > 0) {
        Var v = moreReasonVars.last();
        if(level(v) <= decisionLevel()) break;
        moreReasonVars.pop();
        if(moreReasonWF[v] < moreReasonWFVec.size()-1) moreReasonWFVec.pop();
        moreReasonWF[v] = -1;
    }
}


} // namespace aspino

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
    propagated = 0;
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
            occ[lit.type][lit.id]->push(rule);
        }
    }
    
    vec<Lit> lits;
    for(int i = 0; i < rule->size(); i++) {
        Literal lit = rule->operator[](i);
        lits.push(lit.toLit());
    }
    addClause(lits);
}

void AspSolver::parse(gzFile in_) {
    trace(asp_pre, 1, "Start parsing...");
    Glucose::StreamBuffer in(in_);
    
    for(;;) {
        skipWhitespace(in);
        if(*in == EOF) cerr << "PARSE ERROR! Unexpected EOF." << endl, exit(3);
        int type = parseInt(in);
        if(type == 0) break;
        if(type == 1) parseNormalRule(in);
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
    completion();
    
    trace(asp_pre, 1, "All done.");
}

int AspSolver::getId(int input_id) {
    assert(idmap.size() == occ[0].size() && idmap.size() == occ[1].size() && idmap.size() == occ[2].size() && idmap.size() == occ[3].size() && idmap.size() == occ[4].size());
    while(input_id >= idmap.size()) {
        idmap.push(-1);
        for(int i = Literal::SHEAD; i <= Literal::DNEG; i++) occ[i].push(new vec<vec<Literal>*>());
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
}

lbool AspSolver::solve() {
    lbool res = MaxSatSolver::solve();
    if(res == l_True) {
        nModels++;
        printModel();
    }
    return res;
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
    while(propagated < static_cast<unsigned>(trail.size())) {
        Var v = var(trail[propagated++]);
        if(value(v) == l_True) propagateTrue(v);
        else propagateFalse(v);
    }
}

void AspSolver::propagateTrue(Var v) {
    trace(asp_pre, 10, "Propagate true: " << getName(v));
    
    for(int i = 0; i < occ[Literal::SHEAD][v]->size(); i++)
        onTrueHead(*(occ[Literal::SHEAD][v]->operator[](i)), Literal(v, Literal::SHEAD));
    for(int i = 0; i < occ[Literal::UHEAD][v]->size(); i++)
        onTrueHead(*(occ[Literal::UHEAD][v]->operator[](i)), Literal(v, Literal::UHEAD));
        
    for(int i = 0; i < occ[Literal::POS][v]->size(); i++)
        onTrueBody(*(occ[Literal::POS][v]->operator[](i)), Literal(v, Literal::POS));
    for(int i = 0; i < occ[Literal::DNEG][v]->size(); i++)
        onTrueBody(*(occ[Literal::DNEG][v]->operator[](i)), Literal(v, Literal::DNEG));

    for(int i = 0; i < occ[Literal::NEG][v]->size(); i++)
        onFalseBody(*(occ[Literal::NEG][v]->operator[](i)));
}

void AspSolver::propagateFalse(Var v) {
    trace(asp_pre, 10, "Propagate false: " << getName(v));

    for(int i = 0; i < occ[Literal::SHEAD][v]->size(); i++)
        onFalseHead(*(occ[Literal::SHEAD][v]->operator[](i)), Literal(v, Literal::SHEAD));
    for(int i = 0; i < occ[Literal::UHEAD][v]->size(); i++)
        onFalseHead(*(occ[Literal::UHEAD][v]->operator[](i)), Literal(v, Literal::UHEAD));


    for(int i = 0; i < occ[Literal::NEG][v]->size(); i++)
        onTrueBody(*(occ[Literal::NEG][v]->operator[](i)), Literal(v, Literal::NEG));

    for(int i = 0; i < occ[Literal::POS][v]->size(); i++)
        onFalseBody(*(occ[Literal::POS][v]->operator[](i)));
    for(int i = 0; i < occ[Literal::DNEG][v]->size(); i++)
        onFalseBody(*(occ[Literal::DNEG][v]->operator[](i)));
}

void AspSolver::onTrueHead(vec<Literal>& rule, Literal headAtom) {
    if(rule.size() == 0) return;
    
    for(int i = 0; i < rule.size(); i++) {
        if(rule[i] == headAtom) continue;
        if(rule[i].type == Literal::SHEAD) {
            rule[i].type = Literal::UHEAD;
            vec<vec<Literal>*>& o = *occ[Literal::SHEAD][headAtom.id];
            int k = 0;
            for(int j = 0; j < o.size(); j++) {
                if(o[j] == &rule) continue;
                o[k++] = o[j];
            }
            o.shrink(k-o.size());
            if(o.size() == 0 || (o.size() == 1 && value(rule[i].id) == l_True && !supported[i])) supportInference.push(rule[i].id);
            occ[Literal::UHEAD][headAtom.id]->push(&rule);
        }
    }
}

void AspSolver::onFalseHead(vec<Literal>& rule, Literal headAtom) {
    onTrueBody(rule, headAtom);
}

void AspSolver::onTrueBody(vec<Literal>& rule, Literal bodyAtom) {
    if(rule.size() == 0) return;
    
    bool undefinedBody = false;
    bool hasMustBeTrueBody = false;
    int headAtoms = 0;
    for(int i = 0; i < rule.size(); i++) {
        if(rule[i] == bodyAtom) {
            if(bodyAtom.type != Literal::POS || supported[bodyAtom.id]) {
                rule[i] = rule[rule.size()-1];
                rule.shrink(1);
            }
            continue;
        }
        
        if(value(rule[i].id) == l_Undef && rule[i].body()) undefinedBody = true;
        if(value(rule[i].id) == l_True && rule[i].body()) hasMustBeTrueBody = true;
        if(rule[i].head()) headAtoms++;
    }
    
    if(headAtoms >= 2) return;
    if(headAtoms == 0) {
//        vec<Lit> lits;
//        for(int i = 0; i < rule.size(); i++) lits.push(rule[i].toLit());
//        addClause(lits);
        rule.clear();
        return;
    }
    if(!undefinedBody) {
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
        vec<vec<Literal>*>& heads = *occ[Literal::SHEAD][rule[i].id];
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
    assert(propagated == static_cast<unsigned>(trail.size()));
    
    supportInference.clear();
    for(int i = 0; i < nVars(); i++) {
        if(supported[i]) continue;
        if(value(i) == l_False) continue;
        vec<vec<Literal>*>& heads = *occ[Literal::SHEAD][i];
        if(heads.size() == 0 || (heads.size() == 1 && value(i) == l_True)) supportInference.push(i);
    }
    
    while(supportInference.size() > 0) {
        Var v = supportInference[supportInference.size()-1];
        supportInference.shrink(1);
        
        vec<vec<Literal>*>& heads = *occ[Literal::SHEAD][v];
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
        }
        
        propagate();
    }
}

void AspSolver::completion() {
    trace(asp_pre, 1, "Start completion()...");
    if(!ok) return;
    for(int i = 0; i < nVars(); i++) {
        if(value(i) == l_False) continue;
        if(supported[i]) continue;
        
        vec<vec<Literal>*>& heads = *occ[Literal::SHEAD][i];
        assert_msg(heads.size() > 0, "Atom " << getName(i) << " has " << heads.size() << " definitions.");
        if(heads.size() == 1) {
            assert(!supported[i]);
            assert(value(i) == l_Undef);
            for(int j = 0; j < heads[0]->size(); j++) {
                if(heads[0]->operator[](j).id != i) {
                    addClause(mkLit(i, true), ~heads[0]->operator[](j).toLit());
                }
            }
        }
    }
    
    
    
}

} // namespace aspino

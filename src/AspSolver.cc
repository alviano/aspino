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
    while(propagated < static_cast<unsigned>(trail.size())) {
        Var v = var(trail[propagated++]);
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
    assert(propagated == static_cast<unsigned>(trail.size()));
    
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
        }
        
        propagate();
    }
}

void AspSolver::processComponents() {
    if(!ok) return;
    
    trace(asp_pre, 1, "Start processComponents()...");
    DepGraph dg;
    for(int i = 0; i < nVars(); i++) {
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
    for(int i = 0; i < components.size(); i++) {
        vec<int>& c = components[i];
        if(c.size() == 1) { 
            if(value(i) == l_False) continue;
            if(supported[i]) continue;
            completion(c[0]);
            continue;
        }
        
        trace(asp_pre, 2, "Processing recursive component " << i);
        for(int j = 0; j < c.size(); j++) {
            completion(c[j]);
            
            unfounded.push(c[j]);
            vec<vec<Literal>*>& heads = occ[Literal::SHEAD][j];
            for(int k = 0; k < heads.size(); k++) {
                vec<Literal>* rule = new vec<Literal>();
                heads[k]->moveTo(*rule);
                recursiveRules.push(rule);
                withoutSourcePointer.push(0);
                possibleSourcePointerOf.push(c[j]);
                possibleSourcePointers[c[j]].push(recursiveRules.size()-1);
                for(int idx = 0; idx < rule->size(); idx++) {
                    Literal lit = rule->operator[](idx);
                    if(lit == Literal(c[j], Literal::SHEAD) || lit == Literal(c[j], Literal::UHEAD)) continue;
                    inBody[lit.type == Literal::POS || lit.type ==  Literal::DNEG ? 1 : 0][lit.id].push(recursiveRules.size()-1);
                    if(lit.type == Literal::POS && atom2comp[lit.id] == i) { 
                        inRecBody[lit.id].push(recursiveRules.size()-1); 
                        withoutSourcePointer[withoutSourcePointer.size()-1]++;
                    }
                }
            }
        }
        findSourcePointers();
    }
}

void AspSolver::completion(Var i) {
    trace(asp_pre, 3, "Completion of atom " << i);
    assert(value(i) != l_False);
    assert(!supported[i]);
    vec<vec<Literal>*>& heads = occ[Literal::SHEAD][i];
    assert_msg(heads.size() > 0, "Atom " << getName(i) << " has " << heads.size() << " definitions.");
    trace(asp_pre, 5, "Processing atom " << i << " with " << heads.size() << " possibly supporting rules.");
    if(heads.size() == 1) {
        vec<Literal>& rule = *heads[0];
        assert(!supported[i]);
        assert(value(i) == l_Undef);
        for(int j = 0; j < rule.size(); j++) {
            if(rule[j] != Literal(i, Literal::SHEAD)) {
                addClause(~mkLit(i), ~rule[j].toLit());
            }
        }
        return;
    }
    
    vec<Lit> supp;
    supp.push(~mkLit(i));
    for(int j = 0; j < heads.size(); j++) {
        vec<Literal>& rule = *heads[j];
        assert(rule.size() >= 2);
        Lit aux;
        if(rule.size() == 2)
            aux = rule[0] == Literal(i, Literal::SHEAD) ? ~rule[1].toLit() : ~rule[0].toLit();
        else {
            aux = ~mkLit(nVars());
            newVar();
            vec<Lit> lits;
            lits.push(~aux);
            for(int k = 0; k < rule.size(); i++) {
                if(rule[k] == Literal(i, Literal::SHEAD)) continue;
                lits.push(rule[k].toLit());
                addClause(aux, ~rule[k].toLit());
            }
            addClause(lits);
            assert(propagated == static_cast<unsigned>(trail.size()));
        }
        supp.push(aux);
    }
}

void AspSolver::findSourcePointers() {
    tagCalls++;
    for(int i = 0; i < unfounded.size(); i++) {
        int atom = unfounded[i];
        if(tag[atom] == tagCalls) continue;
        vec<int>& poss = possibleSourcePointers[atom];
        for(int j = 0; j < poss.size(); j++) {
            if(withoutSourcePointer[poss[j]] == 0) {
                vec<int> stackA;
                vec<int> stackR;
                stackA.push(atom);
                stackR.push(poss[j]);
                
                do{
                    int a = stackA[stackA.size()-1];
                    sourcePointer[a] = stackR[stackR.size()-1];
                    tag[a] = tagCalls;
                    stackA.shrink(1);
                    stackR.shrink(1);
                    trace(asp_pre, 10, "New source pointer for atom " << a << ": " << sourcePointer[a]);
                    for(int k = 0; k < inRecBody[a].size(); k++) {
                        if(--withoutSourcePointer[inRecBody[a][k]] == 0 && tag[possibleSourcePointerOf[inRecBody[a][k]]] != tagCalls) {
                            stackA.push(possibleSourcePointerOf[inRecBody[a][k]]);
                            stackR.push(inRecBody[a][k]);
                        }
                    }
                }while(stackA.size() > 0);

                break;
            }
        }
    }
    
    for(int i = 0; i < unfounded.size(); i++) {
        int atom = unfounded[i];
        if(tag[atom] == tagCalls) continue;
        cout << "Unfounded: " << atom << endl; // TODO
        withoutSourcePointer[sourcePointer[atom]] = 0;
    }
}

void AspSolver::clearParsingStructures() {
    for(int i = 0; i < program.size(); i++) delete program[i];
    { vec<vec<Literal>*> tmp; program.moveTo(tmp); }
    
    for(int i = 0; i < 5; i++) { vec<vec<vec<Literal>*> > tmp; occ[i].moveTo(tmp); }
    { vec<bool> tmp; supported.moveTo(tmp); }
    
    { vec<Var> tmp; supportInference.moveTo(tmp); }
}

} // namespace aspino

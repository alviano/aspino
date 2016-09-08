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

#include "TGDsSolver.h"

#include "utils/ParseUtils.h"
#include "core/Dimacs.h"

#include "main.h"

namespace aspino {

TGDsSolver::Listener TGDsSolver::defaultListener;
    
static const int buffer_size = 1048576;
    
class LineStream {
private:    
    gzFile        in;
    int           pos;

public:
    char buf[buffer_size];
    explicit LineStream(gzFile i) : in(i), pos(0) {}

    bool readline() { 
        pos = 0;
        buf[0] = '\0';
        return gzgets(in, buf, sizeof(buf)) != NULL;
    }
    int  operator *  () const { return buf[pos]; }
    void operator ++ ()       { pos++; }
    int  position    () const { return pos; }
};

TGDsSolver::TGDsSolver() : SatSolver() {
    listener = &TGDsSolver::defaultListener;
}

void TGDsSolver::parse(gzFile in_) {
    cout << "c Interpreter is online! Echo mode on." << endl;
    cout << "c Press CTRL+D to flush the input stream. Two times to termiante." << endl;
    cout << "c" << endl;
    cout << "c Start your session with 'p query' followed by a list of literals terminated by 0." << endl;
    cout << "c Add one clause per line in DIMACS format: list of literals terminated by 0." << endl;
    cout << "c Start refutation with 'r' followed by a list of assumption literals terminated by 0." << endl;
    
    LineStream in(in_);
    
    vec<Lit> lits;
    while(in.readline()) {
        cout << "c echo: " << in.buf;
        cout.flush();
        Glucose::skipWhitespace(in);
        assert(*in != EOF);
        if(*in == 'p') {
            if(Glucose::eagerMatch(in, "p query")) {
                readClause(in, *this, candidates);
            }
            else {
                cerr << "PARSE ERROR! Unexpected char: " << static_cast<char>(*in) << endl, exit(3);
            }
        }
        else if(*in == 'r') {
            ++in;
            readClause(in, *this, assumptions);
            cout << "c run: candidates=" << candidates << "; assumptions=" << assumptions << endl;
            disprove();
            cout << "c done: candidates=" << candidates << endl;
            if(candidates.size() == 0) {
                cout << "." << endl;
                break;
            }
        }
        else if(*in == 'c') continue;
        else {
            readClause(in, *this, lits);
            addClause(lits); 
        }
    }
}

void TGDsSolver::disprove() {
    assert(decisionLevel() == 0);
    
    cout << "Number of candidates: " << candidates.size() << endl;
    cout << "Number of assumptions: " << assumptions.size() << endl;
    
//    cout << "VARS: " << nVars() << endl;
//    for(int i = 0; i < clauses.size(); i++) {
//        Clause& clause = ca[clauses[i]];
//        cout << clause << endl;
//    }
//    cout << "ASS: " << assumptions << endl;
//    cout << "CAND: " << candidates << endl;
    
    vec<lbool> flagCandidates(candidates.size(), l_Undef);
    
    for(int i = 0; i < candidates.size(); i++) {
        if(value(candidates[i]) != l_Undef && level(var(candidates[i])) == 0) continue;
        if(flagCandidates[i] != l_Undef) continue;
        
        assumptions.push(~candidates[i]);
        solve_();
        if(status == l_True) {
            for(int j = 0; j < candidates.size(); j++)
                if(value(candidates[j]) != l_True) {
                    flagCandidates[j] = l_False;
                }
        }
        assumptions.pop();
        cancelUntil(assumptions.size());
    }
    
    int j = 0;
    for(int i = 0; i < candidates.size(); i++) {
        if(value(candidates[i]) != l_Undef && level(var(candidates[i])) == 0) flagCandidates[i] = value(candidates[i]);

        if(flagCandidates[i] == l_True) listener->onTrueCandidate(candidates[i], i);
        else if(flagCandidates[i] == l_False) listener->onFalseCandidate(candidates[i], i);
        else candidates[j++] = candidates[i];
    }
    candidates.shrink_(candidates.size()-j);
    
    cancelUntil(0);
}

void TGDsSolver::Listener::onTrueCandidate(Lit candidate, int) {
    cout << "c answer: " << candidate << " is true" << endl;
    cout << "t " << candidate << endl;
}

void TGDsSolver::Listener::onFalseCandidate(Lit candidate, int) {
    cout << "c answer: " << candidate << " is false" << endl;
    cout << "f " << candidate << endl;
}

} // namespace miniasp

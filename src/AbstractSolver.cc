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

#include "AbstractSolver.h"

namespace aspino {

ostream& operator<<(ostream& o, lbool l) {
    if(l == l_True) return o << "T";
    if(l == l_False) return o << "F";
    if(l == l_Undef) return o << "U";
    assert(0);
    return o;
}
    
ostream& operator<<(ostream& o, const Lit& l) {
    return o << (sign(l) ? "~" : "") << (var(l)+1);
}

} // namespace aspino

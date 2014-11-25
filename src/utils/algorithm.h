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

#ifndef __algorithm_h__
#define __algorithm_h__

#include "mtl/Vec.h"

namespace aspino {

template <class T>
void shuffle(vec<T>& v) {
    for(int i = 0; i < v.size(); i++) {
        int j = rand() % v.size();
        T tmp = v[i];
        v[i] = v[j];
        v[j] = tmp;
    }
}
    
} // namespace aspino

#endif

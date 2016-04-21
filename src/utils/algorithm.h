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
void shuffle(vec<T>& v, double random_seed) {
    for(int i = 0; i < v.size(); i++) {
        int j = Glucose::Solver::irand(random_seed, v.size());
        T tmp = v[i];
        v[i] = v[j];
        v[j] = tmp;
    }
}

template <class T>
void reverse(vec<T>& v) {
    for(int i = 0, j = v.size()-1; i < v.size() / 2; i++, j--) {
        T tmp = v[i];
        v[i] = v[j];
        v[j] = tmp;
    }
}

template <class T>
void quickSort(vec<T>& arr, int left, int right) {
    int i = left, j = right;
    T tmp;
    assert((left + right) / 2 >= 0);
    assert_msg((left + right) / 2 < arr.size(), "Accessing element " << (left + right) / 2 << " in array of size " << arr.size());
    int pivot = arr[(left + right) / 2];

    /* partition */
    while (i <= j) {
        while (arr[i] < pivot)
              i++;
        while (arr[j] > pivot)
              j--;
        if (i <= j) {
              tmp = arr[i];
              arr[i] = arr[j];
              arr[j] = tmp;
              i++;
              j--;
        }
    };

    /* recursion */
    if (left < j)
        quickSort(arr, left, j);
    if (i < right)
        quickSort(arr, i, right);
}

template <class T>
void sort(vec<T>& v) {
    if(v.size() > 1) quickSort(v, 0, v.size()-1);
}

} // namespace aspino

#endif

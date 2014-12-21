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

#ifndef __trace_h__
#define __trace_h__

#ifndef TRACE_ON
    #define trace(type, level, msg)
#else

#include <cstdio>
#include <string>

#include <utils/Options.h>

extern Glucose::IntOption option_trace_sat;
extern Glucose::IntOption option_trace_pbs;
extern Glucose::IntOption option_trace_maxsat;

#define trace(type, level, msg) \
    if(option_trace_##type >= level) {\
        std::cerr << "[" << #type << "]";\
        for(int __trace_i__ = 0; __trace_i__ < level; __trace_i__++) std::cerr << " ";\
        std::cerr << msg << std::endl;\
    }

#endif

#endif

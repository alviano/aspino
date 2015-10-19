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

#include "trace.h"

#ifdef TRACE_ON

Glucose::IntOption option_trace_sat("TRACE", "trace-sat", "Set trace level of SAT solving.", 0, Glucose::IntRange(0, INT32_MAX));
Glucose::IntOption option_trace_pbs("TRACE", "trace-pbs", "Set trace level of PBS solving.", 0, Glucose::IntRange(0, INT32_MAX));
Glucose::IntOption option_trace_maxsat("TRACE", "trace-maxsat", "Set trace level of MAXSAT solving.", 0, Glucose::IntRange(0, INT32_MAX));
Glucose::IntOption option_trace_asp("TRACE", "trace-asp", "Set trace level of ASP solving.", 0, Glucose::IntRange(0, INT32_MAX));
Glucose::IntOption option_trace_asp_pre("TRACE", "trace-asp-pre", "Set trace level of ASP solving (preprocessing).", 0, Glucose::IntRange(0, INT32_MAX));
Glucose::IntOption option_trace_fairsat("TRACE", "trace-fairsat", "Set trace level of FairSAT solving.", 0, Glucose::IntRange(0, INT32_MAX));

#endif
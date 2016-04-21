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


#include "main.h"

extern Glucose::EnumOption option_maxsat_strat;
extern Glucose::IntOption option_maxsat_tag;

int main(int argc, char** argv)
{
    premain();
    
//    Glucose::parseOptions(argc, argv, true);
    option_mode = "maxsat";
    option_maxsat_strat = "one-pmres";
    option_maxsat_tag = 8;

    return postmain(argc, argv);
}
/*
 *  Copyright (C) 2017  Mario Alviano (mario@alviano.net)
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

extern Glucose::IntOption option_circumscription_witnessess;

int main(int argc, char** argv)
{
    premain();
//    Glucose::parseOptions(argc, argv, true);
    option_mode = "circumscription";
    
    int j = 0;
    for(int i = 0; i < argc; i++) {
        argv[j] = argv[i];
        char* arg = argv[i];
        if(arg[0] == '-') arg++;
        if(arg[0] == '-') arg++;
        if(strcmp(arg, "help") == 0) {
            cout << "usage: " << argv[0] << " [flags] [input-file]\n\n"
                 << "OPTIONS:\n"
                 << "  -n            = <int32>  [   0 .. imax] (default: 1)\n"
                 << "  --circ-wit    = <int32>  [   0 .. imax] (default: 1)\n"
                 << "  --help        Print help message.\n";
            return 0;
        }
        if(strncmp(arg, "n=", 2) == 0) option_n = stoi(arg + 2);
        else if(strncmp(arg, "circ-wit=", 9) == 0) option_circumscription_witnessess = stoi(arg + 9);
        else if(arg != argv[i]) cerr << "ERROR! Unknown flag " << argv[i] << ". Use '--help' for help. " << endl, exit(3);
        else j++;
    }
    argc = j;

    return postmain(argc, argv);
}
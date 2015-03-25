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

#include "DepGraph.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>

#include <vector>

using std::vector;

namespace aspino {
    
class AdjacencyList: public boost::adjacency_list<>
{
};

DepGraph::DepGraph()
: arcs(*new AdjacencyList()) {
}

DepGraph::~DepGraph() {
    delete &arcs;
}

void DepGraph::add(int from, int to) {
    assert(from != to);
    boost::add_edge(from, to, arcs);
}

void DepGraph::sccs(vec<int>& atom2comp, vec<vec<int> >& components, bool& tight) {
    vector<int> sccs(boost::num_vertices(arcs)), discover_time(boost::num_vertices(arcs));
    vector<boost::default_color_type> color(boost::num_vertices(arcs));
    vector<boost::graph_traits<boost::adjacency_list<> >::vertex_descriptor> root(boost::num_vertices(arcs));
    int n = boost::strong_components(arcs, &sccs[0] , boost::root_map(&root[0]).color_map(&color[0]).discover_time_map(&discover_time[0]));

    assert(n > 0);
    components.clear();
    components.growTo(n);
    tight = true;
    
    for(int i = 0; i < static_cast<int>(sccs.size()); i++) {
        int scc = sccs[i];
        assert(scc < components.size());
    
        atom2comp[i] = scc;
        components[scc].push(i);
        if(components[scc].size() > 1) tight = false;
    }    
}

} // namespace aspino
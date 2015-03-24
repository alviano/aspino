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

void DepGraph::sccs(vector<vector<int> >& components, bool& tight) {
    vector<int> sccs(boost::num_vertices(arcs)), discover_time(boost::num_vertices(arcs));
    vector<boost::default_color_type> color(boost::num_vertices(arcs));
    vector<boost::graph_traits<boost::adjacency_list<> >::vertex_descriptor> root(boost::num_vertices(arcs));
    int n = boost::strong_components(arcs, &sccs[0] , boost::root_map(&root[0]).color_map(&color[0]).discover_time_map(&discover_time[0]));

    assert(n > 0);
    components.clear();
    components.resize(n);
    tight = true;
    
    for(vector<int>::size_type i = 0; i != sccs.size(); i++) {
        int scc = sccs[i];
        assert(scc < static_cast<int>(components.size()));
    
        components[scc].push_back(i);
        if(components[scc].size() > 1) tight = false;
    }    
}

} // namespace aspino
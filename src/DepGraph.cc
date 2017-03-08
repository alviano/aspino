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
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

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

void DepGraph::add(int node) {
    assert(node >= 0);
    boost::add_edge(node + 1, 0, arcs);
}

void DepGraph::add(int from, int to) {
    assert(from != to);
    assert(from >= 0);
    assert(to >= 0);
    boost::add_edge(from + 1, to + 1, arcs);
}

void DepGraph::sccs(vec<int>& atom2comp, vec<vec<int> >& components, bool& tight) {
    vector<int> sccs(boost::num_vertices(arcs)), discover_time(boost::num_vertices(arcs));
    vector<boost::default_color_type> color(boost::num_vertices(arcs));
    vector<boost::graph_traits<boost::adjacency_list<> >::vertex_descriptor> root(boost::num_vertices(arcs));
    int n = boost::strong_components(arcs, &sccs[0] , boost::root_map(&root[0]).color_map(&color[0]).discover_time_map(&discover_time[0]));

    tight = true;
    components.clear();
    if(n == 0) return;
    components.growTo(n-1);
    
    assert(sccs[0] == 0);
    for(int i = 1; i < static_cast<int>(sccs.size()); i++) {
        int scc = sccs[i] - 1;
        assert(scc >= 0);
        assert(scc < components.size());
    
        atom2comp[i-1] = scc;
        components[scc].push(i-1);
        if(components[scc].size() > 1) tight = false;
    }    
}


class WeightedAdjacencyList : public boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::property<boost::edge_weight_t, unsigned> > {};

UndirectedWeightedGraph::UndirectedWeightedGraph()
: arcs(*new WeightedAdjacencyList()) {
}

UndirectedWeightedGraph::~UndirectedWeightedGraph() {
    delete &arcs;
}

void UndirectedWeightedGraph::add(int node) {
    assert(node >= 0);
    boost::add_edge(node + 1, 0, 1, arcs);
}

void UndirectedWeightedGraph::add(int from, int to) {
    assert(from != to);
    assert(from >= 0);
    assert(to >= 0);
    std::pair<WeightedAdjacencyList::edge_descriptor, bool> edge = boost::edge(from + 1, to + 1, arcs);
    if(edge.second) {
        int weight = get(boost::edge_weight_t(), arcs, edge.first);
        boost::put(boost::edge_weight_t(), arcs, edge.first, weight + 1);
        
        edge = boost::edge(to + 1, from + 1, arcs);
        assert(edge.second);
        boost::put(boost::edge_weight_t(), arcs, edge.first, weight + 1);
    }
    else {
        boost::add_edge(from + 1, to + 1, 1, arcs);
        boost::add_edge(to + 1, from + 1, 1, arcs);
    }
}

void UndirectedWeightedGraph::remove(int from, int to) {
    assert(from != to);
    assert(from >= 0);
    assert(to >= 0);
    std::pair<WeightedAdjacencyList::edge_descriptor, bool> edge = boost::edge(from + 1, to + 1, arcs);
    assert(edge.second);
    int weight = get(boost::edge_weight_t(), arcs, edge.first) - 1;
    if(weight == 0) {
        boost::remove_edge(from + 1, to + 1, arcs);
        boost::remove_edge(to + 1, from + 1, arcs);
    }
    else {
        boost::put(boost::edge_weight_t(), arcs, edge.first, weight);
        
        edge = boost::edge(to + 1, from + 1, arcs);
        assert(edge.second);
        boost::put(boost::edge_weight_t(), arcs, edge.first, weight);
    }
}

void UndirectedWeightedGraph::sccs(vec<int>& atom2comp, vec<vec<int> >& components) {
    vector<int> sccs(boost::num_vertices(arcs)), discover_time(boost::num_vertices(arcs));
    vector<boost::default_color_type> color(boost::num_vertices(arcs));
    vector<boost::graph_traits<boost::adjacency_list<> >::vertex_descriptor> root(boost::num_vertices(arcs));
    int n = boost::strong_components(arcs, &sccs[0] , boost::root_map(&root[0]).color_map(&color[0]).discover_time_map(&discover_time[0]));

    components.clear();
    if(n == 0) return;
    components.growTo(n-1);
    
    assert(sccs[0] == 0);
    for(unsigned i = 1; i < sccs.size(); i++) {
        int scc = sccs[i] - 1;
        assert(scc >= 0);
        assert(scc < components.size());
    
        atom2comp[i-1] = scc;
        components[scc].push(i-1);
    }    
}

} // namespace aspino
#include <iostream>
#include <assert.h>
#include <lemon/list_graph.h>
#include <lemon/bfs.h>
#include <lemon/dijkstra.h>
#include <lemon/preflow.h>

using namespace lemon;
using namespace std;

// TODO [EASY] Use conveninece macros to reduce amount of tpes

// TODO [EASY] Maybe subgraph does have a reference to the original graph?
// Like a subgraph, but keep around the orignal graph not to hide e.g. edged going out
template <typename MyGraph>
struct GraphSubset {
	MyGraph &graph;
	SubGraph<MyGraph> &subset;
	GraphSubset(MyGraph& g, SubGraph<MyGraph> &s) : graph(g), subset(s) {};
};

// 1. Notations
// Undirected graphs

template <typename MyGraph>
int deg(const MyGraph &g, const typename MyGraph::Node &n) {
	return countIncEdges(g, n);
}

// How will this work with an adaptor that hides verts?
template <typename MyGraph>
int vol(const MyGraph &g) {
	int volume = 0;
	for(typename MyGraph::NodeIt n(g); n != INVALID; ++n) {
		volume += deg(g, n);
	}
	return volume;
}

template <typename MyGraph>
int vol(const GraphSubset<MyGraph> &gs) {
	int volume = 0;
	for(typename SubGraph<MyGraph>::NodeIt n(gs.subset); n != INVALID; ++n) {
		volume += deg(gs.graph, n);
	}
	return volume;
}

// E(S,T) : edges between two subsets
template <typename GT>
int n_edges_between(const GraphSubset<GT> &S_, const GraphSubset<GT> &T_) {
	assert(&S_.graph == &T_.graph);
	const SubGraph<GT>& S = S_.subset;
	const SubGraph<GT>& T = T_.subset;
	const GT &G = S_.graph;

	int n = 0;
	for(typename GT::EdgeIt e(G); e != INVALID; ++e) {
		typename GT::Node u = G.u(e);
		typename GT::Node v = G.v(e);
		bool crossing = S.status(u) && T.status(v);
		crossing |= S.status(v) && T.status(u);
		n += crossing ? 1 : 0;
	}
	return n;
}

template <typename GT>
int cut_size(const GraphSubset<GT> &S_) {
	const SubGraph<GT>& S = S_.subset;
	const GT &G = S_.graph;

	int n = 0;
	for(typename GT::EdgeIt e(G); e != INVALID; ++e) {
		typename GT::Node u = G.u(e);
		typename GT::Node v = G.v(e);
		bool crossing = S.status(u) && !S.status(v);
		crossing |= S.status(v) && !S.status(u);
		n += crossing ? 1 : 0;
	}
	return n;
}

// conductance of cut S: cut-size(S) / min(vol S, vol comp s)
template <typename GT>
double conductance(const GraphSubset<GT> &S_) {
	const SubGraph<GT>& S = S_.subset;
	const GT &G = S_.graph;

	int d = cut_size(S_);
	int volS = vol(S_);
	int volcompS = vol(G) - volS;

	return ((double)d)/min(volS, volcompS);
}


/*
// This will be very interesting when we do conductances of subsets - need to track what "cuts" are opposite to inside of there.
// COnductace of graph -- total search
template <typename GT>
double conductance(const GT &G) {
	ListGraph::NodeMap<bool> filter(G, false);
	ListGraph::EdgeMap<bool> filterE(G, true);
	SubGraph<ListGraph> g_(g, filter, filterE);
	GraphSubset<ListGraph> gs(g, g_);

	// Base case for empty graph, and also its the max possible
	double min_cond = 1;
	
	vector<GT::NodeIt> nodes;
	for(typename GT::NodeIt n(G); e != INVALID; ++e) {
		nodes.push_back(n);
	}

	// We'll be doing some binary counting here
	bool carry = false;
	int i = 0;
	while(i < nodes.size()) {
		if(carry) {
			if(!filter[nodes[i]]) {
				filter[nodes[i]] = true;
				carry = false;
			} else {
				++i;
				continue;
			}
		}

		double cond = conductance(gs);
		min_cond = min(min_cond, cond);


		// ?
		i = 0;

		// "add one"
		if(!filter[nodes[i]]) {
			filter[nodes[i]] = true;
		} else {
			filter[nodes[i]] = false;
			carry = true;
			++i;
		}
	}

	return min_cond;
}
*/

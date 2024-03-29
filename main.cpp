//  Cut-matching component authored Ludvig Janiuk 2019 as part of individual
//  project at KTH. Rest, by me.
#ifdef __clang__
#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-msc32-c"
#pragma ide diagnostic ignored "cppcoreguidelines-slicing"
#endif

#include <algorithm>
#include <array>
#include <bits/stdc++.h>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>
#include <lemon/core.h>
#include <lemon/edge_set.h>
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <memory>
#include <ostream>
#include <random>
#include <set>
#include <sstream>
#include <vector>

#include "cxxopts.hpp"

#include <limits.h>
#include <omp.h>
#include <sched.h>

using namespace lemon;
using namespace std;
using namespace std::chrono;

using G = ListGraph;
using NodeMapd = typename G::template NodeMap<double>;
using Node = typename G::Node;
using NodeIt = typename G::NodeIt;
using Snapshot = typename G::Snapshot;
using Edge = typename G::Edge;
using EdgeIt = typename G::EdgeIt;
using IncEdgeIt = typename G::IncEdgeIt;
using OutArcIt = typename G::OutArcIt;
using Paths = vector<array<Node, 2>>;
using ArcLookup = ArcLookUp<G>;
template <class T>
using EdgeMap = typename G::template EdgeMap<T>;
using EdgeMapi = EdgeMap<int>;   // LEMON uses ints internally. We might want to
                                 // look into this
using EdgeMapb = EdgeMap<bool>;  // LEMON uses ints internally. We might want to
                                 // look into this
template <class T>
using NodeMap = typename G::template NodeMap<T>;
using NodeMapi = NodeMap<int>;
using NodeMapb = NodeMap<bool>;
using NodeNeighborMap = NodeMap<vector<tuple<Node, int>>>;
using FlowAlgo = Preflow<G, EdgeMapi>;
using Matching = vector<array<Node, 2>>;
using Matchingp = unique_ptr<Matching>;
using Bisection = set<Node>;
using Bisectionp = unique_ptr<Bisection>;
using Cut = set<Node>;
using Cutp = unique_ptr<Cut>;
using CutMap = NodeMap<bool>;

#define MAX_PARALLEL_RECURSIVE_LEVEL 10
#define DEBUG 1

const int NUM_THREADS = 4;
const double MICROSECS = 1000000.0;
const double PHI_UNREACHABLE = 0.00000001;
const double PHI_ACCEPTANCE_TRIMMING = 0.166;
const auto& now = high_resolution_clock::now;

auto
duration_sec(const high_resolution_clock::time_point& start,
             high_resolution_clock::time_point& stop) -> double
{
    return duration_cast<microseconds>(stop - start).count() / MICROSECS;
}


void warn(bool condition, const string& message, double line) {

    if (!condition) {
        cout << "WARNING: pre/post condition failed at line  " << line << endl;
        cout << message << endl;
        cout << "" << endl;
    }

    }

struct InputConfiguration {
    bool load_from_file = false;
    string file_name = "";
    size_t n_nodes_to_generate{};
};

struct Configuration {
    InputConfiguration input;
    bool compare_partition = false;
    string partition_file = "";
    bool seed_randomness = false;
    int seed{};
    int max_rounds{};
    bool output_cut{};
    string output_file;
    bool show_help_and_exit = false;
    bool only_test_expander = false;
    bool decompose_with_tests = false;
    bool use_H_phi_target = false;
    double H_phi_target = 0;
    bool use_G_phi_target = false;
    double G_phi_target = 0;
    // we only break if we find a good enough cut that is also this balanced
    // (has this minside volume)
    bool use_volume_treshold = false;
    double known_phi = 0.;
    double volume_treshold_factor = 1;
    double h_factor = 0;
    double h_ratio = 0;
    int n_nodes_orig = 0;
    int e_edges_orig = 0;
};

struct Logger {
    bool silent = false;
    bool verbose = false;
    ofstream nul;  // UNopened file stream, will act like /dev/null
    Logger()  {};
    auto
    progress() -> decltype(cout)&
    {
        return silent ? nul : cout;
    };
    auto
    debug() -> decltype(cout)&
    {
        return verbose ? cout : nul;
    };
} l;

struct GraphContext {
    G g;
    vector<Node> nodes;
    long num_edges{};
    map<Node, int> orig_degree;

    void
    clear()
    {
        g.clear();
        nodes.clear();
        num_edges = 0;
        orig_degree.clear();
    }
};
using GraphContextp = unique_ptr<GraphContext>;

// I'd say implementing our own adaptor is more effort, we can just do the
// snapshot thing Actually lets just subdivide manually at the start and we dont
// even need to restore.
struct SubdividedGraphContext {
    explicit SubdividedGraphContext(GraphContext& gc)
        : origContext(gc),
          nf(sub_g),
          ef(sub_g),
          n_ref(gc.g, INVALID),
          n_cross_ref(sub_g, INVALID),
          origs(sub_g, false),
          only_splits(sub_g, nf, ef){};

    GraphContext& origContext;
    G sub_g;
    NodeMapb nf;
    EdgeMapb ef;
    NodeMap<Node> n_ref;
    NodeMap<Node> n_cross_ref;
    NodeMap<bool> origs;
    SubGraph<G> only_splits;
    vector<Node> split_vertices;
};

// TODO: What chnages will be necessary?
struct RoundReport {
    size_t index{};
    size_t capacity_required_for_full_flow{};
    double multi_h_conductance{};
    double g_conductance{};
    long volume{};
    bool relatively_balanced{};
    Cutp cut;
};

template <class G>
struct CutStats {
    using Node = typename G::Node;
    using Edge = typename G::Edge;
    using Cut = set<Node>;
    using EdgeIt = typename G::EdgeIt;
    using Bisection = set<Node>;
    size_t crossing_edges = 0;

   private:
    bool is_min_side{};
    size_t min_side = 0;
    size_t cut_volume = 0;
    size_t max_side = 0;
    size_t num_edges = 0;
    auto
    degreesum() -> long
    {
        return num_edges * 2;
    }
    auto
    noncut_volume() -> long
    {
        return degreesum() - cut_volume;
    }

   public:
    CutStats(const G& g, size_t num_vertices, const Cut& cut)
    {
        initialize(g, num_vertices, cut);
    }

    void
    initialize(const G& g, size_t num_vertices, const Cut& cut)
    {
        cut_volume = 0;
        crossing_edges = 0;
        for (EdgeIt e(g); e != INVALID; ++e) {
            ++num_edges;
            if (is_crossing(g, cut, e)) {
                crossing_edges += 1;
}
            // if (any_in_cut(g, cut, e)) cut_volume += 1;
            if (cut.count(g.u(e))) {
                cut_volume += 1;
}
            if (g.u(e) == g.v(e)) {
                continue;
}
            if (cut.count(g.v(e))) {
                cut_volume += 1;
}
        }

        assert(cut.size() <= num_vertices);
        size_t other_size = num_vertices - cut.size();
        min_side = min(cut.size(), other_size);
        max_side = max(cut.size(), other_size);
        is_min_side = cut.size() == min_side;
    }

    static auto
    is_crossing(const G& g, const Bisection& c, const Edge& e) -> bool
    {
        bool u_in = c.count(g.u(e));
        bool v_in = c.count(g.v(e));
        return u_in != v_in;
    }

    static auto
    any_in_cut(const G& g, const Bisection& c, const Edge& e) -> bool
    {
        bool u_in = c.count(g.u(e));
        bool v_in = c.count(g.v(e));
        return u_in || v_in;
    }

    auto
    minside_volume() -> long
    {
        return is_min_side ? cut_volume : noncut_volume();
    }

    auto
    maxside_volume() -> long
    {
        return is_min_side ? noncut_volume() : cut_volume;
    }

    auto
    diff() -> size_t
    {
        return max_side - min_side;
    }

    auto
    num_vertices() -> size_t
    {
        return min_side + max_side;
    }

    auto
    imbalance() -> double
    {
        return diff() * 1. / num_vertices();
    }

    auto
    expansion() -> double
    {
        return min_side == 0 ? 0 : crossing_edges * 1. / min_side;
    }

    auto
    conductance() -> double
    {
        // Q: changed to 1
        return minside_volume() == 0 ? 1
                                     : crossing_edges * 1. / minside_volume();
    }

    void
    print()
    {
        cout << "Edge crossings (E) : " << crossing_edges << endl;
        cout << "cut size: (" << min_side << " | " << max_side << ")" << endl
             << "diff: " << diff() << " (" << imbalance()
             << " of total n vertices)" << endl;
        cout << "Min side: " << min_side << endl;
        // cout << "E/min(|S|, |comp(S)|) = " << expansion() << endl;
        cout << "expansion: " << expansion() << endl;
        cout << "conductance: " << conductance() << endl;
        cout << "cut volume: " << cut_volume << endl;
        cout << "noncut volume: " << noncut_volume() << endl;
    }
};
// Reads the file filename,
// creates that graph in graph g which is assumed to be empty
// In the process fills nodes with each node created at the index of (its id in
// the file minus one) And sets each node's original_ids id to be (its id in the
// file minus one). Of course original_ids must be initialized onto the graph g
// already earlier.
static void
parse_chaco_format(const string& filename, ListGraph& g, vector<Node>& nodes)
{
    assert(nodes.empty());
    l.progress() << "Reading graph from " << filename << endl;
    ifstream file;
    file.open(filename);
    if (!file) {
        cerr << "Unable to read file " << filename << endl;
        exit(1);
    }

    string line;
    stringstream ss;
    getline(file, line);
    ss.str(line);

    unsigned long n_verts;
    unsigned long n_edges;
    ss >> n_verts >> n_edges;
    l.progress() << "Reading a graph with V " << n_verts << "E " << n_edges
                 << endl;
    g.reserveNode(n_verts);
    g.reserveNode(n_edges);

    for (size_t i = 0; i < n_verts; i++) {
        Node n = g.addNode();
        nodes.push_back(n);
    }

    for (size_t i = 0; i < n_verts; i++) {
        getline(file, line);

        Node u = nodes[i];
        istringstream iss(line);
        vector<string> tokens{istream_iterator<string>{iss},
                              istream_iterator<string>{}};
        for (string& str : tokens) {
            size_t v_name = stoi(str);
            // cout << "edge to: " << v_name << "..." ;
            assert(v_name != 0);
            Node v = nodes[v_name - 1];
            // TODO: is this right?
            //Ignore multi-edges unless self-loop
            //If self loop add multiple times
            if  (findEdge(g, u, v) == INVALID || u == v) {
                g.addEdge(u, v);
            }
            //if (u == v)
            //    g.addEdge(v, u);
        }
    }

    /*
    if (n_verts % 2 != 0) {
        l.progress() << "Odd number of vertices, adding extra one." << endl;
        Node n = g.addNode();
        g.addEdge(nodes[0], n);
        nodes.push_back(n);
    }
    */
}

void
generate_large_graph(G& g, vector<Node>& nodes, size_t n_nodes)
{
    assert(n_nodes > 0);
    nodes.reserve(n_nodes);
    for (unsigned long i = 0; i < n_nodes; i++) {
        nodes.push_back(g.addNode());
    }

    g.addEdge(nodes[0], nodes[1]);
    g.addEdge(nodes[1], nodes[2]);
    g.addEdge(nodes[2], nodes[0]);

    unsigned long lim1 = n_nodes / 3;
    unsigned long lim2 = 2 * n_nodes / 3;

    for (unsigned long i = 3; i < lim1; i++) {
        ListGraph::Node u = nodes[i];
        ListGraph::Node v = nodes[0];
        g.addEdge(u, v);
    }
    for (unsigned long i = lim1; i < lim2; i++) {
        ListGraph::Node u = nodes[i];
        ListGraph::Node v = nodes[1];
        g.addEdge(u, v);
    }
    for (unsigned long i = lim2; i < n_nodes; i++) {
        ListGraph::Node u = nodes[i];
        ListGraph::Node v = nodes[2];
        g.addEdge(u, v);
    }
}

void
write_cut(const vector<Node>& nodes, const Cut& cut, const string& file_name)
{
    ofstream file;
    file.open(file_name);
    if (!file) {
        cout << "Cannot open file " << file_name << endl;
        return;
    }

    cout << "Writing partition with " << nodes.size() << " nodes to file "
         << file_name << endl;
    for (const auto& n : nodes) {
        file << (cut.count(n) != 0U ? "1" : "0") << "\n";
    }
    file.close();
}

auto
read_partition_file(const string& filename, const vector<Node>& nodes,
                    vector<set<Node>> partition) -> vector<set<Node>>
{
    ifstream file;
    file.open(filename);
    if (!file) {
        cerr << "Unable to read file " << filename << endl;
        exit(1);
    }

    ifstream is(filename);
    string raw_input;

    int n;
    while( is >> n ) {
        int iter = 0;
        while(getline(file, raw_input)) {
            partition.emplace_back();
            istringstream iss(raw_input);
            int m;
            while (iss >> m) {
                cout << m << endl;
                partition[iter].insert(nodes[m]);
            }
            cout << raw_input << endl;
            iter++;
        }
    }
    return partition;
}

void
initGraph(GraphContext& gc, const InputConfiguration& config)
{
    if (config.load_from_file) {
        parse_chaco_format(config.file_name, gc.g, gc.nodes);
    }
    else {
        l.debug() << "Generating graph with " << config.n_nodes_to_generate
                  << " nodes." << endl;
        generate_large_graph(gc.g, gc.nodes, config.n_nodes_to_generate);
    }

    for (NodeIt n(gc.g); n != INVALID; ++n) {
        for (IncEdgeIt e(gc.g, n); e != INVALID; ++e) {
            gc.orig_degree[n] += 1;
}
    }
    gc.num_edges = countEdges(gc.g);
}

// For some reason lemon returns arbitrary values for flow, the difference is
// correct tho
inline auto
flow(const ArcLookUp<G>& alp, const unique_ptr<Preflow<G, EdgeMapi>>& f, Node u,
     Node v) -> int
{
    return f->flow(alp(u, v)) - f->flow(alp(v, u));
}

void
print_end_round_message(int i)
{
    l.debug() << "======================" << endl;
    l.progress() << "== End round " << i << " ==" << endl;
    l.debug() << "======================" << endl;
}

template <typename GG>
void
print_matching(GG& g, const Matchingp& m, decltype(cout)& stream)
{
    for (auto& e : *m) {
        stream << "(" << g.id(e[0]) << ", " << g.id(e[1]) << "), ";
    }
    stream << endl;
}

void
print_cut(const Bisection& out_cut, decltype(cout)& stream)
{
    for (Node n : out_cut) {
        stream << G::id(n) << ", ";
    }
    stream << endl;
}

void
print_graph(G& g, decltype(cout)& stream)
{
    stream << "Printing a graph" << endl;
    stream << "Vertices: " << countNodes(g) << ", Edges: " << countEdges(g)
           << endl;
    stream << "==" << endl;
    for (NodeIt n(g); n != INVALID; ++n) {
        stream << G::id(n) << ", ";
    }
    stream << "\n==" << endl;
    for (EdgeIt e(g); e != INVALID; ++e) {
        stream << G::id(e) << ": " << G::id(g.u(e)) << " - " << G::id(g.v(e))
               << "\n";
    }
    stream << endl;
}

// Actually copies the graph.
void
createSubdividedGraph(SubdividedGraphContext& sgc)
{
    graphCopy(sgc.origContext.g, sgc.sub_g).nodeCrossRef(sgc.n_cross_ref).run();
    graphCopy(sgc.origContext.g, sgc.sub_g)
        .nodeRef(sgc.n_ref)
        .nodeCrossRef(sgc.n_cross_ref)
        .run();
    G& g = sgc.sub_g;
    for (NodeIt n(g); n != INVALID; ++n) {
        sgc.origs[n] = true;
    }
    vector<Edge> edges;
    for (EdgeIt e(g); e != INVALID; ++e) {
        if (g.u(e) != g.v(e)) {
            edges.push_back(e);
}
    }

    for (NodeIt n(g); n != INVALID; ++n) {
        sgc.only_splits.disable(n);
    }

    for (auto& e : edges) {
        if (g.u(e) == g.v(e)) {
            continue;
}
        Node u = g.u(e);
        Node v = g.v(e);
        g.erase(e);

        Node s = g.addNode();
        sgc.origs[s] = false;
        sgc.only_splits.enable(s);
        g.addEdge(u, s);
        if (u != v) {
            g.addEdge(s, v);
}

        // if (u !=v)
        sgc.split_vertices.push_back(s);
    }
}

struct CutMatching {
    const Configuration& config;
    GraphContext& gc;
    SubdividedGraphContext sgc;
    default_random_engine& random_engine;
    // vector<unique_ptr<RoundReport>> sub_past_rounds;
    vector<unique_ptr<RoundReport>> sub_past_rounds;
    vector<Matchingp> matchings;
    vector<Matchingp> sub_matchings;
    bool reached_H_target = false;
    // Input graph
    CutMatching(GraphContext& gc, const Configuration& config_,
                default_random_engine& random_engine_)
        : config(config_), gc(gc), sgc{gc}, random_engine(random_engine_)
    {
        assert(static_cast<int>(gc.nodes.size()) % 2 == 0);
        assert(!gc.nodes.empty()) ;
        assert(connected(gc.g));

        createSubdividedGraph(sgc);
    };

    // During the matching step a lot of local setup is actually made, so it
    // makes sense to group it inside a "matching context" that exists for the
    // duration of the mathing step
    struct MatchingContext {
        // This NEEDS to be the whole graph
        G& g;
        Node s;
        Node t;
        EdgeMapi capacity;
        CutMap cut_map;
        Snapshot snap;  // RAII

        explicit MatchingContext(G& g_)
            : g(g_), capacity(g_), cut_map(g_), snap(g_)
        {
        }

        ~MatchingContext() { snap.restore(); }

        auto
        touches_source_or_sink(Edge& e) -> bool
        {
            return g.u(e) == s || g.v(e) == s || g.u(e) == t || g.v(e) == t;
        }

        // Fills given cut pointer with a copy of the cut map
        auto
        extract_cut() -> Cutp
        {
            Cutp cut(new Cut);
            for (NodeIt n(g); n != INVALID; ++n) {
                if (n == s || n == t) {
                    continue;
}
                if (cut_map[n]) {
                    cut->insert(n);
}
            }
            return cut;
        }

        void
        reset_cut_map()
        {
            for (NodeIt n(g); n != INVALID; ++n) {
                cut_map[n] = false;
            }
        }
    };

    struct MatchResult {
        Cutp cut_from_flow;
        size_t capacity;  // First capacity (minumum) that worked to get full
                          // flow thru
    };

    static inline void
    extract_path_fast(const G&  /*g*/, const unique_ptr<Preflow<G, EdgeMapi>>&  /*f*/,
                      NodeNeighborMap& flow_children, Node u_orig,
                      Node t,  // For assertsions
                      array<Node, 2>& out_path)
    {
        out_path[0] = u_orig;
        Node u = u_orig;
        int i = 0;
        while (true) {
            i++;
            auto& vv = flow_children[u];
            assert(!vv.empty());
            auto& tup = vv.back();
            Node v = get<0>(tup);
            --get<1>(tup);

            if (get<1>(tup) == 0) {
                flow_children[u].pop_back();
}

            if (flow_children[v].empty()) {
                assert(v == t);
                assert(u != u_orig);

                out_path[1] = u;
                break;
            }

            u = v;
        }
    }

    void
    decompose_paths_fast(const MatchingContext& mg,
                         const unique_ptr<FlowAlgo>& f, Paths& out_paths)
    {
        f->startSecondPhase();
        EdgeMapi subtr(mg.g, 0);
        NodeNeighborMap flow_children(mg.g, vector<tuple<Node, int>>());
        out_paths.reserve(countNodes(mg.g) / 2);

        // Calc flow children (one pass)
        ArcLookup alp(mg.g);
        for (EdgeIt e(mg.g); e != INVALID; ++e) {
            if (mg.g.u(e) == mg.g.v(e)) {
                continue;
}
            Node u = mg.g.u(e);
            Node v = mg.g.v(e);
            long e_flow = flow(alp, f, u, v);
            if (e_flow > 0) {
                flow_children[u].push_back(tuple(v, e_flow));
            }
            else if (e_flow < 0) {
                flow_children[v].push_back(tuple(u, -e_flow));
            }
        }

        for (IncEdgeIt e(mg.g, mg.s); e != INVALID; ++e) {
            if (mg.g.u(e) == mg.g.v(e)) {
                continue;
}
            assert(mg.g.u(e) == mg.s || mg.g.v(e) == mg.s);
            Node u = mg.g.u(e) == mg.s ? mg.g.v(e) : mg.g.u(e);

            out_paths.push_back(array<Node, 2>());
            extract_path_fast(mg.g, f, flow_children, u, mg.t,
                              out_paths[out_paths.size() - 1]);
        }
    }

    // Works for sub too, with the assumption that mg.g realy is the whole graph
    static void
    run_min_cut(const MatchingContext& mg, unique_ptr<FlowAlgo>& p) 
    {
        p = std::make_unique<Preflow<G, EdgeMapi>>(mg.g, mg.capacity, mg.s, mg.t);
        auto start2 = now();
        p->runMinCut();  // Note that "startSecondPhase" must be run to get
                         // flows for individual verts
        auto stop2 = now();
        l.progress() << "flow: " << p->flowValue() << " ("
                     << duration_sec(start2, stop2) << " s)" << endl;
    }

    // This should work well for sub too
    static void
    set_matching_capacities(MatchingContext& mg, size_t cap) 
    {
        for (EdgeIt e(mg.g); e != INVALID; ++e) {
            if (mg.touches_source_or_sink(e)) {
                continue;
}
            mg.capacity[e] = cap;
        }
    }

    auto
    bin_search_flows(MatchingContext& mg, unique_ptr<FlowAlgo>& p,
                     unsigned long flow_target) const -> MatchResult
    {
        auto start = now();
        unsigned long cap = 1;
        // for (; cap < mg.static_cast<int>(gc.nodes.size()); cap *= 2) {
        for (; cap < flow_target * 2; cap *= 2) {
            l.progress() << "Cap " << cap << " ... " << flush;
            set_matching_capacities(mg, cap);
            run_min_cut(mg, p);

            // bool reachedFullFlow = p->flowValue() == mg.static_cast<int>(gc.nodes.size()) / 2;
            bool reachedFullFlow = static_cast<unsigned long>(p->flowValue()) >= flow_target;
            if (reachedFullFlow) {
                l.debug()
                    << "We have achieved full flow, but half this capacity "
                       "didn't manage that!"
                    << endl;
}

            // So it will always have the mincutmap of "before"
            // mincuptmap is recomputed too many times of course but whatever
            // If we reached it with cap 1, already an expander I guess?
            // In this case this was never done even once, so we have to do it
            // before breaking
            if (!reachedFullFlow || cap == 1) {
                mg.reset_cut_map();
                p->minCutMap(mg.cut_map);
            }

            if (reachedFullFlow) {
                break;
}
        }

        // Not we copy out the cut
        MatchResult result{mg.extract_cut(), cap};

        auto stop = now();
        l.progress() << "Flow search took (seconds) "
                     << duration_sec(start, stop) << endl;

        return result;
    }

    void
    decompose_paths(const MatchingContext& mg, const unique_ptr<FlowAlgo>& p,
                    vector<array<Node, 2>>& paths)
    {
        decompose_paths_fast(mg, p, paths);
    }

    template <typename GG>
    void
    make_sink_source(GG& g, MatchingContext& mg, const set<Node>& cut) const
    {
        mg.s = g.addNode();
        mg.t = g.addNode();
        int s_added = 0;
        int t_added = 0;
        for (typename GG::NodeIt n(g); n != INVALID; ++n) {
            if (n == mg.s) {
                continue;
}
            if (n == mg.t) {
                continue;
}
            Edge e;
            if (cut.count(n)) {
                e = g.addEdge(mg.s, n);
                s_added++;
            }
            else {
                e = g.addEdge(n, mg.t);
                t_added++;
            }
            mg.capacity[e] = 1;
        }
        int diff = s_added - t_added;
        assert(-1 <= diff && diff <= 1);
    }

    // Actually, cut player gets H
    // Actually Actually, sure it gets H but it just needs the matchings...
    // TODO: Ok so can we just call this with split_only and matchings of those?
    template <typename GG, typename M>
    auto
    cut_player(const GG& g, const vector<unique_ptr<M>>& given_matchings,
               double& h_multi_cond_out) -> Bisectionp
    {
        l.debug() << "Running Cut player" << endl;
        typename GG::template NodeMap<double> probs(g);
        vector<Node> all_nodes;

        // uniform_int_distribution<int> uniform_dist(0, 1);
        for (typename GG::NodeIt n(g); n != INVALID; ++n) {
            all_nodes.push_back(n);
        }
        uniform_int_distribution<int> uniform_dist(0, 1);
        for (typename GG::NodeIt n(g); n != INVALID; ++n) {
            probs[n] = uniform_dist(random_engine) ? 1.0 / all_nodes.size()
                                                   : -1.0 / all_nodes.size();
        }

        size_t num_vertices = all_nodes.size();

        ListEdgeSet H(g);
        ListEdgeSet H_single(g);
        for (const unique_ptr<M>& m : given_matchings) {
            for (auto& e : *m) {
                Node u = e[0];
                Node v = e[1];
                double avg = probs[u] / 2 + probs[v] / 2;
                probs[u] = avg;
                probs[v] = avg;

                H.addEdge(u, v);
                // Updating H_single
                if (findEdge(H_single, u, v) == INVALID) {
                    assert(findEdge(H_single, v, u) == INVALID);
                    H_single.addEdge(u, v);
                }
            }
        }
        shuffle(all_nodes.begin(), all_nodes.end(), random_engine);
        sort(all_nodes.begin(), all_nodes.end(),
             [&](Node a, Node b) { return probs[a] < probs[b]; });

        size_t size = all_nodes.size();
        // With subdivisions, won't be this way longer
        // assert(size % 2 == 0);
        all_nodes.resize(size / 2);
        auto b = std::make_unique<Bisection>(all_nodes.begin(), all_nodes.end());
        l.debug() << "Cut player gave the following cut: " << endl;
        print_cut(*b, l.debug());

        // So how does it give output?
        // Ok it assigns h_outs, but actually also returns Bisectionp
        auto hcs = CutStats<decltype(H)>(H, num_vertices, *b);
        l.progress() << "H conductance: " << hcs.conductance()
                     << ", num cross: " << hcs.crossing_edges << endl;
        h_multi_cond_out = hcs.conductance();
        auto hscs = CutStats<decltype(H_single)>(H_single, num_vertices, *b);
        l.progress() << "H_single conductance: " << hscs.conductance()
                     << ", num cross: " << hscs.crossing_edges << endl;

        return b;
    }

    // returns capacity that was required
    // Maybe: make the binsearch an actual binsearch
    // TODO: Let listedgeset just be 2-arrays of nodes. Lemon is getting in the
    // way too much. But also what is assigned in MtchResult?
    auto
    matching_player(const set<Node>& bisection, Matching& m_out) -> MatchResult
    {
        MatchingContext mg(gc.g);
        make_sink_source(mg.g, mg, bisection);

        unique_ptr<FlowAlgo> p;
        MatchResult mr = bin_search_flows(mg, p, static_cast<int>(gc.nodes.size()) / 2);

        decompose_paths(mg, p, m_out);

        // Now how do we extract the cut?
        // In this version, in one run of the matching the cut is strictly
        // decided. We just need to decide which one of them. Only when we
        // change to edge will the cut need to be explicitly extracted. Rn the
        // important thing is to save cuts between rounds so I can choose the
        // best.

        return mr;
    }

    // returns capacity that was required
    // Maybe: make the binsearch an actual binsearch
    // TODO: Let listedgeset just be 2-arrays of nodes. Lemon is getting in the
    // way too much. But also what is assigned in MtchResult?
    auto
    sub_matching_player(const set<Node>& bisection,
                        vector<array<Node, 2>>& m_out) -> MatchResult
    {
        MatchingContext mg(sgc.sub_g);
        make_sink_source(sgc.only_splits, mg, bisection);
        // TODO: so have s, t been created on the big graph?
        // cout << id << endl;
        // cout << countNodes(sgc.sub_g) << endl;

        unique_ptr<FlowAlgo> p;
        MatchResult mr = bin_search_flows(mg, p, sgc.split_vertices.size() / 2);

        decompose_paths(mg, p, m_out);

        // Now how do we extract the cut?
        // In this version, in one run of the matching the cut is strictly
        // decided. We just need to decide which one of them. Only when we
        // change to edge will the cut need to be explicitly extracted. Rn the
        // important thing is to save cuts between rounds so I can choose the
        // best.

        return mr;
    }

    auto
    volume_treshold() -> long
    {
        return config.volume_treshold_factor * gc.num_edges;
    }
    // IS this right?
    auto
    sub_volume_treshold() -> long
    {
        // return config.volume_treshold_factor * sgc.origContext.nodes.size();
        if (config.h_ratio > 0) {
            return sgc.origContext.num_edges * config.h_ratio;
        }  {
            return (double(sgc.origContext.num_edges) /
                    (10. * config.volume_treshold_factor *
                     pow(log(sgc.origContext.num_edges), 2)));
}
    }

    /*
    // Ok lets attack from here
    // Theres a lot of risk for problems with "is this a cut in the orig graph
    or in the splits? unique_ptr<RoundReport> one_round() {
    unique_ptr<RoundReport> report = make_unique<RoundReport>();

        // So this needs to be a bisection of splitnodes, I feel this could be
    very opaque.
        // We'd need a subgraph that is just the splitnodes
        Bisectionp bisection = cut_player(gc.g, matchings,
    report->multi_h_expansion);

        // Matching on splitnodes, but we need to also save the actual cut...
        Matchingp matchingp(new Matching());
        MatchResult mr = matching_player(*bisection, *matchingp);
        matchings.push_back(move(matchingp));

        //c.cuts.push_back(move(mr.cut_from_flow));
        report->index = sub_past_rounds.size();
        report->capacity_required_for_full_flow = mr.capacity;
        report->cut = move(mr.cut_from_flow);
        auto cs = CutStats<G>(gc.g, static_cast<int>(gc.nodes.size()), *report->cut);
        report->g_expansion = cs.expansion();
        l.progress() << "G cut expansion " << report->g_expansion << endl;
        report->volume = cs.minside_volume();
        l.progress() << "G cut minside volume " << cs.minside_volume() << endl;
        l.progress() << "G cut maxside volume " << cs.maxside_volume() << endl;
        report->relatively_balanced = report->volume > volume_treshold();
        return move(report);
    }
  */
    // Ok lets attack from here
    // Theres a lot of risk for problems with "is this a cut in the orig graph
    // or in the splits?
    auto
    sub_one_round() -> unique_ptr<RoundReport>
    {
        unique_ptr<RoundReport> report = make_unique<RoundReport>();

        // Bisectionp sub_bisection = cut_player(sgc.only_splits, sub_matchings,
        // h_multi_out_sub);
        Bisectionp sub_bisection = cut_player(sgc.only_splits, sub_matchings,
                                              report->multi_h_conductance);
        // Well ok, it's doing the first random thing well.
        // TODO: test on rest...

        Matchingp smatchingp(new Matching());
        MatchResult smr = sub_matching_player(*sub_bisection, *smatchingp);
        sub_matchings.push_back(move(smatchingp));

        // TODO: ocmputation of cut at the end is_ wrong...
        // c.cuts.push_back(move(mr.cut_from_flow));
        report->index = sub_past_rounds.size();
        report->capacity_required_for_full_flow = smr.capacity;
        report->cut = make_unique<Cut>();
        for (auto& n : *(smr.cut_from_flow)) {
            // if(sgc.n_cross_ref[n] != INVALID) {
            report->cut->insert(sgc.n_cross_ref[n]);
            //}
        }
        auto cs = CutStats<G>(sgc.origContext.g, sgc.origContext.nodes.size(),
                              *report->cut);
        // report->g_expansion = cs.expansion();
        // l.progress() << "SUBG cut expansion " << report->g_expansion << endl;
        report->g_conductance = cs.conductance();
        /*
        if (report->g_conductance == 1) {
            cout << "LYING" << endl;
        }
        */
        l.progress() << "SUBG cut conductance: " << report->g_conductance
                     << endl;
        report->volume = cs.minside_volume();
        l.progress() << "SUBG cut minside volume " << cs.minside_volume()
                     << endl;
        l.progress() << "SUBG cut maxside volume " << cs.maxside_volume()
                     << endl;
        // Q: >, >=?
        report->relatively_balanced = report->volume >= sub_volume_treshold();
        //cout << "Volume is: " << report->volume << " Threshold is: " <<
        //    sub_volume_treshold() << endl;
        return report;
        // Thing is, the cut is not a "materialized cut" in the subdiv graph.
        // But the stats we want are the implied cut on the orig graph.
    }

    /*
    // Stopping condition
    bool should_stop() {
        int i = sub_past_rounds.size();
        if(i == 0) return false;
        if(i >= config.max_rounds && config.max_rounds != 0) return true;

        const auto& last_round = sub_past_rounds[sub_past_rounds.size() - 1];
        if(config.use_H_phi_target && last_round->multi_h_expansion >=
    config.H_phi_target) { cout << "H Expansion target reached, this will be
    case 1 or 3. According to theory, this means we probably won't find a better
    cut. That is, assuming you set H_phi right. " "If was used together with
    G_phi target, this also certifies the input graph is a G_phi expander unless
    there was a very unbaanced cut somewhere, which we will proceed to look
    for." << endl; reached_H_target = true; return true;
        }

        if(config.use_G_phi_target)
        if(last_round->g_expansion >= config.G_phi_target) {
            if(config.use_volume_treshold && last_round->relatively_balanced) {
                cout << "CASE2 G Expansion target reached with a cut that is
    relatively balanced. Cut-matching game has found a balanced cut as good as
    you wanted it."
                     << endl;
                return true;
            }

            if(!config.use_volume_treshold) {
                cout << "G Expansion target reached. Cut-matching game has found
    a cut as good as you wanted it. Whether it is balanced or not is up to you."
                     << endl;
                return true;
            }
        }
    }
     */

    // Stopping condition
    auto
    sub_should_stop() -> bool
    {
        int i = sub_past_rounds.size();
        if (i == 0 || i == 1) {
            return false;
}
        if (i >= config.max_rounds && config.max_rounds != 0) {
            return true;
}

        const auto& last_round = sub_past_rounds[sub_past_rounds.size() - 1];

        // if (config.use_volume_treshold && (last_round->relatively_balanced))
        // {
        //    cout << "balanced cut" << endl;
        //
        //    return true;
        //}
        if (config.use_H_phi_target &&
            last_round->multi_h_conductance >= config.H_phi_target) {
            cout
                << "H Conductance target reached, this will be case 1 or 3. "
                   "According to theory, this means we probably won't find a "
                   "better "
                   "cut. That is, assuming you set H_phi right. "
                   "If was used together with G_phi target, this also "
                   "certifies the "
                   "input graph is a G_phi expander unless there was a very "
                   "unbaanced cut somewhere, which we will proceed to look for."
                << endl;
            reached_H_target = true;
            return true;
        }

        if (!(*last_round->cut).empty() &&
            last_round->g_conductance < config.G_phi_target) {
            // return true;
            if (last_round->relatively_balanced) {
                cout << "CASE2 G Expansion target reached with a cut that is "
                        "relatively balanced. Cut-matching game has found a "
                        "balanced "
                        "cut as good as you wanted it."
                     << endl;
                return true;
            }
        }

        return false;
    }

    void
    run()
    {
        while (!sub_should_stop()) {
            // sub_past_rounds.push_back(one_round());
            sub_past_rounds.push_back(sub_one_round());
            print_end_round_message(sub_past_rounds.size() - 1);
        }
    }
};

// TODO: Make cut always the smallest (maybe)
// TODO (In edge version) Implement breaking-logik for unbalance
// om vi hittar phi-cut med volym obanför treshold
// Om vi hittar phi-cut med volum under treshold, så ingorerar vi det och kör p
// och sen om vi når H, då definieras det bästa som phi-cuttet med högsta volym
auto
create_options() -> cxxopts::Options
{
    cxxopts::Options options(
        "executable_name",
        "Individual project implementation of thatchapon's paper to find graph partitions. Currently only cut-matching game. \
                             \nRecommended usage: \n\texecutable_name -s -f ./path/to/graph -o output.txt\
                             \nCurrently only running a fixed amount of rounds is supported, but the more correct \
                             \nversion of running until H becomes an expander is coming soon.\
                             ");
    options.add_options()("h,help", "Show help")(
        "H_phi",
        "Phi expansion treshold for the H graph. Recommend to also set -r=0. ",
        cxxopts::value<double>()->implicit_value("10.0"))(
        "G_phi",
        "Phi expansion target for the G graph. Means \"what is a good enough "
        "cut?\" Recommended with -r=0. This is the PHI from the paper. ",
        cxxopts::value<double>()->implicit_value("0.8"))(
        "vol",
        "Volume treshold. Only used if G_phi is used. Will be multiplied by "
        "number of edges, so to require e.g. minimum 10% volume, write 0.1.",
        cxxopts::value<double>()->implicit_value("1"))(
        "known_phi",
        "use explicit value for g phi in some tests",
        cxxopts::value<double>()->implicit_value("0.0"))(
        "only_test_expander",
        "test trimming (don't decompose)",
        cxxopts::value<bool>()->implicit_value("false"))(
        "decompose_with_tests",
        "perform more correct assertions (slower)",
        cxxopts::value<bool>()->implicit_value("false"))(
        "f,file", "File to read graph from", cxxopts::value<std::string>())(
        "r,max-rounds", "Number of rounds after which to stop (0 for no limit)",
        cxxopts::value<long>()->default_value("25"))(
        "s,seed", "Use a seed for RNG (optionally set seed manually)",
        cxxopts::value<int>()->implicit_value("1337"))(
        "p,partition", "Partition file to compare with",
        cxxopts::value<std::string>())(
        "o,output",
        "Output computed cut into file. The cut is written as the vertices of "
        "one side of the cut.",
        cxxopts::value<std::string>())(
        "n,nodes",
        "Number of nodes in graph to generate. Should be even. Ignored if -f "
        "is "
        "set.",
        cxxopts::value<long>()->default_value("100"))(
        "h_factor", "factor of push relabel height",
        cxxopts::value<double>()->implicit_value("1"))(
        "h_ratio", "use balance ratio instead",
        cxxopts::value<double>()->implicit_value("0.0"))(
        "S,Silent", "Only output one line of summary at the end")(
        "v,verbose",
        "Debug; Whether to print nodes and cuts Does not include "
        "paths. Produces a LOT of output on large graphs.")(
        "d,paths", "Debug; Whether to print paths");
    return options;
}

void
parse_options(int argc, char** argv, Configuration& config)
{
    auto cmd_options = create_options();
    auto result = cmd_options.parse(argc, argv);

    if (result.count("help") != 0U) {
        config.show_help_and_exit = true;
        cout << cmd_options.help() << endl;
    }
    if (result.count("H_phi") != 0U) {
        config.use_H_phi_target = true;
        config.H_phi_target = result["H_phi"].as<double>();
    }
    if (result.count("G_phi") != 0U) {
        config.use_G_phi_target = true;
        config.G_phi_target = result["G_phi"].as<double>();
    }
    if (result.count("vol") != 0U) {
        config.use_volume_treshold = true;
        config.volume_treshold_factor = result["vol"].as<double>();
    }
    if (result.count("file") != 0U) {
        config.input.load_from_file = true;
        config.input.file_name = result["file"].as<string>();
    }
    if (result.count("nodes") != 0U) {
        assert(!config.input.load_from_file);
}
    config.input.n_nodes_to_generate = result["nodes"].as<long>();
    if (result.count("max-rounds") != 0U) {
        config.max_rounds = result["max-rounds"].as<long>();
}
    if (result.count("verbose") != 0U) {
        l.verbose = result["verbose"].as<bool>();
}
    if (result.count("Silent") != 0U) {
        l.silent = result["Silent"].as<bool>();
}

    if (result.count("seed") != 0U) {
        config.seed_randomness = true;
        config.seed = result["seed"].as<int>();
    }

    if (result.count("output") != 0U) {
        config.output_cut = true;
        config.output_file = result["output"].as<string>();
    }
    if (result.count("partition") != 0U) {
        config.compare_partition = true;
        config.partition_file = result["partition"].as<string>();
    }
    if (result.count("h_factor") != 0U) {
        config.h_factor = result["h_factor"].as<double>();
    }
    if (result.count("h_ratio") != 0U) {
        config.h_ratio = result["h_ratio"].as<double>();
    }

    if (result.count("only_test_expander") != 0U) {
        config.only_test_expander = true; //result["only_test_expander"].as<bool>();
    }
    if (result.count("decompose_with_tests") != 0U) {
        config.decompose_with_tests = true; //result["decompose_with_tests"].as<bool>();
    }
    if (result.count("known_phi") != 0U) {
        config.known_phi = result["known_phi"].as<double>();
    }


}

using DG = ListDigraph;
using DGNode = typename DG::Node;
template <class T>
using DGNodeMap = typename DG::template NodeMap<T>;
using DGNodeIt = typename DG::NodeIt;
using Arc = typename DG::Arc;
template <class T>
using ArcMap = typename DG::template ArcMap<T>;
using DGOutArcIt = typename DG::OutArcIt;



auto
cut_complement(vector<Node> nodes, const set<Node>& cut) -> set<Node>
{
    set<Node> R;

    std::set_difference(nodes.begin(), nodes.end(), cut.begin(), cut.end(),
                        std::inserter(R, R.end()));

    return R;
}

auto
digraph_from_graph(G& g, DG& dg) -> ListDigraph*
{
    for (NodeIt n(g); n != INVALID; ++n) {
        dg.addNode();
}

    for (EdgeIt a(g); a != INVALID; ++a) {
        dg.addArc(DG::nodeFromId(G::id(g.u(a))),
                  DG::nodeFromId(G::id(g.v(a))));
        dg.addArc(DG::nodeFromId(G::id(g.v(a))),
                  DG::nodeFromId(G::id(g.u(a))));
    }

    return &dg;
}

auto
digraph_from_graph(G& g, DG& dg, ArcMap<Arc>& reverse_arc) -> ListDigraph*
{
    for (NodeIt n(g); n != INVALID; ++n) {
        dg.addNode();
}

    for (EdgeIt a(g); a != INVALID; ++a) {
        Arc a1 = dg.addArc(DG::nodeFromId(G::id(g.u(a))),
                           DG::nodeFromId(G::id(g.v(a))));
        Arc a2 = dg.addArc(DG::nodeFromId(G::id(g.v(a))),
                           DG::nodeFromId(G::id(g.u(a))));
        reverse_arc[a1] = a2;
        reverse_arc[a2] = a1;
    }

    return &dg;
}



struct node_label_pair {
    DGNode node;
    int    label{};

    auto operator<(const node_label_pair& other) const -> bool {
        return label > other.label;
    }
};

struct flow_instance {
    flow_instance(DG& dg_, double phi)
        : edge_flow(dg_, 0.0),
          edge_cap(dg_, 2.0 / phi),
          reverse_arc(dg_),
          node_flow(dg_, 0.0),
          node_cap(dg_, 0.0),
          initial_mass(dg_, 0.0),
          node_label(dg_, 0){};

    unsigned long h{};
    double phi{};
    int n{};
    int e{};

    ArcMap<double> edge_flow;
    ArcMap<double> edge_cap;
    ArcMap<Arc> reverse_arc;

    DGNodeMap<double> node_flow;
    DGNodeMap<double> node_cap;
    DGNodeMap<double> initial_mass;
    DGNodeMap<unsigned long> node_label;

    set<DGNode> A;
    set<DGNode> R;

    priority_queue<node_label_pair> active_nodes;
    set<DGNode> trimmed_last_round;
};

auto
OutDegree(DG& dg, DGNode n) -> int
{
    int s = 0;
    for (DGOutArcIt e(dg, n); e != INVALID; ++e) {
        s += 1;
    }
    return s;
}

// Q: should iterate over edges in cut
auto
level_cut(DG& dg, flow_instance& fp, vector<list<DGNode>>& level_queue) -> set<DGNode>
{
    int cut_index = level_queue.size() - 1;
    int s_vol = 0;
    int z_i = 0;

    set<DGNode> trimmed_nodes;
    int e = fp.e;
    unsigned long n_upper_nodes = 0;
    cout << "local flow: start level cut" << endl;

    for (int i = level_queue.size() - 1; i > 0; i--) {
        cut_index = i;
        // Q: ? Assert false?
        // if (level_queue[i].size() == 0)
        //    continue;
        s_vol = 0;
        z_i = 0;
        if (!level_queue[i].empty()) {
            cout << "local flow: level queue at level: " << i
                 << " has size: " << level_queue[i].size() << endl;
}
        // else
        //    continue;
        n_upper_nodes += level_queue[i].size();
        for (auto& n : level_queue[i]) {
            if (fp.R.count(n) != 0U) {
                assert(false);
                continue;
            }
            //assert(fp.A.count(n));
            for (DGOutArcIt e(dg, n); e != INVALID; ++e) {
                // Q: is this right?
                // if (!fp.A.count(dg.source(e)) || !fp.A.count(dg.target(e)))
                //    continue;
                s_vol++;
                if (fp.node_label[n] - 1 == fp.node_label[dg.target(e)]) {
                    ;
}
                {
                    z_i++;
                }
            }
        }
        if (z_i <= 5. * s_vol * log(e) /
                       level_queue.size()) {  // Q: fp.h in denominator?
            cut_index = i;
            break;
        }
    }
    cout << "local flow: s_vol: " << s_vol << "log(e)" << log(e) << "e: " << e
         << "fp.h:  " << fp.h << endl;
    cout << "local flow: cut thresh is: " << 5.0 * s_vol * log(e) / fp.h
         << endl;
    cout << "local flow: cut index is: " << cut_index
         << " A size: " << fp.A.size() << endl;
    //q: why?
    warn(cut_index > 0, "level cut should be larger than 0", __LINE__);
    trimmed_nodes.clear();

    cout << "local flow: n upper nodes: " << n_upper_nodes << endl;
    assert(n_upper_nodes <= fp.A.size());
    if (n_upper_nodes == 0) {
        return trimmed_nodes;
}

    // if (n_upper_nodes < fp.A.size()/2) {
    if (true) {
        cout << "local flow: cut on upper side" << endl;
        for (int i = level_queue.size() - 1; i >= cut_index; i--) {
            for (auto& n : level_queue[i]) {
                // Q: remove all these as they will be very slow. Or add debug.
                // assert(fp.R.count(n) == 0 && trimmed_nodes.count(n) == 0);
                fp.R.insert(n);
                trimmed_nodes.insert(n);
            }

            level_queue[i].clear();
        }
    }

    unsigned long old_A = fp.A.size();
    fp.A.clear();

    for (auto& lst : level_queue) {
        for (auto& n : lst) {
            assert(fp.R.count(n) == 0 && trimmed_nodes.count(n) == 0);
            fp.A.insert(n);
        }
    }
    cout << "local flow: trimmed nodes size: " << trimmed_nodes.size() << endl;
    cout << "local flow: retained nodes size: " << fp.A.size() << endl;
    assert(trimmed_nodes.size() + fp.A.size() == old_A &&
           "trimmed nodes and new A size is old A size");

    return trimmed_nodes;
}

void
adjust_flow_source(DG& dg, flow_instance& fp, const set<DGNode>& trimmed_nodes)
{
    assert(trimmed_nodes.size() <= fp.R.size());
    for (auto& n : trimmed_nodes) {
        //assert(!fp.A.count(n));
        for (DGOutArcIt e(dg, n); e != INVALID; ++e) {
            //assert(!fp.A.count(n));
            if (fp.A.count(dg.target(e)) != 0U) {
                fp.node_flow[dg.target(e)] += 2.0/fp.phi - fp.edge_flow[e];
                fp.initial_mass[dg.target(e)] += 2.0 / fp.phi; // - fp.edge_flow[e];
                // Q: should not be relevant
                // fp.edge_flow[fp.reverse_arc[e]] -= 2/fp.phi;
            }
        }
    }
    }

void
uf(DG& dg, flow_instance& fp, const Configuration& conf)
{
    bool excess_flow;
    (void)excess_flow;
    //assert(fp.A.size() > 0 && fp.R.size() > 0 && fp.A.size() >= fp.R.size());
    assert(fp.A.size() >= fp.R.size());
    // Set up level queues
    vector<list<DGNode>> level_queue;
    for (unsigned long i = 0; i < fp.h; i++) {
        level_queue.emplace_back();
    }

    level_queue.clear();

    //Q: fp.h
    for (unsigned long i = 0; i < fp.h && i < fp.A.size(); i++) {
        level_queue.emplace_back();
    }

    for (unsigned long i = 0; i < fp.h && i < fp.A.size(); i++) {
        level_queue[i].clear();
        assert(level_queue[i].empty());
    }

    for (auto& n : fp.A) {
        level_queue[0].push_back(n);
        assert(fp.node_flow[n] >= 0);
    }

    cout << "local flow: start unit flow" << endl;

    if (fp.h < level_queue.size()) {
        fp.h = level_queue.size();
}
    
    double pushed_flow = 0.0;
    int flow_iter = 0;
unit_start:
    assert(level_queue.size() <= fp.h);
    int cut_size_before = fp.A.size();
    while (!fp.active_nodes.empty()) {
        //int lowest_label = fp.h;
        DGNode n;
        //Q: make priority queue
        node_label_pair nlp = fp.active_nodes.top();
        n = nlp.node;
        //cout << "local flow: working on: " << fp.node_flow[n] << endl;

        excess_flow = true;
        assert(fp.node_flow[n] > 0);
        assert(fp.node_label[n] < fp.h - 1);
        assert(fp.node_flow[n] > fp.node_cap[n]);  // push-relabel
        assert(fp.h >= level_queue.size());

        for (DGOutArcIt e(dg, n); e != INVALID; ++e) {
            if (dg.source(e) == dg.target(e)) {
                continue;
}
            if ((fp.A.count(dg.source(e)) == 0U) || (fp.A.count(dg.target(e)) == 0U)) {
                // cout << "local flow: try to push inside R" << endl;
                continue;
            }
            assert(n == dg.source(e));

            if (fp.node_label[n] > fp.node_label[dg.target(e)] + 1) {
                assert(fp.edge_flow[e] >= 2. / fp.phi);
}

            if (fp.node_label[n] - fp.node_label[(dg.target(e))] == 1 &&
                fp.edge_cap[e] - fp.edge_flow[e] > 0) { //&& 
                //fp.node_flow[dg.target(e)] - fp.node_cap[dg.target(e)] < 0) {  // push
                //Source excess flow, sink has remaining capacity
                //Q: Why not?
                //assert(fp.node_flow[dg.target(e)] - fp.node_cap[dg.target(e)] <
                //       0);
                double ex_flow = fp.node_flow[n] - fp.node_cap[n];
                double res_flow = 2. / fp.phi - fp.edge_flow[e];
                double deg_min_ex = fp.node_cap[n] * fp.phi / 2.0 -
                                    (fp.node_cap[n] - fp.node_flow[n]);
                double phi = min(ex_flow, min(res_flow, deg_min_ex));
                assert(phi > 0);
                fp.node_flow[n] -= phi;
                fp.node_flow[dg.target(e)] += phi;
                assert(n != dg.target(e));
                fp.edge_flow[e] += phi;
                fp.edge_flow[fp.reverse_arc[e]] -= phi;
                //Q: is inactive or active and was not before
                if (fp.node_flow[n] <=
                    fp.node_cap[n] &&
                    fp.node_flow[n] + phi >
                    fp.node_cap[n]) {
                    //if (conf.decompose_with_tests)
                    //    assert(find(fp.active_nodes.begin(), fp.active_nodes.end(),
                    //                n) != fp.active_nodes.end());
                    //fp.active_nodes.erase(n);   
                    //Q: implicit delete
                    fp.active_nodes.pop();
                }
                if (fp.node_flow[dg.target(e)] >
                    fp.node_cap[dg.target(e)] &&
                    fp.node_flow[dg.target(e)] - phi <=
                    fp.node_cap[dg.target(e)] && 
                    fp.node_label[dg.target(e)] < fp.h -1) {
                    //if (conf.decompose_with_tests)
                    //    assert(find(fp.active_nodes.begin(), fp.active_nodes.end(),
                    //                dg.target(e)) == fp.active_nodes.end());
                    node_label_pair new_active_pair;
                    new_active_pair.node = dg.target(e);
                    new_active_pair.label = fp.node_label[dg.target(e)];
                    fp.active_nodes.push(new_active_pair);
                }
                pushed_flow += phi;
                //fp.active_nodes.push(nlp);
                //if (flow_iter % 1000 == 0)  {
                //    cout << "local flow: trimming progress, pushed flow: " << pushed_flow << endl;
                //}
                goto unit_start;
            }
        }

        // relabel

        if (conf.decompose_with_tests) {
            for (DGOutArcIt e(dg, n); e != INVALID; ++e) {
                //Q: costly check, only debug
                if (n == dg.target(e) || (fp.A.count(dg.target(e)) == 0U)) {
                    continue;
}
                assert(n == dg.source(e));
                if (fp.edge_cap[e] - fp.edge_flow[e] > 0) {
                    if (fp.node_label[n] > fp.node_label[dg.target(e)]) {
                        cout << "local flow: first label not smaller" << DG::id(n)
                            << " : " << DG::id(dg.target(e)) << endl;
                        cout << "local flow: " << fp.node_label[n] << " : "
                            << fp.node_label[dg.target(e)] << endl;
                        assert(false);
                    }
                }
                // assert (fp.node_label[n] <= fp.node_label[dg.target(e)]);
            }
        }

        int lv = fp.node_label[n];
        //assert(find(level_queue[lv].begin(), level_queue[lv].end(), n) !=
        //       level_queue[lv].end());
        unsigned long c1 = level_queue[lv + 1].size();
        level_queue[lv + 1].push_back(n);
        unsigned long c2 = level_queue[lv].size();
        level_queue[lv].remove(n);
        fp.node_label[n] = fp.node_label[n] + 1;
        nlp.label = nlp.label + 1;
        assert(level_queue[lv + 1].size() == c1 + 1);
        assert(level_queue[lv].size() == c2 - 1);
        if (fp.node_label[n] >= fp.h - 1) {
            fp.active_nodes.pop();
        }   //implicit deletion

            //fp.active_nodes.erase(n);

        //goto unit_start;
    }

    flow_iter++;

    cout << "local_flow: finish run of uf" << endl;
    cout << "local_flow: cut size before: " << cut_size_before << endl;
    cout << "local_flow: complement size before: " << fp.R.size() << endl;
    set<DGNode> trimmed_nodes = level_cut(dg, fp, level_queue);
    adjust_flow_source(dg, fp, trimmed_nodes);
    //Q: another assertion
    double vol_prune = 0.;
    int cut_edges = 0;
    //Q: going to dg, we have doubled volume, self-loops, cut edges
    //Works for known expanders, of course.
    for (auto& n : trimmed_nodes) {
        for (DGOutArcIt e(dg, n) ; e != INVALID; ++e) {
            vol_prune += 1;
            if (fp.A.count(dg.target(e)) != 0U) {
                cut_edges += 1;  
}
        }
    }
    cout << "vol prune: " << vol_prune << "  flow iter 8 / fp phi: " << flow_iter * 8 / fp.phi << endl; 

    warn(vol_prune <= flow_iter * 8 / fp.phi, "volume of cut too large", __LINE__);
    warn(cut_edges <= 4 * flow_iter, "too many cut edges", __LINE__);

    cout << "local flow: cut_size_after" << fp.A.size() << endl;
    cout << "local flow: complemen  t siz after" << fp.R.size() << endl;

    fp.active_nodes = priority_queue <node_label_pair>();
    for (auto& n : fp.A) {
        if (fp.node_flow[n] >
            fp.node_cap[n] &&
            fp.node_label[n] < fp.h - 1) {
                node_label_pair next_nlp;
                next_nlp.label = fp.node_label[n];
                next_nlp.node  = n;
                fp.active_nodes.push(next_nlp);
            }
    }

    //this is bad? check instead if one push or relabel was done.
    if (fp.active_nodes.empty()) { //(cut_size_before == fp.A.size()) {
        cout << "unit flow done. pushed ... flow:" << pushed_flow << endl;
        assert(level_queue[fp.h - 1].empty());
        for (auto& n : fp.A) {
            //Q: this should hold
            assert(fp.node_flow[n] <= fp.node_cap[n] || fp.node_label[n] == fp.h - 1);
        }
        return;
    }


    //fp.h = fp.A.size();

    goto unit_start;
}

auto
slow_trimming(GraphContext& gc, const Configuration&  /*conf*/, set<Node> cut, double phi) -> set<Node>
{
    // Q: why? is this necessary or assumption?
    // assert(cut.size() >= static_cast<int>(gc.nodes.size()) / 2);
    DG sg;
    DGNodeMap<Node> sg_to_orig(sg);
    NodeMap<DGNode> orig_to_sg(gc.g);
    ArcMap<int> cap(sg, 0.);
    NodeMap<int> orig_degree(gc.g, 0);

    // reverse definition
    NodeMap<bool> in_cut(gc.g, false);
    for (auto& n : cut) {
        in_cut[n] = true;
    }

    for (NodeIt n(gc.g); n != INVALID; ++n) {
        if (in_cut[n]) {
            DGNode sn = sg.addNode();
            sg_to_orig[sn] = n;
            orig_to_sg[n] = sn;
        }
    }

    DGNode s = sg.addNode();
    DGNode t = sg.addNode();
    int cut_edges = 0;
    for (EdgeIt e(gc.g); e != INVALID; ++e) {
        Node u = gc.g.u(e);
        Node v = gc.g.v(e);
        DGNode sg_u = orig_to_sg[u];
        DGNode sg_v = orig_to_sg[v];
        orig_degree[u] += 1;
        orig_degree[v] += 1;
        if (u == v) {
            continue;
}
        if (in_cut[u] && in_cut[v]) {
            Arc e1 = sg.addArc(sg_u, sg_v);
            Arc e2 = sg.addArc(sg_v, sg_u);
            cap[e1] = 2. / phi;
            cap[e2] = 2. / phi;
        }
        else if (in_cut[u] && !in_cut[v]) {
            Arc e = sg.addArc(s, sg_u);
            cap[e] = 2. / phi;
            cut_edges++;
        }
        else if (!in_cut[u] && in_cut[v]) {
            Arc e = sg.addArc(s, sg_v);
            cap[e] = 2. / phi;
            cut_edges++;
        }
        else {
            continue;
        }
    }
    for (NodeIt n(gc.g); n != INVALID; ++n) {
        if (in_cut[n]) {
            Arc e = sg.addArc(orig_to_sg[n], t);
            cap[e] = orig_degree[n];
            assert(cap[e] > 0);
        }
}

    Preflow<DG> preflow(sg, cap, s, t);
    preflow.init();
    preflow.run();
    cout << "local flow: cut size before: " << cut.size() << endl;

    Bfs<DG> bfs(sg);
    bfs.run(s);

    cout << "IDs for s, t are... " << DG::id(s) << ", " << DG::id(t) << endl;
    for (DGNodeIt n(sg); n != INVALID; ++n) {
        if (n != t && n != s) {
            ;  // cout << "local flow: Reached n: " << gc.g.id(sg_to_orig[n]) <<
}
}
               // " ? "
               // << bfs.reached(n) << endl;
    // assert(bfs.reached(t));
    if (!bfs.reached(t)) {
        cout << "local flow: WARNING - A not connected" << endl;
        // assert (false && "require A connected for trimming");
        // assert(!connected(sg));
    }

    cut.clear();

    for (DGNodeIt n(sg); n != INVALID; ++n) {
        //for (DGOutArcIt e(sg, n); e != INVALID; ++e) {
        //    out_flow += preflow.flow(e);
        //}
        //if (out_flow != 0) {
        //    cout << "local flow: arc flow" << out_flow << " n: " << sg.id(n)  << endl;
        //}
        // cout << preflow.minCut(n) << " : " << sg.id(n) << endl;
        if (!preflow.minCut(n) && n != s && n != t) {  //) // && n != s && n != t)
            cut.insert(sg_to_orig[n]);
}
            //cout << "local flow: slow trim: " << endl;
    }
    for (DGNodeIt n(sg); n != INVALID; ++n) {
        ;  // cout << (preflow.
    }

    cout << "local flow: cut size after " << cut.size() << endl;
    cout << "local flow: maximum preflow: " << preflow.flowValue() << endl;
    //Q: * 2 after reduction?
    cout << "local flow: sought preflow: " << 2 * cut_edges / phi << endl;
    cout << "local flow: (phi = ) " << phi << "cut edges " << cut_edges << endl;
    // Sought == maximum, we achieved our goal.

    return cut;
}

auto
trimming(GraphContext& gc, const Configuration& conf, set<Node> cut, double phi) -> set<Node>
{
    // assert (phi >= 0.00001);
    if (cut.size() == gc.nodes.size()) {
        return cut;
    } 
    assert(cut.size() >= gc.nodes.size() / 2);
    DG dg;
    digraph_from_graph(gc.g, dg);

    flow_instance fp = flow_instance(dg, phi);
    fp.n = static_cast<int>(gc.nodes.size());
    fp.e = gc.num_edges;
    set<Node>   R_ug;
    set<DGNode> R;
    set<DGNode> A;

    // Setup larger side A and complement
    for (auto& n : cut) {
        A.insert(DG::nodeFromId(G::id(n)));
}

    std::set_difference(gc.nodes.begin(), gc.nodes.end(), cut.begin(),
                        cut.end(), std::inserter(R_ug, R_ug.end()));
    // cut.swap(R_ug);

    for (auto& n : R_ug) {
        R.insert(DG::nodeFromId(G::id(n)));
}

    assert(R.size() <= A.size());

    fp.A = A;
    fp.R = R;
    fp.h = 40 * log(2 * gc.num_edges) / phi;
    fp.h = A.size();
    fp.phi = phi;
    // Wasteful but still mlogm
    // Crossing edge sources

    // set sink capacities
    // Q: equals degree?
    for (auto& n : fp.A) {
        fp.node_cap[n] = 0.0;
        fp.node_flow[n] = 0.0;
        assert(R.count(n) == 0);
        for (DGOutArcIt e(dg, n); e != INVALID; ++e)
            if (A.count(dg.target(e)) && A.count(dg.source(e))) {
                fp.node_cap[n] += 1;
                fp.edge_cap[e] = 2.0 / phi;
            }
        fp.node_label[n] = 0;
        fp.node_flow[n] = 0.0;
        assert(fp.node_cap[n] >= 0.0);
    }
    // for (ListDigraph::ArcIt e(DG); e != INVALID; ++e) {
    //    fp.edge_cap[e] = 2/phi;
    //}
    // assert(R.size() > 0);
    for (auto& n : fp.R) {
        fp.node_flow[n] = 0.0;
        for (DGOutArcIt e(dg, n); e != INVALID; ++e) {
            if (!R.count(dg.target(e))) {
                fp.edge_cap[e] = 2. / phi;
                fp.node_flow[n] -= 2. / phi;
                fp.edge_flow[e] = 2. / phi;
                
                fp.node_flow[dg.target(e)] += 2. / phi;
                fp.edge_flow[fp.reverse_arc[e]] = -2. / phi;
            }
            else {
                fp.edge_cap[e] = 0;
                fp.edge_flow[e] = 0.0;
            }
        }
        fp.node_label[n] = 0;
        fp.node_cap[n] = 0.0;
    }

    for (auto& n : fp.A) {
        if (fp.node_flow[n] - min(fp.node_flow[n], fp.initial_mass[n]) >
            fp.node_cap[n]) {
                node_label_pair nlp;
                nlp.label = 0;
                nlp.node  = n;
                fp.active_nodes.push(nlp);
            }
    }

    if (fp.active_nodes.size() == 0) {
        cout << "local flow: nothing to do" << endl;
        return cut;
    }
    // assert (fp.active_nodes.size() > 0);

    assert(fp.A.size() >= fp.R.size());
    uf(dg, fp, conf);

    // dummy
    cut.clear();
    for (auto& n : fp.A)
        cut.insert(gc.g.nodeFromId(dg.id(n)));
    return cut;
}

void
graph_from_cut(GraphContext& g, GraphContext& sg, const set<Node>& cut)
{
    map<Node, Node> reverse_map = map<Node, Node>();
    set<Node> all_nodes = set<Node>();
    set<Node> complement_nodes = set<Node>();
    sg.num_edges = 0;
    assert(sg.nodes.empty());

    for (auto n : cut) {
        Node x = sg.g.addNode();
        sg.nodes.push_back(x);
        // Maybe not necessary
        reverse_map[n] = x;
        sg.orig_degree[x] = g.orig_degree[n];
    }
    map<Node, int> self_loops_added_to_node;
    for (NodeIt n(sg.g); n != INVALID; ++n) {
        self_loops_added_to_node[n] = 0;
    }
    //map<Node, bool> added;
    for (const auto& n : cut) {
        for (IncEdgeIt a(g.g, n); a != INVALID; ++a) {
            if (cut.count(g.g.target(a)) > 0 && cut.count(g.g.source(a)) > 0 &&
                G::id(g.g.source(a)) < G::id(g.g.target(a))) {
                //assert(reverse_map[n] !=
                //       reverse_map[g.g.target(a)]);  // no self -loop
                assert(sg.g.id(reverse_map[g.g.source(a)]) < static_cast<int>(sg.nodes.size()) &&
                       sg.g.id(reverse_map[g.g.target(a)]) < static_cast<int>(sg.nodes.size()));
                sg.g.addEdge(reverse_map[g.g.source(a)],
                             reverse_map[g.g.target(a)]);
                sg.num_edges++;
            }
        }
    }
    for (const auto& n : cut) {
        Node nn = reverse_map[n];
        int c = 0;
        for (IncEdgeIt a(sg.g, nn); a != INVALID; ++a) {
            c += 1;
        }
        assert(c <= sg.orig_degree[nn]);
        for (int i = 0; i < sg.orig_degree[nn] - c; i++) {
            sg.g.addEdge(nn, nn);
        }
    }
    sg.num_edges = countEdges(sg.g);

    #if DEBUG == 1
    for (NodeIt n(sg.g); n != INVALID; ++n) {
        int c = 0;
        int self_loops = 0;
        for (IncEdgeIt e(sg.g, n); e != INVALID; ++e) {
            c += 1;
            if (G::id(sg.g.source(e)) == G::id(sg.g.target(e))) {
                self_loops += 1;
            }
        }
        if (sg.orig_degree[n] != c - self_loops/2) {
            cout << "orig degree: " << sg.orig_degree[n] << " " << "new degree: " << c << endl;
            cout << "self loops: " << self_loops << endl;
            cout << "adjusted: " << sg.orig_degree[n] << " == " << c - self_loops/2 << endl;
            assert(false && "subgraph degrees not equal to original degree");
        }
        assert(sg.orig_degree[n] == (c - self_loops/2));
    }
    #endif

    }

auto
graph_from_cut(GraphContext& g, GraphContext& sg, set<Node> cut,
               map<Node, Node> old_map, bool complement = false) -> map<Node, Node>
{
    map<Node, Node> new_map = map<Node, Node>();
    map<Node, Node> reverse_map = map<Node, Node>();
    set<Node> all_nodes = set<Node>();
    set<Node> complement_nodes = set<Node>();
    sg.num_edges = 0;
    assert(sg.nodes.empty());

    if (complement) {
        for (NodeIt n(g.g); n != INVALID; ++n) {
            all_nodes.insert(n);
        }
        // cut.clear();
        std::set_difference(
            all_nodes.begin(), all_nodes.end(), cut.begin(), cut.end(),
            std::inserter(complement_nodes, complement_nodes.end()));
        cut.swap(complement_nodes);
    }

    for (auto n : cut) {
        Node x = sg.g.addNode();
        sg.nodes.push_back(x);
        // Maybe not necessary
        new_map[x] = old_map[n];
        reverse_map[n] = x;
        sg.orig_degree[x] = g.orig_degree[n];
    }
    // sort(sg.nodes.begin(), sg.nodes.end());

    map<Node, int> self_loops_added_to_node;
    for (NodeIt n(sg.g); n != INVALID; ++n) {
        self_loops_added_to_node[n] = 0;
    }
    for (const auto& n : cut) {
        for (IncEdgeIt a(g.g, n); a != INVALID; ++a) {

            if (cut.count(g.g.target(a)) > 0 && cut.count(g.g.source(a)) > 0 &&
                G::id(g.g.source(a)) < G::id(g.g.target(a))) {
                //assert(reverse_map[n] !=
                //       reverse_map[g.g.target(a)]);  // no self -loop
                assert(sg.g.id(reverse_map[g.g.source(a)]) < static_cast<int>(sg.nodes.size()) &&
                       sg.g.id(reverse_map[g.g.target(a)]) < static_cast<int>(sg.nodes.size()));
                sg.g.addEdge(reverse_map[g.g.source(a)],
                             reverse_map[g.g.target(a)]);
                sg.num_edges++;
            }
        }
    }
    for (const auto& n : cut) {
        Node nn = reverse_map[n];
        int c = 0;
        for (IncEdgeIt a(sg.g, nn); a != INVALID; ++a) {
            c += 1;
        }
        assert(c <= sg.orig_degree[nn]);
        for (int i = 0; i < sg.orig_degree[nn] - c; i++) {
            sg.g.addEdge(nn, nn);
        }
    }
    sg.num_edges = countEdges(sg.g);
    #if DEBUG == 1
    for (NodeIt n(sg.g); n != INVALID; ++n) {
        int c = 0;
        int self_loops = 0;
        for (IncEdgeIt e(sg.g, n); e != INVALID; ++e) {
            c += 1;
            if (G::id(sg.g.source(e)) == G::id(sg.g.target(e))) {
                self_loops += 1;
            }
        }
        if (sg.orig_degree[n] != c - self_loops/2) {
            cout << "orig degree: " << sg.orig_degree[n] << " " << "new degree: " << c << endl;
            cout << "self loops: " << self_loops << endl;
            cout << "adjusted: " << sg.orig_degree[n] << " == " << c - self_loops/2 << endl;
            assert(false && "subgraph degrees not equal to original degree");
        }
        //assert(sg.orig_degree[n] == c - self_loops/2);
    }
    #endif
    return new_map;
}

auto
find_connected_components(GraphContext& g) -> vector<set<Node>>
{
    NodeMap<bool> visited_nodes(g.g, false);
    vector<set<Node>> labels;

    Bfs<ListGraph> bfs(g.g);
    bfs.init();

    // bfs.init();
    // bfs.start();

    unsigned long n_visited = 0;
    unsigned long cc = 0;
    Node source;
    Node y;

    // Q:
    unsigned long test_check = 0;
    for (NodeIt n(g.g); n != INVALID; ++n) {
        test_check = test_check + 1;
}
    assert(test_check == g.nodes.size());
    while (n_visited < g.nodes.size()) {
        labels.emplace_back();

        for (NodeIt n(g.g); n != INVALID; ++n) {
            // cout << "Checking node n: " << g.g.id(n) << endl;
            // Start with some unvisited node
            if (!visited_nodes[n]) {
                y = n;
                break;
            }
            // continue;
            // goto return_labels;
        }

        // cout << "Choosing node: " << g.g.id(y) << endl;
        bfs.addSource(y);
        // bfs.run();
        assert(y != INVALID);
        assert(visited_nodes[y] == false);

        // starting node y too
        n_visited++;
        labels[cc].insert(y);
        visited_nodes[y] = true;
        Node s_saved = y;
        while (!bfs.emptyQueue() && y != INVALID) {
            if (!visited_nodes[y]) {
                labels[cc].insert(y);
                visited_nodes[y] = true;
                n_visited++;
            }
            y = bfs.processNextNode();
        }
        assert(y != INVALID);
        if (!visited_nodes[y]) {
            labels[cc].insert(y);
            visited_nodes[y] = true;
            n_visited++;
        }
        assert(labels[cc].count(s_saved) == 1);
        assert(!labels[cc].empty());
        cc++;
    }
    unsigned long tot = 0;
    for (auto& l : labels) {
        tot = tot + l.size();
}
    cout << "N visited " << n_visited << "actual: " << tot << " g node size "
         << g.nodes.size() << " cc:  " << cc << endl;
    // assert (n_visited == g.nodes.size());
    assert(cc <= g.nodes.size());
    assert(tot == g.nodes.size());
    cout << "n components: " << cc << " should equal " << labels.size() << endl;
    return labels;
}

auto
connected_component_from_cut(GraphContext& gc_orig, set<Node> A) -> set<Node>
{
    set<Node> R;

    std::set_difference(gc_orig.nodes.begin(), gc_orig.nodes.end(), A.begin(),
                        A.end(), std::inserter(R, R.end()));
    #if DEBUG == 1
    assert(connected(gc_orig.g));
    vector<set<Node>> test_cc = find_connected_components(gc_orig);
    assert(test_cc.size() == 1);
    #endif
    // Q: just for now?
    // assert(A.size() >= R.size());
    GraphContext A_sg;
    GraphContext R_sg;

    map<Node, Node> gc_to_gc;
    for (auto& n : gc_orig.nodes)
        {gc_to_gc[n] = n;}

    map<Node, Node> A_to_orig = graph_from_cut(gc_orig, A_sg, A, gc_to_gc);
    map<Node, Node> R_to_orig = graph_from_cut(gc_orig, R_sg, R, gc_to_gc);

    vector<set<Node>> components_A = find_connected_components(A_sg);
    vector<set<Node>> components_R = find_connected_components(R_sg);

    // bool longest(const std::vector<int>& lhs, const std::vector<int>& rhs) {
    //    return lhs.size() < rhs.size();
    //}
    cout << "orig graph: n: " << gc_orig.nodes.size() << endl;

    // auto comp = [](int x, int y){ return x < y; };
    // auto set_comp ()(set<Node> &a, set<Node> &b ) { return a.size() <
    // b.size();
    // };
    auto comp = [](set<Node>& a, set<Node>& b) -> bool {
        return a.size() < b.size();
    };

    while (components_A.size() > 1 || components_R.size() > 1) {
        cout << "components not connected: n components in A: "
             << components_A.size() << " n in R: " << components_R.size()
             << endl;

        if (components_A.size() > 1) {
            for (auto c = components_A.begin(); c != components_A.end(); c++) {
                // R.insert((*c).begin(), (*c).end());
                if (c == std::max_element(components_A.begin(),
                                          components_A.end(), comp)) {
                    continue;
                }
                for (auto& n : *c) {
                    R.insert(A_to_orig[n]);
}
            }
            A.clear();
            for (auto& n : *std::max_element(components_A.begin(),
                                             components_A.end(), comp)) {
                A.insert(A_to_orig[n]);
}
        }
        else if (components_R.size() > 1) {
            for (auto c = components_R.begin(); c != components_R.end(); c++) {
                if (c == std::max_element(components_R.begin(),
                                          components_R.end(), comp)) {
                    continue;
                }
                for (auto& n : *c) {
                    A.insert(R_to_orig[n]);
}
            }
            R.clear();
            for (auto& n : *std::max_element(components_R.begin(),
                                             components_R.end(), comp)) {
                R.insert(R_to_orig[n]);
}
        }

        assert(A.size() + R.size() == gc_orig.nodes.size());
        GraphContext A_sg_;
        GraphContext R_sg_;
        A_to_orig = graph_from_cut(gc_orig, A_sg_, A, gc_to_gc);
        R_to_orig = graph_from_cut(gc_orig, R_sg_, R, gc_to_gc);
        assert(A_sg_.nodes.size() + R_sg_.nodes.size() == gc_orig.nodes.size());

        components_A = find_connected_components(A_sg_);
        components_R = find_connected_components(R_sg_);
        assert(components_A.size() == 1 || components_R.size() == 1);
    }
    assert(A.size() + R.size() == gc_orig.nodes.size());

    if (A.size() >= R.size()) {
        return A;
}
    return R;
}

//Q: not used. should be min element
auto
cut_from_cm(GraphContext& gc, const Configuration& config) -> set<Node>
{
    cout << "init cut-matching" << endl;
    default_random_engine random_engine = config.seed_randomness
                                              ? default_random_engine(
                                                    config.seed)
                                              : default_random_engine(
                                                    random_device()());
    CutMatching cm(gc, config, random_engine);
    cout << "Start cut-matching" << endl;
    cm.run();
    cout << "Finish cut-matching" << endl;
    assert(!cm.sub_past_rounds.empty());
    auto& best_round = *max_element(
        cm.sub_past_rounds.begin(), cm.sub_past_rounds.end(),
        [](auto& a, auto& b) { return a->g_conductance < b->g_conductance; });

    return *(best_round->cut);
}

auto
standalone_conductance(GraphContext& gc, const set<Node>& cut) -> double
{
    int edges_interior = 0;
    int cross_edges = 0;
    int complement_edges = 0;
    for (G::ArcIt e(gc.g); e != INVALID; ++e) {
        if ((cut.count(gc.g.source(e)) != 0U) && (cut.count(gc.g.target(e)) != 0U)) {
            edges_interior++;
        } else if (((cut.count(gc.g.source(e)) != 0U) && (cut.count(gc.g.target(e)) == 0U)) ||
                 ((cut.count(gc.g.source(e)) == 0U) && (cut.count(gc.g.target(e)) != 0U))) {
            cross_edges++;
        } else {
            complement_edges++;
}
    }
    cout << "local flow: interior e: " << edges_interior
         << " cross: " << cross_edges << " comp: " << complement_edges << endl;
    cout << "local flow: conductance: "
         << (1.0 * cross_edges) /
                (cross_edges + min(edges_interior, complement_edges))
         << endl;
    return (1.0 * cross_edges) /
           (cross_edges + min(edges_interior, complement_edges));
}

struct cm_result {
    set<Node> best_cut;
    set<Node> last_cut;
    double best_conductance{};
    double last_conductance{};
    bool reached_H_target{};
    bool best_relatively_balanced{};
    bool last_relatively_balanced{};
    double best_volume_ratio{};
    double last_volume_ratio{};
    int best_cut_crossing_edges{};
};

void
run_cut_matching(GraphContext& gc, Configuration& config, cm_result& cm_res)
{
    assert(connected(gc.g));
    assert(static_cast<int>(gc.nodes.size()) > 2);

    Node temp_node;
    Edge temp_edge;

    bool added_node = false;
    if (static_cast<int>(gc.nodes.size()) % 2 != 0) {
        added_node = true;
        temp_node = gc.g.addNode();
        temp_edge = gc.g.addEdge(gc.nodes[static_cast<int>(gc.nodes.size()) - 2], temp_node);
        gc.nodes.push_back(temp_node);
        gc.num_edges++;
    }

    cout << "init cut-matching" << endl;

    default_random_engine random_engine = config.seed_randomness
                                              ? default_random_engine(
                                                    config.seed)
                                              : default_random_engine(
                                                    random_device()());
    CutMatching cm(gc, config, random_engine);
    cout << "Start cut-matching" << endl;
    cm.run();
    cout << "Finish cut-matching" << endl;
    assert(!cm.sub_past_rounds.empty());
    auto& best_round = *min_element(
        cm.sub_past_rounds.begin(), cm.sub_past_rounds.end(),
        [](auto& a, auto& b) { return a->g_conductance < b->g_conductance; });
    // Q: hack
    auto& last_round = *min_element(
        cm.sub_past_rounds.end() - 1, cm.sub_past_rounds.end(),
        [](auto& a, auto& b) { return a->g_conductance == b->g_conductance; });

    if (added_node) {
        (*(best_round->cut)).erase(temp_node);
        (*(last_round->cut)).erase(temp_node);
        gc.g.erase(temp_edge);
        gc.g.erase(temp_node);
        gc.nodes.pop_back();
        gc.num_edges--;
    }

    cout << "The best with highest expansion was found on round"
         << best_round->index << endl;
    cout << "Best cut sparsity: " << endl;
    auto& best_cut = best_round->cut;
    CutStats<G> cs = CutStats<G>(gc.g, static_cast<int>(gc.nodes.size()), *best_cut);
    cs.print();
    int crossing_edges = cs.crossing_edges;

    cm_result cms;

    // cut = cm.reached_H_target == true ? cut = (*(best_round->cut)) :
    // (*(last_round->cut));

    cm_res.best_cut = *best_round->cut;
    cm_res.last_cut = *last_round->cut;
    cm_res.best_conductance = best_round->g_conductance;
    cm_res.last_conductance = last_round->g_conductance;
    cm_res.reached_H_target = cm.reached_H_target;
    cm_res.best_relatively_balanced = best_round->relatively_balanced;
    cm_res.last_relatively_balanced = last_round->relatively_balanced;
    cm_res.best_volume_ratio = best_round->volume;
    cm_res.last_volume_ratio = best_round->volume;
    cm_res.best_cut_crossing_edges = crossing_edges;
}

auto
test_subgraph_expansion(GraphContext& gc, Configuration config, set<Node> cut,
                        double acceptance_ratio) -> bool
{
    cm_result cm_g;
    cm_result cm_sg;

    // Practically unreachable
    double saved_phi = config.G_phi_target;
    config.G_phi_target = PHI_UNREACHABLE;
    run_cut_matching(gc, config, cm_g);

    GraphContext sg_dummy;
    graph_from_cut(gc, sg_dummy, cut);
    if (!connected(sg_dummy.g)) {
        cut = connected_component_from_cut(gc, cut);
    }
    GraphContext sg;
    graph_from_cut(gc, sg, cut);

    config.G_phi_target = saved_phi * PHI_ACCEPTANCE_TRIMMING;
    run_cut_matching(gc, config, cm_sg);

    config.G_phi_target = saved_phi;
    cout << "g expansion: " << cm_g.best_conductance
         << " sg_expansion: " << cm_sg.best_conductance << endl;
    return cm_sg.best_conductance >= cm_g.best_conductance * acceptance_ratio;
}

auto
cut_volume(GraphContext& gc, const set<Node>& cut) -> int
{
    int cut_volume = 0;

    for (EdgeIt e(gc.g); e != INVALID; ++e) {

        if (cut.count(gc.g.u(e)) != 0U) {
            cut_volume += 1;
}
        if (gc.g.u(e) == gc.g.v(e)) {
            continue;
}
        if (cut.count(gc.g.v(e)) != 0U) {
            cut_volume += 1;
}
    }

    cout << "gc.num_edges " << gc.num_edges << " edges and cut volume "
         << cut_volume << endl;

    assert(cut.size() <= gc.nodes.size());
    assert(cut_volume <= gc.num_edges * 2);
    return cut_volume;
}

struct decomp_trace {
    int depth = 0;
    int size = 0;
    double cut_vol_ratio = -1;
};

struct decomp_stats {
    int n_decomps = 0;
    int n_cm = 0;
    int n_trims = 0;
    double time_in_cm = 0.;
    double time_in_fl = 0.;
    double time_in_fl_tests = 0.;
    vector<decomp_trace> traces;
};

auto
decomp(GraphContext& gc, Configuration config,
       map<Node, Node> map_to_original_graph,
       vector<map<Node, Node>> node_maps_to_original_graph, decomp_stats& stats,
       decomp_trace trace) -> vector<map<Node, Node>>
{
    trace.depth += 1;
    trace.size = static_cast<int>(gc.nodes.size());
    stats.n_decomps += 1;
    cout << "####" << stats.n_decomps << endl;

    if (gc.nodes.empty()) {
        return node_maps_to_original_graph;
}

    if (static_cast<int>(gc.nodes.size()) == 1) {
        node_maps_to_original_graph.push_back(map_to_original_graph);
        return node_maps_to_original_graph;
    }

    if (!(connected(gc.g))) {
        vector<set<Node>> labels = find_connected_components(gc);
        int node_cnt = 0;
        for (auto& sg_cut : labels) {
            cout << "decomp on component: n: " << sg_cut.size() << endl;
            /*
            for (auto& n : sg_cut)
                cout << gc.g.id(map_to_original_graph[n]) << " ";
            */
            cout << endl;
            GraphContext sg;
            map<Node, Node> comp_map = graph_from_cut(gc, sg, sg_cut,
                                                      map_to_original_graph);
            vector<map<Node, Node>> decomp_map;
            vector<map<Node, Node>> empty_map;
            decomp_map = decomp(sg, config, comp_map, empty_map, stats, trace);
            //~sg.g();
            node_maps_to_original_graph.insert(
                node_maps_to_original_graph.end(), decomp_map.begin(),
                decomp_map.end());
            node_cnt = node_cnt + static_cast<int>(sg.nodes.size());
        }
        return node_maps_to_original_graph;
    }

    #if DEBUG == 1
    assert(connected(gc.g));
    #endif

    // TOFIX check expansion
    if (static_cast<int>(gc.nodes.size()) == 2) {
        map<Node, Node> map_1 = {
            {G::nodeFromId(0), map_to_original_graph[G::nodeFromId(0)]}};
        map<Node, Node> map_2 = {
            {G::nodeFromId(1), map_to_original_graph[G::nodeFromId(1)]}};
        node_maps_to_original_graph.push_back(map_1);
        node_maps_to_original_graph.push_back(map_2);
        return node_maps_to_original_graph;
    }

    cm_result cm_res;
    auto start = now();
    run_cut_matching(gc, config, cm_res);
    auto stop = now();
    stats.time_in_cm += duration_sec(start, stop);
    stats.n_cm += 1;

    bool balanced = false;
    bool cut_is_good = false;
    set<Node> cut;


    cut = cm_res.best_cut;

    balanced = cm_res.best_relatively_balanced;
    cut_is_good = cm_res.best_conductance < config.G_phi_target;


    // Q: should be an improvement. But this will affect balance. re-calculate?
    //TOOD: use this heuristic
    cut = connected_component_from_cut(gc, cut);

    int cut_vol = cut_volume(gc, cut);

    double h = double(gc.num_edges) / (10. *
                                       pow(log(gc.num_edges), 2));

    double upper = (1. * gc.num_edges - h) / (1. * gc.num_edges);
    double lower = 1. - upper;

    if ((((1. - config.h_ratio) * gc.num_edges) >= (cut_vol / 2.)) &&
        ((cut_vol / 2.) >= (config.h_ratio * gc.num_edges))) {
        balanced = true;
    }
    else if ((config.h_ratio == 0.) && ((gc.num_edges - h) >= (cut_vol / 2.)) &&
             ((cut_vol / 2.) >= h)) {
        balanced = true;
    }
    else {
        balanced = false;
    }

    cout << "cut vol: " << cut_vol << " gc num edges " << gc.num_edges
         << " balanced?: " << balanced << endl;
    cout << "h ratio " << config.h_ratio << " "
         << " h: " << h << " |E| " << gc.num_edges << endl;
    cout << "upper/lower balance threshold" <<  upper << " / " << lower << endl;


    trace.cut_vol_ratio = 1.0 * cut_vol / (2. * gc.num_edges); 
    trace.cut_vol_ratio = trace.cut_vol_ratio < 0.5? 1 - trace.cut_vol_ratio: trace.cut_vol_ratio;
    trace.cut_vol_ratio = cut_is_good ? trace.cut_vol_ratio : -1;

    assert(trace.cut_vol_ratio <= 1.0);

    if (!cut_is_good) {
        cout << "CASE1 NO Goodenough cut (timeout), G certified expander."
             << endl;
        assert(trace.cut_vol_ratio == -1);
        node_maps_to_original_graph.push_back(map_to_original_graph);

        #if DEBUG == 1
        //GraphContext debug_sg_1;
        //graph_from_cut(gc, debug_sg_1, cut);
        assert(connected(gc.g));
        #endif
    }

    // break early due to balanced and good cut
    else if (cut_is_good && balanced) {
        assert(!cut.empty() != static_cast<int>(gc.nodes.size()));
        int t = NUM_THREADS;
        int edge_count = 0;
	    bool done_recurse = false;
        (void) done_recurse;
#pragma omp parallel for num_threads(t) schedule(dynamic, 1)
        for (int i = 0; i < 2; i++) {
            int thread_num = omp_get_thread_num();
            int cpu_num = sched_getcpu();
            printf("parallel: Thread %3d is running on CPU %3d\n", thread_num,
                   cpu_num);

            GraphContext A;

            bool swtch = i == 1;
            map<Node, Node> new_map = graph_from_cut(
                gc, A, cut, map_to_original_graph, swtch);
            // TODO: do this but adjust for self-loops.
            /*
            if (i == 0) {
                cout << "A.num edges" << A.num_edges << endl;
                cout << "Cut vol /2 " << cut_vol / 2<< endl;
                cout << "best cut crossing" << cm_res.best_cut_crossing_edges << endl;
                assert (A.num_edges == cut_vol/ 2);
            }
            */
            edge_count += A.num_edges;
            cout << "ec" << edge_count << " i " << i << endl;
            // assert (A.nodes.size() == cut.size());

            vector<map<Node, Node>> empty_map;
            vector<map<Node, Node>> decomp_map = decomp(
                A, config, new_map, empty_map, stats, trace);

#pragma omp critical
            node_maps_to_original_graph.insert(
                node_maps_to_original_graph.end(), decomp_map.begin(),
                decomp_map.end());
            if (i == 1) {
                done_recurse = true;
}
            #if DEBUG == 1
            //Doesn't have to be true
            ; //assert(connected(A.g));
            #endif
        }

    }

    // check best cut found so far
    else if (!balanced && cut_is_good) {
        assert(!cut.empty());

        if (cut.size() < gc.nodes.size() / 2) {
            set<Node> t;
            std::set_difference(gc.nodes.begin(), gc.nodes.end(), cut.begin(),
                                cut.end(), std::inserter(t, t.begin()));
            cut = t;
        }

        cout << "Call local flow with cut of size: " << cut.size() << " / "
             << static_cast<int>(gc.nodes.size()) << endl;
        cout << "local flow: conductance before: " << endl;
        standalone_conductance(gc, cut);
        set<Node> cut_saved = cut;
        cut = slow_trimming(gc, config, cut, config.G_phi_target);
        auto start2 = now();
        set<Node> real_trim_cut = trimming(gc, config, cut_saved,
                                           config.G_phi_target);
        auto stop2 = now();
        stats.time_in_fl += duration_sec(start2, stop2);
        stats.n_trims += 1;

        cout << "After local flow, cut is reduced to n nodes: " << cut.size()
             << endl;

        cout << "After local flow (real), cut is reduced to n nodes: "
             << real_trim_cut.size() << endl;


        GraphContext V_over_A;
        GraphContext A;
        map<Node, Node> R_map = graph_from_cut(gc, V_over_A, real_trim_cut,
                                               map_to_original_graph, true);
        map<Node, Node> A_map = graph_from_cut(gc, A, real_trim_cut,
                                               map_to_original_graph, false);
 
        #if DEBUG == 1
        if (real_trim_cut.size() < cut.size() && config.decompose_with_tests) { 
            auto start3 = now();
            bool sg_is_expander = test_subgraph_expansion(gc, config, real_trim_cut,
                                                        PHI_ACCEPTANCE_TRIMMING);


            warn(connected(A.g), "subgraph is not connected after trimming!", __LINE__);
            warn(sg_is_expander, "sg is not an expander after trimming!", __LINE__);
            auto stop3 = now();
            stats.time_in_fl_tests += duration_sec(start3, stop3);
        }
        #endif


        cut = real_trim_cut;
        assert(A.nodes.size() + V_over_A.nodes.size() == gc.nodes.size());
        //Q: is this a real error?
        warn(!V_over_A.nodes.empty(), "non trimmed side is 0!", __LINE__);
        warn(cut.empty(), "trimmed component has size 0", __LINE__);

        node_maps_to_original_graph.push_back(A_map);
        vector<map<Node, Node>> empty_map;
        vector<map<Node, Node>> cuts_empty_map;
        vector<map<Node, Node>> decomp_map = decomp(
            V_over_A, config, R_map, cuts_empty_map, stats, trace);
        #if DEBUG == 1
        assert(connected(A.g));
        #endif
        node_maps_to_original_graph.insert(node_maps_to_original_graph.end(),
                                           decomp_map.begin(),
                                           decomp_map.end());
    }

    else {  // provisional round limit
        cout << "CASE(4) no good cut, round break, G certified nothing."
             << endl;
        node_maps_to_original_graph.push_back(map_to_original_graph);
        assert(false);
    }
    stats.traces.push_back(trace);
    return node_maps_to_original_graph;
}

void
butterfly_test(const Configuration& config, int ls, int rs, int conn_ls, int conn_rs)
{
    GraphContext gc;
    for (int i = 0; i < ls + rs + 1; i++) {
        Node x = gc.g.addNode();
        gc.nodes.push_back(x);
    }
    for (int i = 0; i < ls; i++) {
        for (int j = i + 1; j < ls; j++) {
            gc.g.addEdge(G::nodeFromId(i), G::nodeFromId(j));
            gc.num_edges++;
        }
        if (i < conn_ls) {
            gc.g.addEdge(G::nodeFromId(i), G::nodeFromId(ls + rs));
            gc.num_edges++;
        }
    }
    for (int i = ls; i < ls + rs; i++) {

        for (int j = i + 1; j < ls + rs; j++) {
            gc.g.addEdge(G::nodeFromId(i), G::nodeFromId(j));
            gc.num_edges++;
        }
        if (i < conn_rs + ls) {
            gc.g.addEdge(G::nodeFromId(i), G::nodeFromId(ls + rs));
            gc.num_edges++;
        }
    }

    for (NodeIt n(gc.g); n != INVALID; ++n) {
        int ec = 0;
        for (OutArcIt e(gc.g, n); e != INVALID; ++e) {
            ec++;
            break;
        }
        assert(ec);
    }

    set<Node> cut;
    for (int i = 0; i < ls + rs; i++) {
        cut.insert(G::nodeFromId(i));
    }
    assert(connected(gc.g));

    assert(cut.size() == gc.nodes.size() - 1);
    GraphContext sg;
    graph_from_cut(gc, sg, cut);

    assert(!connected(sg.g));

    set<Node> new_cut = slow_trimming(gc, config, cut,
                                      1. / pow(static_cast<int>(gc.nodes.size()), 2));
    GraphContext sg_2;
    graph_from_cut(gc, sg_2, new_cut);
    assert(connected(sg_2.g));
    new_cut = trimming(gc, config, cut, 1. / pow(static_cast<int>(gc.nodes.size()), 2));
    for (auto& n : new_cut) {
        cout << G::id(n) << endl;
    }
    GraphContext sg_3;
    graph_from_cut(gc, sg_3, new_cut);
    assert(connected(sg_3.g));
}

void
test_connect_subgraph(GraphContext& gc, Configuration conf, double phi)
{
    GraphContext sg;
    GraphContext sg_slow;
    GraphContext sg_fast;
    vector<int> indices;
    set<Node> cut;

    conf.G_phi_target = PHI_UNREACHABLE;

    /*
    for (int i = 0; i < static_cast<int>(gc.nodes.size()); i++) {
        indices.push_back(i);
    }
    */
    

    // Q: use parameter or set dummy
    cm_result cm_res;
    if (conf.known_phi == 0.0) {
        run_cut_matching(gc, conf, cm_res);
        phi = cm_res.best_conductance * 0.8;
    }
    else {
        phi = conf.known_phi;
        cm_res.best_conductance = phi;
    }


    //double ratio_to_remove = 1./static_cast<int>(gc.nodes.size());

    double h_ratio = conf.h_ratio;
    if (h_ratio <= 0.) {
        h_ratio = (double(gc.num_edges) /
                    (10. * conf.volume_treshold_factor *
                     pow(log(gc.num_edges), 2))) / static_cast<int>(gc.nodes.size());
    cout << "local flow:, set h_ratio to " << h_ratio << endl;
    }


    uniform_real_distribution<double> random_fraction(0, h_ratio);
    random_device rd;
    default_random_engine r_engine(rd());

    int iter = 0;
    do {
        sg.clear();
        cut.clear();
        indices.clear();

        double node_fraction_to_remove = random_fraction(r_engine);
        cout << node_fraction_to_remove << " fraction" << endl;

        for (int i = 0; i < static_cast<int>(gc.nodes.size()); i++) {
            indices.push_back(i);
        }

        shuffle(indices.begin(), indices.end(),
                default_random_engine(random_device()()));

        indices.erase(indices.begin(),
                      indices.begin() + static_cast<int>(gc.nodes.size()) * node_fraction_to_remove);

        for (auto& i : indices) {
            cut.insert(G::nodeFromId(i));
        }

        cout << indices[0] << " indices 0" << endl;
        cout << "nodes in cut" << indices.size() << endl;

        assert(gc.nodes.size() >= cut.size() && cut.size() > 0);

        graph_from_cut(gc, sg, cut);
        iter++;

    } while (iter < 1); //while (connected(sg.g));


    //phi = 1. / (gc.num_edges);

    //assert(!connected(sg.g));
    cout << "nodes orig" << static_cast<int>(gc.nodes.size()) << endl;
    cout << "nodes cut" << cut.size() << endl;
    set<Node> scut = slow_trimming(gc, conf, cut, phi);
    set<Node> scut2 = slow_trimming(gc, conf, scut, phi);
    set<Node> fcut = trimming(gc, conf, cut, phi);

    graph_from_cut(gc, sg_slow, scut2);
    graph_from_cut(gc, sg_fast, fcut);
    // assert(sg_slow.num_edges == sg_fast.num_edges > 0);
    // assert(static_cast<int>(gc.nodes.size()) >= sg_slow.nodes.size() > 0);

    cout << "connected slow?" << connected(sg_slow.g) << endl;
    cout << "connected fast?" << connected(sg_fast.g) << endl;
    assert(connected(sg_slow.g));
    assert(connected(sg_fast.g));
    // Q: clear sg

    cm_result cm_res_slow;
    cm_result cm_res_fast;
    if (sg_slow.nodes.size() > 1) {
        run_cut_matching(sg_slow, conf, cm_res_slow);
    } else {
        cm_res_slow.best_conductance = 1;        
}
    if (sg_fast.nodes.size() > 1) {
        run_cut_matching(sg_fast, conf, cm_res_fast);
    } else {
        cm_res_fast.best_conductance = 1;
}
    
    double h = double(gc.num_edges) /
                    (10. * 1. *
                     pow(log(1. * gc.num_edges), 2));

    cout << "===" << endl;
    if (connected(sg.g)) {
        cm_result cm_res_sub;
        run_cut_matching(sg, conf, cm_res_sub);
        cout << "exptest: lowest conductance after (after cut, before trim) <-- " << cm_res_sub.best_conductance << " cut size " << static_cast<int>(sg.nodes.size()) << "/" << static_cast<int>(gc.nodes.size()) << endl;
    }
    cout << "phi used (0.8 * best before * 0.16) " << phi << endl;
    cout << "exptest: graph connected after cut? " << connected(sg.g) << endl;
    cout << "exptest: lowest conductance before       <-- " << cm_res.best_conductance      << " cut size " << indices.size()       << "/" << static_cast<int>(gc.nodes.size()) << " (in theory, we should only remove h = " << h << " nodes)" << endl;
    cout << "exptest: lowest conductance after (slow) <-- " << cm_res_slow.best_conductance << " cut size " << sg_slow.nodes.size() << "/" << static_cast<int>(gc.nodes.size()) << endl;
    cout << "exptest: lowest conductance after (fast) <-- " << cm_res_fast.best_conductance << " cut size " << sg_fast.nodes.size() << "/" << static_cast<int>(gc.nodes.size()) << endl; 
    cout << "===" << endl;
}

void
test_expander_ratio(GraphContext& gc, Configuration conf, double phi_ratio_target)
{

    conf.G_phi_target = PHI_UNREACHABLE;

    // Q: use parameter or set dummy
    cm_result cm_res;
    run_cut_matching(gc, conf, cm_res);
    //Graph is a phi expander. Target phi * ratio target
    double phi = cm_res.best_conductance;

    double phi_thres = phi;

    set<Node> non_expander_cut;

    bool break_var = false;
    int some_big_int = 100000;
    //Should parallelize the function call, instead.
#pragma omp parallel for  num_threads(NUM_THREADS) schedule(dynamic, 1)
    for (int i = 0; i < some_big_int; i += 1) {
        GraphContext sg;
        set<Node> cut;
        sg.clear();
        cut.clear();

        vector<int> indices;
        for (int i = 0; i < static_cast<int>(gc.nodes.size()); i++) {
            indices.push_back(i);
        }

        shuffle(indices.begin(), indices.end(),
                default_random_engine(random_device()()));
        indices.erase(indices.begin(),
                      indices.begin() + static_cast<int>(gc.nodes.size()) * conf.h_ratio);

        for (auto& i : indices) {
            cut.insert(G::nodeFromId(i));
        }

        assert(gc.nodes.size() >= cut.size() && cut.size() > 0);

        graph_from_cut(gc, sg, cut);
        //Test connectivity elsewhere
        if (!connected(sg.g)) {
            continue;
}

        cm_result cm_res;
        run_cut_matching(sg, conf, cm_res);
        phi = cm_res.best_conductance;

        i = 0;
        if (phi_thres * phi_ratio_target > phi) {
            break_var = true;
}
            non_expander_cut = cut;
        
        if (!break_var) {
            continue;
}
    }

#pragma omp taskwait

    set<Node> cut = non_expander_cut;
    for (auto& n : cut) {
        cout << G::id(n) << " ";
    }
    cout << endl;

    cout << "nodes orig" << static_cast<int>(gc.nodes.size()) << endl;
    cout << "nodes cut" << cut.size() << endl;
    set<Node> scut  = slow_trimming(gc, conf, cut, phi);
    set<Node> scut2 = slow_trimming(gc, conf, scut, phi);
    set<Node> fcut  = trimming(gc, conf, cut, phi);

    GraphContext sg_slow;
    GraphContext sg_fast;
    graph_from_cut(gc, sg_slow, scut2);
    graph_from_cut(gc, sg_fast, fcut);
    // assert(sg_slow.num_edges == sg_fast.num_edges > 0);
    // assert(static_cast<int>(gc.nodes.size()) >= sg_slow.nodes.size() > 0);

    conf.G_phi_target = phi * phi_ratio_target;
    cm_result cm_res_slow;
    cm_result cm_res_fast;
    run_cut_matching(sg_slow, conf, cm_res);
    run_cut_matching(sg_fast, conf, cm_res);
    phi = cm_res.best_conductance;

    cout << "lowest phi slow?" << cm_res_slow.best_conductance << endl;
    cout << "lowest phi fast?" << cm_res_slow.best_conductance << endl;
    assert(cm_res_slow.best_conductance >= phi * phi_ratio_target);
    assert(cm_res_fast.best_conductance >= phi * phi_ratio_target);
}



auto random_walk_distribution(GraphContext& gc, const set<Node>& cut, long walk_length, int n_trials, int check_convergence_every_n, double tolerance) -> long {
    
    //Q: fix this
    assert(cut.size() == gc.nodes.size());
    if (cut.size() == 1){
        return 1;
    }
    assert (cut.size() > 1);

    uniform_int_distribution<int> random_int(0, static_cast<int>(gc.nodes.size()) - 1);
    random_device rd;
    default_random_engine r_engine(rd());

    //double node_fraction_to_remove = random_fraction(r_engine);
    //cout << node_fraction_to_remove << " fraction" << endl;

    map<Node, int> counts;
    map<Node, vector<Node>> nbors_lst;
    map<Node, unsigned long> nbors;

    vector<double> mixing_rates;

    int tot_c = 0;
    for (NodeIt n(gc.g); n != INVALID; ++n) {
        counts[n] = 0;
        assert (0 <= gc.g.id(n) && gc.g.id(n) < static_cast<int>(gc.nodes.size()));
        nbors_lst[n] = vector<Node>();
        int c = 0;
        for (IncEdgeIt e(gc.g, n); e != INVALID; ++e) {
            assert(gc.g.id(gc.g.target(e)) < static_cast<int>(gc.nodes.size()) && gc.g.id(gc.g.target(e)) >= 0);
            nbors_lst[n].push_back(gc.g.target(e));
            c++;
        }
        nbors[n] = c;
        tot_c += c;
        assert(nbors[n] == nbors_lst[n].size());
    }
    assert (tot_c == gc.num_edges * 2);
    assert (nbors.size() == nbors_lst.size() && nbors_lst.size() == gc.nodes.size());
    int sum_c = 0;
    double old_l1_dist = 1;
    double l1_dist = 1;
    long n_steps = 0;
    for (int i = 0; i < n_trials; i++) {
        n_steps = 0;
        int start_node_index = random_int(r_engine);

        Node curr_node = gc.nodes[start_node_index];
        for (int j = 0; j < walk_length; j++) {
            assert (nbors[curr_node] >= 1);
            if (nbors_lst[curr_node][0] == curr_node) { assert(nbors[curr_node] > 1);
}
            assert(gc.g.id(curr_node) < static_cast<int>(nbors.size()));
            uniform_int_distribution<int> random_nbor_index(0, nbors[curr_node] - 1);
            //if (nbors[curr_node] <= 1) continue;
            random_device rd_;
            default_random_engine r_engine_(rd_());
            unsigned long next_index = random_nbor_index(r_engine_);

            assert (next_index < nbors_lst[curr_node].size());

            counts[curr_node]++;
            assert (next_index <= nbors_lst[curr_node].size());
            curr_node = nbors_lst[curr_node][next_index];
            long saved_steps = n_steps;
            n_steps++;
            //overflows
            assert(n_steps - 1 == saved_steps);
            sum_c++;
 
            double assert_check_max = 0.;
            if (j % check_convergence_every_n == 0) {
                l1_dist = 0.;
                double uniform_expectation = 1./static_cast<int>(gc.nodes.size());
                assert (uniform_expectation > 0 && uniform_expectation <= 1. / (1. * static_cast<int>(gc.nodes.size())));
                double total_visited_nodes = 1. * n_steps;
                double mean_degree = 2. * gc.num_edges / (1. * static_cast<int>(gc.nodes.size()));
                double norm_term = (1. * total_visited_nodes) / (1. *static_cast<int>(gc.nodes.size()));
                assert(uniform_expectation > 0 && total_visited_nodes > 0 && mean_degree > 0 && norm_term > 0);
                for (NodeIt n(gc.g); n != INVALID; ++n) {
                    //double d = abs(1./static_cast<int>(gc.nodes.size()) - ((1. * counts[n]) / (n_trials * walk_length))) * (((1.* gc.num_edges)/ (1.* static_cast<int>(gc.nodes.size()))) / ( 1.*nbors[gc.g.nodeFromId(i)]));
                    double vertex_norm = mean_degree / (1. * nbors[n]);
                    //this can be positive.
                    double ratio_visited_this_node = (1. * counts[n]) / total_visited_nodes;
                    //cout << ratio_visited_this_node << endl;
                    double d_s = uniform_expectation - ((ratio_visited_this_node * vertex_norm)); //  / norm_term);
                    assert_check_max = max(ratio_visited_this_node * vertex_norm, assert_check_max);
                    double d = abs(d_s);

                    assert ((1.*counts[n])/(1.*n_trials)/(1.*walk_length) <= 1);

                    l1_dist += d;
                }
                assert(assert_check_max >= uniform_expectation);

                double delta = static_cast<int>(gc.nodes.size()) * (old_l1_dist - l1_dist) / double(check_convergence_every_n);
                mixing_rates.push_back(delta);
                old_l1_dist = l1_dist;
            }

            if (l1_dist <= tolerance) {

                break;
            }
        }
    }

    double mixings_sum = 0;
    for (auto& m: mixing_rates) {
        mixings_sum += m;
    }

    return n_steps;
}



void
test_expander(GraphContext& gc, const Configuration& conf, double  /*phi*/)
{



    set<Node> full_graph_cut;
    for (NodeIt n(gc.g); n!= INVALID; ++n) {
        full_graph_cut.insert(n);
    }

    //int walk_length = log(static_cast<int>(gc.nodes.size()));
    //int n_trials = 100 * static_cast<int>(gc.nodes.size()) / walk_length;
    int n_trials = 1;
    long walk_length = LONG_MAX;
    //walk_length * n_trials = 100 * gc nodes.size;
    int check_every_n = 1000;
    double tolerance = 0.1;

    if (conf.partition_file.empty()) {
        cout << "start random walk of length " << walk_length << " with n trials: " << n_trials << endl;
        int n_steps = random_walk_distribution(gc, full_graph_cut, walk_length,n_trials, check_every_n, tolerance);
        cout << "steps to stationary convergence:" << n_steps << " steps per node: " << 1. * n_steps / static_cast<int>(gc.nodes.size()) << endl;
        cout << "if == " << walk_length << ", did not converge in time" << endl;
    } else {
        vector<set<Node>> partition;
        partition = read_partition_file(conf.partition_file, gc.nodes, partition);

        int assert_count_nodes = 0;
        for (auto& p : partition) {
            GraphContext sg; 

            graph_from_cut(gc, sg, p);
            assert(connected(sg.g));
            set<Node> dummy_cut;
            for (auto& n: sg.nodes) {dummy_cut.insert(n);}
            int TRIALS = 10;
            int sum_res = 0;
            for (int i = 0; i < TRIALS; i++) {
                long n_steps = random_walk_distribution(sg, dummy_cut, walk_length,n_trials, check_every_n, tolerance);
                sum_res += n_steps;
            }
            int mean_res = (1. * sum_res)/TRIALS;
            cout << "steps to stationary convergence;" << mean_res << ";steps per node;" << mean_res / p.size() << ";partition size;" << p.size() << endl; 
            assert_count_nodes += p.size();

        }
        assert (assert_count_nodes == static_cast<int>(gc.nodes.size()));
    }
    
    }




auto
main(int argc, char** argv) -> int
{

    omp_set_nested(1);
    omp_set_num_threads(NUM_THREADS);

    Configuration config = Configuration();
    parse_options(argc, argv, config);

    if (config.show_help_and_exit) {
        return 0;
}

    GraphContext gc;
    initGraph(gc, config.input);
    //assert(connected(gc.g));
    config.n_nodes_orig = static_cast<int>(gc.nodes.size());
    config.e_edges_orig = gc.num_edges;

    cout << "n: " << static_cast<int>(gc.nodes.size())
         << " e: " << gc.num_edges << endl;
    map<Node, Node> map_to_original_graph;
    for (NodeIt n(gc.g); n != INVALID; ++n) {
        map_to_original_graph[n] = n;
}

    vector<set<Node>> cuts = vector<set<Node>>();
    vector<map<Node, Node>> node_maps_to_original_graph =
        vector<map<Node, Node>>();

    map<Node, int> degree;
    for (NodeIt n(gc.g); n != INVALID; ++n) {
        int d = 0;
        for (IncEdgeIt e(gc.g, n); e != INVALID; ++e) {
            d++;
}
        degree[n] = d;
    }

    // main
    if (static_cast<int>(config.only_test_expander) == 1) {
        test_expander(gc, config, config.G_phi_target);
        return 0;
    }

    decomp_stats stats;
    decomp_trace trace;


    auto start = now();
    vector<map<Node, Node>> cut_maps = decomp(gc, config, map_to_original_graph,
                                              node_maps_to_original_graph,
                                              stats, trace);
    auto stop = now();

    #if DEBUG == 1
    for (auto& cut: cut_maps) {
        set<Node> cut_;
        for (auto& m: cut) {
            cut_.insert(m.second);
        }
        GraphContext sg;
        graph_from_cut(gc, sg, cut_);
        assert(connected(sg.g));
    }
    cout << "all subgraphs verified as connected" << endl;
    #endif

    cout << "Done decomp" << endl;
    cout << "output:" << endl;

    vector<int> all_nodes;

    unsigned long n_singletons = 0;
    unsigned long all_nodes_count = 0;
    vector<vector<Node>> cuts_node_vector;
    for (const auto& m : cut_maps) {
        if (m.size() == 1) {
            for (const auto& c : m) {
                //We never saw it before
                assert(count(all_nodes.begin(), all_nodes.end(),
                             gc.g.id(c.second)) == 0);
                //all_nodes.push_back(gc.g.id(c.second));
            }
            //all_nodes_count++;
            n_singletons++;
            //continue;
        }
        //Comment out above to count singletons correctly below and also output explicit singletones
        cuts_node_vector.emplace_back();
        for (const auto& c : m) {
            assert(count(all_nodes.begin(), all_nodes.end(),
                gc.g.id(c.second)) == 0);
            all_nodes_count++;
            // cout << gc.g.id(c.second) << " ";
            all_nodes.push_back(G::id(c.second));
            cuts_node_vector[cuts_node_vector.size() - 1].push_back(c.second);
        }
    }

    assert(gc.nodes.size() == all_nodes.size());
    assert(all_nodes_count == all_nodes.size());

    // vector<vector<int>> fog(cut_maps.size(), vector<int>(cut_maps.size(),
    // 0));
    vector<double> node_ratio_edges_inside;

    double coms_volume = 0;
    double internal_volume = 0;
    double max_unbalance = 0.;
    int max_size_cluster = 1;
    int min_size_cluster = 1;
    for (auto & i : cuts_node_vector) {
        // for (const auto &m : cuts_node_vector) {
        double edges_outside_cluster = 0.0;
        double edges_inside_cluster = 0.0;
        int all_edges = 0;
        for (const auto& n : i) {
            for (IncEdgeIt e(gc.g, n); e != INVALID; ++e) {
                all_edges += 1;
                if (count(i.begin(), i.end(), gc.g.target(e)) > 0 && count(i.begin(), i.end(), gc.g.source(e)) > 0) {
                    edges_inside_cluster++;
                    internal_volume++;

                }
                else {edges_outside_cluster += 1.0; coms_volume += 1.0;}
        }   }

        //For singletons, we are not double counting
        
        if (i.size() == 1) {
            //edges_inside_cluster++;
            //internal_volume +=1;
            edges_inside_cluster+=1;
            //coms_volume++;
            ;  
        }


        cout << "edges inside cluster/total cluster edges;"
             << edges_inside_cluster << ";" << all_edges << endl;
        node_ratio_edges_inside.push_back(edges_inside_cluster /
                                          static_cast<double>(all_edges));
        max_size_cluster = max(max_size_cluster, all_edges);
        min_size_cluster = min(min_size_cluster, all_edges);
        max_unbalance = max(double(max_unbalance), double((1.*max_size_cluster)/(1.*min_size_cluster)));

    }

    cout << "singletons: " << n_singletons << endl;
    internal_volume = internal_volume / 2;
    coms_volume = coms_volume / 2.;
    cout << "internal volume: " << internal_volume << " coms_volume: " << coms_volume << " num edges: " << gc.num_edges << endl;
    assert(internal_volume + coms_volume == gc.num_edges);

    int i = 0;
    int n_clusters = 0;
    for (auto& c : cuts_node_vector) {
        n_clusters++;
        cout << "decomp cluster " << i << ";";
        for (auto& n : c) {
            cout << G::id(n) << ";";
        }
        cout << endl;
        i++;
    }


    for (auto& t : stats.traces) {
        cout << "vol cut/vol graph ratio;" << t.cut_vol_ratio << endl;
    }
    for (auto& t : stats.traces) {
        cout << "decomp depth;" << t.depth << endl;
    }

    i = 0;

    //NOT a good reason to abort!
    for (const auto& r : node_ratio_edges_inside) {
        cout << "inside cluster vol/total vol cluster;" << r << endl;
        double goal = pow(config.G_phi_target, stats.traces[i].depth);
        //assert (0.0 <= goal && goal <= 1.0);
        //assert ((0.0 < goal && goal < 1.0) || config.G_phi_target == 0 || config.G_phi_target == 1);
        warn(r >= goal, "did not meet conductance goal?", __LINE__);
        i++;
    }

    cout << "graph;" << config.input.file_name << endl;
    cout << "nodes;" << static_cast<int>(gc.nodes.size()) << endl;
    cout << "edges;" << gc.num_edges << endl;
    cout << "g_phi target;" << config.G_phi_target << endl;
    cout << "h_phi target;" << config.H_phi_target << endl;
    cout << "unbalance threshold;" << config.h_ratio << endl;
    cout << "time in cm;" << stats.time_in_cm << endl;
    cout << "time in fl;" << stats.time_in_fl << endl;
    cout << "time testing fl;" << stats.time_in_fl_tests << endl;
    cout << "n local flows;" << stats.n_trims << endl;
    cout << "n cut matches;" << stats.n_cm << endl;
    cout << "n decomps;" << stats.n_decomps << endl;
    cout << "n singletons;" << n_singletons << endl;
    cout << "coms volume;"   << coms_volume << endl;
    cout << "n clusters;" << n_clusters << endl;
    cout << "total time;" << duration_sec(start, stop) << endl;
    cout << "cluster unbalance;" << max_unbalance << endl;
    cout << "output end;" << endl;

    string out = config.input.file_name + "cut.txt";
    ofstream file;
    file.open(out);
    if (!file) {
        cout << "Cannot open file ";  // << OUTPUT_FILE << endl;
        return 1;
    }

    for (const auto& m : cut_maps) {
        for (const auto& c : m) {
            file << G::id(c.second) << " ";
        }
        file << endl;
    }
    file.close();

    return 0;
}
#ifdef __clang__
#pragma clang diagnostic pop
#endif
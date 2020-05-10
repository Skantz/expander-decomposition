import random

import networkx as nx



def cluster_graph(n_clusters, size_clusters, n_crossing_edges=None):



    g = nx.Graph(undirected=True)

    nodes_by_cluster = [[size_clusters*j + i for i in range(size_clusters)] for j in range(n_clusters )]

    for cluster in nodes_by_cluster:
        g.add_nodes_from(cluster)

    for i in range(n_clusters):
        for n1 in nodes_by_cluster[i]:
            for n2 in nodes_by_cluster[i]:
                if n1 != n2:
                    g.add_edge(n1, n2)

    #random.shuffle(cluster_indices)
    #1-indexed
    for i in range(n_clusters - 1):
        g.add_edge(i*size_clusters, (i + 1) * size_clusters)


    print(len(g), n_clusters * size_clusters)
    assert(len(g) == n_clusters * size_clusters)

    return g


def expander_graph(p_nodes): 
    g = nx.chordal_cycle_graph(p_nodes)
    return g

def margulis_graph(n):
    """construct 8-regular marglis graph, |V| = n^2"""
    g = nx.MultiGraph(undirected=True)
    g.add_nodes_from([i for i in range(n**2)])

    adjacency_dict = {}

    #for x in range(n):
    #    for y in range(n):
    #        assert ((x, y) not in adjacency_dict)
    #        adjacency_dict[(x,y)] = []
    #        for pm in [-1, 1]:
    #            adjacency_dict[(x, y)] += [(((x + pm * 2*y)) % n, y)]
    #            adjacency_dict[(x, y)] += [(((x + pm * (2*y + 1)) % n, y))]
    #            adjacency_dict[(x, y)] += [(x, (y + pm * 2*x) % n)]
    #            adjacency_dict[(x, y)] += [(x, (y + pm * (2*x + 1)) % n)]
    for x in range(n):
        for y in range(n):
            assert ((x, y) not in adjacency_dict)
            adjacency_dict[(x,y)] = []
            for pm in [-1, 1]:
                adjacency_dict[(x, y)] += [((x + pm * 2*y) % n, y)]
                adjacency_dict[(x, y)] += [((x + pm * (2*y + 1)) % n, y)]
                adjacency_dict[(x, y)] += [(x, (y + pm * 2*x) % n)]
                adjacency_dict[(x, y)] += [(x, (y + pm * (2*x + 1)) % n)]
    
    for k in adjacency_dict:
        assert(len(adjacency_dict[k]) == 8)
        edges_t = []
        for i, v in enumerate(adjacency_dict[k]):
            s, t = k[0] * n + k[1], v[0]*n + v[1]
            #one-way adjacency representation
            #including self-loops
            if t > s and s in g.neighbors(t):
                continue

            edges_t.append((s, t, dict(k=i)))
            #print(len(list(g.neighbors(k[0]*n + k[1]))))
        g.add_edges_from(edges_t)
        assert(len(list(g.neighbors(s))) <= 8)
        #print(g[s])
        #assert(len(list(g.neighbors(k[0]*n + k[1]))) == 8))

    #assert 8-regular
    #for n in g:
    #    print(len(list(g.neighbors(n))))
    #    #assert(len(list(g.neighbors(n))) == 8)
    return g


def random_d_regular(d, n):
    g = nx.random_regular_graph(d, n)
    return g

def write_graph(G, f):
    G = nx.convert_node_labels_to_integers(G, first_label=1)
    f.write(str(len(G.nodes())) +" "+ str(len(G.edges()))+"\n")
    for line in nx.generate_adjlist(G):
        f.write(line.partition(' ')[2])
        f.write("\n")

f = open("random_3_regular_5000.graph", "w+")
write_graph(random_d_regular(3, 5000), f)

f = open("chordal_cycle_graph_2000.graph", "w+")
write_graph(expander_graph(2000), f)

f = open("margulis_100000.graph", "w+")
write_graph(margulis_graph(int(100000**0.5)), f)

f = open("nxmargulis_10000.graph", "w+")
write_graph(nx.margulis_gabber_galil_graph(int(10000**0.5)), f)

f = open("cluster_graph_50_15.graph", "w+")
write_graph(cluster_graph(5, 50, 10), f)

f = open("cluster_graph_5_5.graph", "w+")
write_graph(cluster_graph(5, 5, 5), f)

f = open("complete10.graph", "w+")
write_graph(nx.complete_graph(10), f)

f = open("complete100.graph", "w+")
write_graph(nx.complete_graph(100), f)

f = open("complete1000.graph", "w+")
write_graph(nx.complete_graph(1000), f)

#f = open("complete10000.graph", "w+")
#write_graph(nx.complete_graph(10000), f)

#f = open("complete100000.graph", "w+")
#write_graph(nx.complete_graph(100000), f)

f = open("barbell4-4.graph", "w+")
write_graph(nx.barbell_graph(4, 0), f)


f = open("barbell10-10.graph", "w+")
write_graph(nx.barbell_graph(10, 1), f)

f = open("barbell25-25.graph", "w+")
write_graph(nx.barbell_graph(25, 1), f)

f = open("barbell100-100.graph", "w+")
write_graph(nx.barbell_graph(100, 1), f)

f = open("barbell1000-1000.graph", "w+")
write_graph(nx.barbell_graph(1000, 1), f)

#f = open("barbell10000-10000.graph", "w+")
#write_graph(nx.barbell_graph(10000, 10000), f)

#f = open("barbell100000-100000.graph", "w+")
#write_graph(nx.barbell_graph(100000, 100000), f)



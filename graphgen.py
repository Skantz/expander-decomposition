import random

import networkx as nx



def cluster_graph(n_clusters, size_clusters, n_crossing_edges=None):

    if not n_crossing_edges:
        n_crossing_edges = n_clusters - 1

    cluster_indices = [i for i in range(n_clusters)]

    g = nx.Graph(undirected=True)

    nodes_by_cluster = [[size_clusters*j + i for i in range(size_clusters)] for j in range(n_clusters)]
    for cluster in nodes_by_cluster:
        g.add_nodes_from(cluster)
    for cluster in nodes_by_cluster:
        for i in range(size_clusters):
            for j in range(i + 1, size_clusters):
                g.add_edge(cluster[i], cluster[j])

    #random.shuffle(cluster_indices)
    #1-indexed
    for i in range(n_clusters - 1):
        g.add_edge(cluster_indices[i]*size_clusters + 1, cluster_indices[i + 1]*size_clusters + 1)


    e_to_add = n_crossing_edges - n_clusters + 1
    while e_to_add > 0:
        n1, n2, c1, c2 = random.randint(0, size_clusters - 1), random.randint(0, size_clusters - 1), \
                         random.randint(0, n_clusters - 1), random.randint(0, n_clusters - 1)
        o1, o2 = random.randint(0, size_clusters - 1), random.randint(0, size_clusters - 1)
        if (c1 != c2 and not g.has_edge(n1 + c1*size_clusters + o1, n2 + c2*size_clusters + o2)):
            g.add_edge(n1 + c1*size_clusters + o1, n2 + c2*size_clusters + o2)
            e_to_add -= 1


    return g


def expander_graph(p_nodes): 
    g = nx.chordal_cycle_graph(p_nodes)
    return g

def write_graph(G, f):
    G = nx.convert_node_labels_to_integers(G, first_label=1)
    f.write(str(len(G.nodes())) +" "+ str(len(G.edges()))+"\n")
    for line in nx.generate_adjlist(G):
        f.write(line.partition(' ')[2])
        f.write("\n")

f = open("chordal_cycle_graph_2000.graph", "w+")
write_graph(expander_graph(2000), f)

f = open("cluster_graph_10_5.graph", "w+")
write_graph(cluster_graph(5, 200, 20), f)

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

f = open("barbell100-100.graph", "w+")
write_graph(nx.barbell_graph(100, 1), f)

f = open("barbell1000-1000.graph", "w+")
write_graph(nx.barbell_graph(1000, 1), f)

#f = open("barbell10000-10000.graph", "w+")
#write_graph(nx.barbell_graph(10000, 10000), f)

#f = open("barbell100000-100000.graph", "w+")
#write_graph(nx.barbell_graph(100000, 100000), f)



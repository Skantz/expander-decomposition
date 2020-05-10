import sys
import networkx as nx
import matplotlib.pyplot as plt

graph_file = sys.argv[1]
cut_file   = sys.argv[2]

G = nx.Graph()
with open(graph_file) as f_g:
    first_line = f_g.readline().split()
    n_nodes, n_edges = int(first_line[0]), int(first_line[1])
    iter = 0
    G.add_nodes_from([i for i in range(n_nodes)])
    for line in f_g:
        line_split = line.split()
        if len(line_split) <= 0:
            iter += 1
            continue
        base_node = int(line_split[0])
        for neighbor in line_split:
            G.add_edge(iter, int(neighbor) - 1)
        iter += 1


colors = [i for i in range(n_nodes)]
with open(cut_file) as f_c:
    c = 0
    cnt = 0
    for line in f_c:
        print(line.split())
        for n_str in line.split():
            colors[int(n_str)] = c
            cnt += 1
        c += 1
    assert(cnt == n_nodes)

print(colors)
print(len(colors))
print(len(G))

nx.draw_networkx(G, with_labels=False, node_color=colors, cmap=plt.get_cmap('Reds'), pos=nx.spring_layout(G))
plt.show()


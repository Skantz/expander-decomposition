import networkx as nx
import argparse

parser = argparse.ArgumentParser("Generate graphs")

parser.add_argument('--n_cliques', dest='n_cliques', type=int, nargs=1, help='number of cliques')
parser.add_argument('--clique_size', dest='clique_size', type=int, nargs=1, help='size of cliques')
parser.add_argument('--output_file', dest='output_file', type=str, nargs=1, help="output file relative path")

args = parser.parse_args()

n_cliques, clique_size, output_file = args.n_cliques[0], args.clique_size[0], args.output_file[0]

print(n_cliques, clique_size, output_file)

G = nx.Graph(directed=False)

G.add_nodes_from([i for i in range(n_cliques*clique_size)])
#G.add_edges_from([(i, i+1) for i in range(n_cliques*clique_size-1) if  (i+1) % clique_size != 0])
for c in range(0, n_cliques**2, n_cliques):
    G.add_edges_from([(c+i, c+j) for i in range(clique_size) for j in range(i + 1, clique_size)])
    G.add_edges_from([(clique_size * i, clique_size * (i + 1)) for i in range(n_cliques - 1)])

with open(output_file, 'w') as f:
    f.write(str(G.number_of_nodes()) + " " + str(G.number_of_edges()) + "\n")
    for n in range(G.number_of_nodes()):
        f.write(" ".join([str(v) for v in G[n]]) + "\n")

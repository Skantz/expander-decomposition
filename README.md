# expander-decomposition

Graph partitioning software based on https://arxiv.org/abs/1812.08958

Decomposes graphs with cluster expansion guarantees.


USE

install lemon

Compiles with
g++ main.cpp -std=c++17 -lemon -O3 -fopenmp

example use

./a.out  --G_phi=0.01 --H_phi=0.4 --vol=1 --h_ratio=0.01 -f random_3_regular_5000.graph

G_phi is our target conductance, H_phi a convergence threshold for cut-matching.

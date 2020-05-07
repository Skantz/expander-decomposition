# expander-decomposition
decompose graph with cluster expansion guarantee
https://arxiv.org/abs/1812.08958


install lemon

g++ main.cpp -std=c++17 -lemon -O3 -fopenmp

example use

./a.out  --G_phi=0.01 --H_phi=0.4 --vol=1 --h_ratio=0.01 -f random_3_regular_5000.graph --only_test_expander=true --known_phi=0.03 -S

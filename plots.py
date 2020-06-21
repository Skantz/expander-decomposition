
import numpy as np
from sklearn.preprocessing import normalize
from sklearn.metrics import pairwise_distances
"""
n = 3
matrix = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
matrix = np.array(matrix)

uniform = [ 1/n for _ in range(n)]
uniform = np.array(uniform)

walk    = [1] + [0 for _ in range(n - 1)]
walk    = np.array(walk)

matrix = normalize(matrix, axis=1, norm='l1')

distance = pairwise_distances((matrix @ walk).reshape(-1, 1), uniform.reshape(-1, 1))
distance = pairwise_distances((matrix @ (matrix @ (matrix @ walk))).reshape(-1, 1), uniform.reshape(-1, 1))
steps = 0
threshold = 0.01
while np.sum(pairwise_distances(walk.reshape(-1, 1), uniform.reshape(-1, 1))) > threshold:
    walk = matrix @ walk
    print(walk)
    print( np.sum(pairwise_distances(walk.reshape(-1, 1), uniform.reshape(-1, 1))))
    steps += 1

print("steps", steps)
"""

def test_convergence(matrix, threshold):
    #np.set_printoptions(threshold=np.inf)    
    n = len(matrix)
    assert(len(matrix[0]) == n)
    print("n:", n)
    #matrix = np.array(matrix)
    e = sum([sum(row) for row in matrix])
 
    print("e", e)

    colsums = [sum([matrix[i][j] for j in range(n)]) for i in range(n)]
    uniform = np.array([colsums[i]/e for i in range(n)])
    #uniform = np.array([1] + [0 for _ in range(n - 1)]) #np.array([sum([e for e in matrix[i]])/e for i in range(n)])
    assert(0.99 <=sum(uniform) <= 1.01)
    walk = np.array([1] + [0 for _ in range(n - 1)]) #Seed?
    #walk = np.array([1./n for _ in range(n)])

    sums = [sum(matrix[i]) for i in range(n)]
    for s in sums:
        assert(s > 0)
    matrix_with_drift = [[e/sums[j] for i, e in enumerate(subl)] for j, subl in enumerate(matrix)]

    assert(len(matrix_with_drift) == n)
    assert(len(matrix_with_drift[0]) == n)

    for i in range(n):
        matrix_with_drift[i][i] = sum(matrix_with_drift[i])
        assert(matrix_with_drift[i][i] > 0)
        cnt = sum(matrix_with_drift[i])
        matrix_with_drift[i] = [m/cnt for m in matrix_with_drift[i]]
        assert(0.99 <= sum(matrix_with_drift[i]) <= 1.01)

    matrix = matrix_with_drift
    print(matrix[:5])
    print("matrix to array")
    #matrix = np.array([[e if i !=j else np.count_nonzero(subl == 1) for i, e in enumerate(subl)] for j, subl in enumerate(matrix)])
    matrix = np.array(matrix)
    print("normalize")
    #matrix = normalize(matrix, axis=1, norm='l1')
    #print(matrix)
    steps = 0
    dist = np.sum(pairwise_distances(walk.reshape(-1, 1), uniform.reshape(-1, 1))) / n
    print("start iterate")
    print(np.sum(walk))
    while dist > threshold and steps < 20000:
        walk = np.matmul(matrix.T, walk.T)
        assert(0.99 < sum(walk) < 1.01)
        steps += 1
        dist = np.linalg.norm(walk - uniform, 1) # np.sum(pairwise_distances(walk.reshape(-1, 1) - uniform.reshape(-1, 1)) / n 

        if steps % 100 == 0:
            print(walk)
            #print(dist)

        print(steps)
        print("sum walk", np.sum(walk))
        print("dist", dist)
        ##print("uniform", uniform)
    print("sum uniform", np.sum(uniform))
    return steps


# In[2]:


import numpy as np

def graph_and_cut_to_numpy(gf, cf):

    graph = {}
    with open(gf, "r") as f:
        first_row = f.readline().strip("\n").split(" ")
        while [] in first_row:
            first_row.remove([])
        n, e = [int(e) for e in first_row]
        for i, row in enumerate(f.readlines()):
            graph[i] = []
            for v in row.strip("\n").split(" "):
                if v == "":
                    continue
                assert(int(v))
                graph[i].append(int(v))

    clusters = []
    with open(cf, "r") as f:
        for i, row in enumerate(f.readlines()):
            clusters.append([int(e) for e in row.strip("\n").split(" ") if e != ""])

    return graph, clusters

#graph, cuts = graph_and_cut_to_numpy("graphs/add20.graph", "old_results/add20.graphcut.txt")

graph, cuts = graph_and_cut_to_numpy("graphs/add32.graph", "old_results/add32.graphcut.txt")

subgs = []

for cut in cuts:

    v_set = set()
    for el in cut:
        v_set.add(int(el))

    #print(v_set)
    n = len(v_set)

    subg = {k:[] for k in cut}
    print("first element in cut", cut[0])
    for k in graph:
        for v in graph[k]:
            assert(int(v))
            if v in v_set and k in v_set:
                subg[k].append(v)


    print(" ".join([str((e, subg[e])) for e in subg])[:10])
    #rekey
    rekey_dict = {v:i for i, v in enumerate(subg)}
    rekey_reverse= {i:v for i, v in enumerate(subg)}

    new_subg = [[0 for _ in range(n)] for _ in range(n)]
    for k in subg: 
        assert(len(subg[k]) > 0)
        for v in subg[k]:
            new_subg[rekey_dict[k]][rekey_dict[v]] = 1
            new_subg[rekey_dict[v]][rekey_dict[k]] = 1    
    assert(len(new_subg) == n)
    assert(len(new_subg[0]) == n)
    test_convergence(new_subg, 0.1)

#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import glob
from os import path
import os
import subprocess
import sys

np.seterr(all='raise')


def metis_partition_file_converter(inp, out):
    import sys

    #inp = sys.argv[1]
    #out = sys.argv[2]

    partitions = {}

    def add(el, idx):
        if idx not in partitions:
            partitions[idx] = [el]
        else:
            partitions[idx].append(el)


    with open(inp, "r+") as f:
        el = 0
        for line in f:
            idx = int(line.strip())
            add(el, idx)
            el += 1

    with open(out, "w+") as f:
        idxs = list(partitions.keys())
        idxs.sort()
        assert(idxs[0] <= idxs[-1])
        for ix in idxs:
            f.write(" ".join([str(e) for e in partitions[ix]]) + "\n")




import numpy as np
from sklearn.preprocessing import normalize
from sklearn.metrics import pairwise_distances
import scipy.sparse as sparse

def test_convergence(matrix, threshold):
    #np.set_printoptions(threshold=np.inf)    
    #n = len(matrix)
    n = matrix.shape[0]
    #assert(len(matrix[0]) == n)

    #e = sum([sum(row) for row in matrix])
    e = matrix.sum()

    print("n;", n)
    print("e;", e)
    #matrix = np.array(matrix)

 
    #colsums = [sum([matrix[i][j] for j in range(n)]) for i in range(n)]
    colsums = matrix.sum(axis=1)
    #uniform = np.array([colsums[i]/e for i in range(n)], dtype=np.float64)
    uniform = np.zeros(shape=(n,))
    #print(colsums)
    for i in range(n):
        uniform[i] = colsums[i]/e
    
    #uniform = np.array([1] + [0 for _ in range(n - 1)]) #np.array([sum([e for e in matrix[i]])/e for i in range(n)])
    assert(0.99 <= uniform.sum() <= 1.01)

    #scipy matmul numpy : A.dot(v)

    walk = np.array([1.] + [0 for _ in range(n - 1)], np.float64) #Seed?
    #walk = np.array([1./n for _ in range(n)])

    #sums = [sum(matrix[i]) for i in range(n)]
    sums = matrix.sum(axis=0)
    assert((sums != 0).any())

    for s in sums:
        pass #assert(s > 0)
    
    #No longer sparse ?
    for i in range(n):
        for j in range(n):
            if matrix[i, j] != 0:
                matrix[i, j] = matrix[i, j] / int(sums[:, j])

    #matrix_with_drift = [[e/sums[j] for i, e in enumerate(subl)] for j, subl in enumerate(matrix)]

    #assert(len(matrix) == n)
    #assert(len(matrix[0]) == n)


    for i in range(n):
        matrix[i, i] = matrix[i, :].sum()
        assert(matrix[i, i] > 0)
        cnt = matrix[i, :].sum()
        matrix[i, :] /= cnt

        assert(0.99 <= matrix[i, :].sum() <= 1.01)

    #matrix = np.array(matrix, np.float64)
    matrix = sparse.coo_matrix(matrix)

    dist = np.linalg.norm(walk - uniform, 1)
    #dist  = np.sum(pairwise_distances(walk.reshape(-1, 1), uniform.reshape(-1, 1))) / n
    steps = 0
    while dist > threshold and steps < 10000: #and steps < len(matrix):
        walk = matrix.T.dot(walk.T)
        assert(0.99 < sum(walk)    < 1.01)
        assert(0.99 < sum(uniform) < 1.01)
        steps += 1
        dist = np.linalg.norm(walk - uniform, 1)
        
    """
    try:
        s_matrix     = sparse.coo_matrix(matrix)
        eigenvecs, _ = sparse.linalg.eigs(s_matrix, k=2)

        print("eigenvecs sample;",     eigenvecs[:5])
        print("first two eigvenvals;", eigenvecs[eigenvecs.argsort()[-2:][::-1]])
    except:
        print("Eigenvalue calculation failed")
    """

    print("dist;", dist)

    return steps



import numpy as np

def graph_and_cut_to_numpy(gf, cf):

    graph = {}
    with open(gf, "r") as f:
        first_row = f.readline().strip("\n").split(" ")
        while [] in first_row:
            first_row.remove([])
        n, e = [int(e) for e in first_row]
        for i, row in enumerate(f.readlines()):
            graph[i + 1] = []
            for v in row.strip("\n").split(" "):
                if v == "":
                    continue
                assert(int(v))
                graph[i + 1].append(int(v))
                try:
                    graph[v].append(int(i + 1))
                except KeyError:
                    graph[v] = [int(i + 1)]

    clusters = []
    #TOFIX: clusters are 0 index
    with open(cf, "r") as f:
        n_counter = 0
        e_counter = 0
        for i, row in enumerate(f.readlines()):
            clusters.append([int(e) + 1 for e in row.strip("\n").split(" ") if e != ""])
            #print(clusters[-1])
            n_counter += 1
        e_counter += len(clusters[-1])


    return graph, clusters


def test(graph_f, cut_f):
    graph, cuts = graph_and_cut_to_numpy(graph_f, cut_f)
    #graph, cuts = graph_and_cut_to_numpy("graphs/whitaker3.graph", "old_results/whitaker3.graph.part.94.chaco")
    subgs = []

    for clidx, cut in enumerate(cuts):

        if len(cut) == 1:
            print("Singleton, continue")
            continue
        v_set = set()
        for el in cut:
            v_set.add(int(el))

        #print(v_set)
        n = len(v_set)

        subg = {k:[] for k in cut}
        #print("first element in cut", cut[0])
        for k in graph:
            for v in graph[k]:
                if v in v_set and k in v_set:
                    subg[k].append(v)


        rekey_dict    = {v:i for i, v in enumerate(sorted(list(subg.keys())))}
        rekey_reverse = {i:v for i, v in enumerate(sorted(list(subg.keys())))}

        success = True
        #new_subg = [[0 for _ in range(n)] for _ in range(n)]
        new_subg = sparse.dok_matrix((n, n))

        for k in subg: 

            if len(subg[k]) == 0:
                print("WARNING - disconnected cluster;", clidx, "number of nodes in subgraph is;", len(cut))
                success = False
                break

            for v in subg[k]:
                #new_subg[rekey_dict[k]][rekey_dict[v]] = 1
                #new_subg[rekey_dict[v]][rekey_dict[k]] = 1
                new_subg[rekey_dict[k], rekey_dict[v]] = 1
                new_subg[rekey_dict[v], rekey_dict[k]] = 1

        if not success:
            continue


        steps = test_convergence(new_subg, 0.01)

        print("Cluster", clidx, "of size: ", n, "took n steps to converge;", steps)

        import math
        print("steps/ log2 nodes;", steps/math.log2(n))

        #Scipy and numpy can't calculate eigenvalues quick
        
    n = len(graph)
    if True: #n <= 5000
        trivial_subg = sparse.dok_matrix((n, n))
        for i, k in enumerate(graph):
            for j, _ in enumerate(graph[k]):
                #trivial_subg[i][j] = 1
                trivial_subg[i, j] = 1

        steps = test_convergence(trivial_subg, 0.01)
        print("Whole graph took n steps to converge:", steps)
        import math
        print("steps/ log2 nodes", steps/math.log2(n))
        


import glob


def run_decomp(graph_files):
    for graph in graph_files:
        print("decomp on;", graph)
        print(graph)

        #run decomposition on all graphs with time limit
        out_path = graph + ".out"
        #        get_ipython().system('time -p timeout 15m ./a.out  --G_phi=0.01 --H_phi=0.4 --vol=1 --h_ratio=0. -f $graph | tee "$out_path"')
        #get_ipython().system('time -p timeout 15m ./a.out  --G_phi=0.01 --H_phi=0.4 --vol=1 --h_ratio=0. -f $graph')
        os.system('time -p timeout 15m ./a.out  -S --G_phi=0.01 --H_phi=0.4 --vol=1 --h_ratio=0. -f ' + graph + " | tee " + out_path)

        #TODO return code
        if not glob.glob(graph + "cut.txt"):
            print("decomp on", graph, "did not finish in time")
            continue

        #how many clusters, k, did we get?
        with open(out_path) as f:
            lns = f.readlines()
            cluster_line = next(l for l in lns if "n clusters" in l)
            cluster_line.strip("\n")
            n_clusters = int(cluster_line.split(";")[1])

        if n_clusters <= 1:
            print("No cut found, continue")
            continue

        print("n clusters found", n_clusters)
        #run metis on k
        metis_stdio_path = graph + ".out.metis"
        #TOFIX: why does this freeze?
        #get_ipython().system('timeout 10s /usr/bin/gpmetis -ufactor=1000 $graph $n_clusters -contig >  $metis_stdio_path')
        #get_ipython().system('timeout 10s /usr/bin/gpmetis -ufactor=1000 $graph $n_clusters -contig')
        #subprocess.run('timeout 5s /usr/bin/gpmetis -ufactor=1000 ' + graph + ' ' + str(n_clusters) + '-contig', shell=True) #stdout=subprocess.std, stderr=subprocess.PIPE)
        #subprocess.Popen('timeout 5s /usr/bin/gpmetis -ufactor=1000 ' + graph + ' ' + str(n_clusters) + '-contig', shell=True)#.communicate()[0] 
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

        #metis_file = graph + ".out.metis"
        #Metis does not always return k parts
        #metis_file = graph + ".part." + str(n_clusters)
        #metis_file = glob.glob(graph + ".part.*")[0]
        decomp_file = graph + "cut.txt"
        
        print("Converting metis file")
        #metis_partition_file_converter(metis_file, metis_file)
        print("Conversion done")

        #rw_graph      = graph + ".row_whole"
        #rw_file_ours  = graph + ".rw_ours"
        #rw_file_metis = graph + ".rw_metis"

        #TOFIX This should be saving to file!
        print("decomp ours")
        test(graph, decomp_file)
        #print("decomp metis")
        #test(graph, metis_file)        

        bname = os.path.basename(graph).split(".")[0]
        #get_ipython().system('mkdir -p results/"$bname"')
        os.system('mkdir -p results/"$bname"')
        for f in glob.glob("".join(graph.split(".")[:-1]) + "*"):
            if f != graph:
                #get_ipython().system('mv $f results/"$bname"/')
                os.system('mv $f results/"$bname"/')



import glob


#g1s = list(glob.glob("synthetic/*"))
#run_decomp(g1s)
blacklist = ["graphs/144.graph"]
g2s = list(glob.glob("graphs/*"))
g2s.sort()
g2s = [s for s in g2s if s not in blacklist]
run_decomp(g2s)


# In[ ]:





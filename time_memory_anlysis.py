import sys
sys.path.append('clust')
import tree
import time
import random
import numpy as np
from tqdm import tqdm

num2dna = ['A', 'T', 'G', 'C']

def channel(txSeqs, pi, pd, ps):
    rxSeqs = []
    simerror = [0, 0, 0]
    for _, tx in enumerate(txSeqs):
        i, rx = 0, ''
        while i < len(tx):
            if randsel(pi): # Insertion
                rx += num2dna[np.random.randint(0, 4)]
                simerror[0] += 1
                continue      
            if not randsel(pd): # Deletion
                if randsel(ps): # Substitution
                    tmp = num2dna.copy()
                    tmp.remove(tx[i])
                    rx += tmp[np.random.randint(0, 3)]
                    simerror[2] += 1
                else:
                    rx += tx[i]
            else:
              simerror[1] += 1
            i += 1
        rxSeqs.append(rx)
    return rxSeqs, simerror

def randsel(p):
    sel = np.random.rand() < p
    return sel

def gen_random_dataset_with_size(data_size, seq_len, pi, pd, ps, cluster_size = [50, 200]):
    data = []
    count = 0
    while True:
        repeat = random.randint(cluster_size[0], cluster_size[1])
        tx = ''
        for _ in range(seq_len):
            tx += num2dna[random.randint(0, 3)]
        rxs, _ = channel([tx for _ in range(repeat)], pi, pd, ps)
        data.extend(rxs)
        count += len(rxs)
        if count > data_size:
            break
    data = data[:data_size]
    return data

def cluster_without_return(data, tree_depth):
    trie = tree.new_tree(tree_depth)
    clust_ind = 0
    for seq in tqdm(data):
        seq = seq[:tree_depth]
        if len(seq) == tree_depth:
            result = tree.quick_search(trie, seq, 3, 2)
            if result.label == 0:
                clust_ind += 1
                tree.insert(trie, seq, clust_ind)

if __name__=='__main__':
    print("Time Test...")
    for i in range(5, 8):
        data_size = 10**i
        data_len = 10
        print('Number of reads:%d...'%data_size)
        print('Generate data...')
        data = gen_random_dataset_with_size(data_size, data_len, 0.05, 0.05, 0.05)
        begin = time.time()
        cluster_without_return(data, data_len)
        end = time.time()
        print('Time %d'%(end-begin))

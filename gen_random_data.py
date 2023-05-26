import random
import numpy as np
from clust import tree
num2dna = ['A', 'T', 'G', 'C']
dna2num = {
    'A': '00',
    'T': '01',
    'G': '10',
    'C': '11',
    }
pattern = ['AAA', 'TTT', 'GGG', 'CCC']

def channel(txSeqs, pi, pd, ps):
    rxSeqs = []
    for _, tx in enumerate(txSeqs):
        i, rx = 0, ''
        while i < len(tx):
            if randsel(pi): # Insertion
                rx += num2dna[np.random.randint(0, 4)]
                continue      
            if not randsel(pd): # Deletion
                if randsel(ps): # Substitution
                    tmp = num2dna.copy()
                    tmp.remove(tx[i])
                    rx += tmp[np.random.randint(0, 3)]
                else:
                    rx += tx[i]
            i += 1
        rxSeqs.append(rx)
    return rxSeqs

def randsel(p):
    sel = np.random.rand() < p
    return sel

def lfsr(state, mask):
    result = state
    nbits = mask.bit_length() - 0b001
    while True:
        result = result << 0b001
        xor = result >> nbits
        if xor != 0:
            result ^= mask
        yield result

def screen_homopolymers(data, pattern):
    for w in pattern:
        if w in data:
            return 0
    return 1 # pass

def screen_gc(data):
    gc = (data.count('G') + data.count('C') + 0.0) / len(data)
    if gc < 0.4 or gc > 0.6:
        return 0
    return 1 # pass

def screen_edit(index, trie, trie_th):
    if trie_th == 0:
        align = tree.quick_search(trie, index, trie_th, 0)
    else:
        align = tree.search(trie, index, trie_th)
    if align.label > 0:
        return 0
    else:
        return 1 # pass

def gen_random_data(n_cluster, n_repeat, extra, pi, pd, ps, dist=0):
    state = 0b001
    mask = 0b100000000000000000000000011000101
    psnr = lfsr(state, mask)
    data = []
    clInds = []
    trie = tree.new_tree(16)
    random.seed(255)
    np.random.seed(255)
    for i in range(n_cluster):
        seed = '{:032b}'.format(next(psnr))
        payload = ''
        for j in range(extra):
            payload += num2dna[random.randint(0, 3)]

        index = ''
        for j in range(16):
            index += num2dna[int(seed[j * 2:j * 2 + 2], 2)]
        invalid = not(screen_homopolymers(index, pattern) and screen_gc(index) and screen_edit(index, trie, dist))
        while invalid:
            seed = '{:032b}'.format(next(psnr))
            index = ''
            for j in range(16):
                index += num2dna[int(seed[j * 2:j * 2 + 2], 2)]
            invalid = not(screen_homopolymers(index, pattern) and screen_gc(index) and screen_edit(index, trie, dist))
        tree.insert(trie, index, i+1)
        tx = index + payload
        rxs = channel([tx for _ in range(n_repeat)], pi, pd, ps)
        data.extend(rxs)
        clInds.extend([i + 1 for _ in range(n_repeat)])
        print(i)
    return data, clInds

if __name__ == '__main__':

    data_file = 'testdata/randomReads_10000_dist_0.txt'
    clust_file = 'testdata/randomTags_10000_dist_0.txt'
    data, tags = gen_random_data(10000, 20, 10, 0.01, 0.01, 0.01, dist=0)

    # shuffle indexes
    inds = [i for i in range(len(data))]
    np.random.seed(255)
    np.random.shuffle(inds)

    with open(data_file, 'w') as f:
        f.truncate(0)
        for i, ind in enumerate(inds):
            f.write('%d %s\n'%(i+1, data[ind]))

    with open(clust_file, 'w') as f:
        f.truncate(0)
        for i, ind in enumerate(inds):
            f.write('%d,%d\n'%(i+1, tags[ind]))
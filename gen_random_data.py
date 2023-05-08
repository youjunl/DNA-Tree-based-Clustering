import random
import numpy as np

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

def gen_random_data(n_cluster, n_repeat, extra, pi, pd, ps):
    state = 0b001
    mask = 0b100000000000000000000000011000101
    psnr = lfsr(state, mask)
    data = []
    clInds = []
    dna_dict = {}
    for i in range(n_cluster):
        seed = '{:032b}'.format(next(psnr))
        tx = ''
        for j in range(16):
            tx += num2dna[int(seed[j * 2:j * 2 + 2], 2)]
        for j in range(extra):
            tx += num2dna[random.randint(0, 3)]

        while not (screen_homopolymers(tx, pattern) and screen_gc(tx) and tx[:16] not in dna_dict.keys()):
            seed = '{:032b}'.format(next(psnr))
            tx = ''
            for j in range(16):
                tx += num2dna[int(seed[j * 2:j * 2 + 2], 2)]
            for j in range(extra):
                tx += num2dna[random.randint(0, 3)]
        dna_dict[tx[:16]] = True
        rxs = channel([tx for _ in range(n_repeat)], pi, pd, ps)
        data.extend(rxs)
        clInds.extend([i + 1 for _ in range(n_repeat)])
    print(len(dna_dict.keys()))
    return data, clInds

if __name__ == '__main__':

    data_file = 'testdata/randomReads_10000.txt'
    clust_file = 'testdata/randomTags_10000.txt'
    data, tags = gen_random_data(10000, 20, 10, 0.01, 0.01, 0.01)

    # shuffle indexes
    inds = [i for i in range(len(data))]
    np.random.shuffle(inds)

    with open(data_file, 'w') as f:
        f.truncate(0)
        for i, ind in enumerate(inds):
            f.write('%d %s\n'%(i+1, data[ind]))

    with open(clust_file, 'w') as f:
        f.truncate(0)
        for i, ind in enumerate(inds):
            f.write('%d,%d\n'%(i+1, tags[ind]))
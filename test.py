import random
from reedsolo import RSCodec
import numpy as np
import random
import matplotlib.pyplot as plt
from clust.seq_kmeans import SeqKmeans, SoftSeqKmeans
from clust.tree import Trie
import time
from collections import defaultdict
# Target
def lfsr(state, mask):
  result = state
  nbits = mask.bit_length()-1
  while True:
    result = (result << 1)
    xor = result >> nbits
    if xor != 0:
        result ^= mask
    yield result

def screen_homopolymers(data, pattern):
    for w in pattern:
        if w in data:
            return False
    return True

def channel(txSeqs, pi, pd, ps):
    rxSeqs = []
    simerror = [0, 0, 0]
    for tx in txSeqs:
        i, rx = 0, ''
        while i < len(tx):
            if randsel(pi): # Insertion
                rx += num2dna[random.randint(0, 3)]
                simerror[0] += 1
                continue      
            if not randsel(pd): # Deletion
                if randsel(ps): # Substitution
                    tmp = num2dna.copy()
                    tmp.remove(tx[i])
                    rx += tmp[random.randint(0, 2)]
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

def dna_to_int_array(dna_str):
    #convert a string like ACTCA to an array of ints like [10, 2, 4]
    s = ''
    for ch in dna_str:
        s += '{0:02b}'.format(int(dna_to_num[ch]))
    
    data = [int(s[t:t+8],2) for t in range(0,len(s), 8)]
    return data

def int_array_to_dna(data):
    s = ''
    bin_data = ''
    for num in data:
        bin_data += '{0:08b}'.format(num)
    for i in range(0, len(bin_data), 2):
        s += num_to_dna[int(bin_data[i:i+2], 2)]
    return s

def clust(tree, dnaData):
    test_num = 0
    tree_threshold = 6
    h_drift = 6
    indexList = []
    dna_number = 0
    timer = time.time()
    error = np.array([0, 0, 0])
    core_set = dict()
    for i, seq in enumerate(dnaData):
        test_num += 1
        align = tree.fuzz_fin(seq, tree_threshold)
        if sum(align[1]) < h_drift:
            indexList.append((clInds[i], align[0]))
            error = error + np.array(align[1])
        else:
            dna_number += 1
            tree.insert(seq, dna_number)
            core_set[dna_number] = seq
            indexList.append((clInds[i], dna_number)) # Fix Clover's problem

    # Output Result
    t = (time.time() - timer)
    print('time spent:', t)
    sortedList = sorted(indexList, key = lambda k: k[0])
    print('Estimated channel')
    print(simerror)
    print(error)
    prob = error / len(dnaData) / len(dnaData[0])
    print(prob)
    with open('result.txt', 'w') as f:
        f.truncate(0)
        for t in sortedList:
            f.write(str(t[0]) + ' ' + str(t[1]) + '\n')
    return indexList

def computeAccuracy(indexList):
    f = open('refInd.txt', 'r').readlines()
    refind = [int(t.strip()) for t in f]
    algind = indexList
    # caculate freqs
    clustNum = defaultdict(int)
    results = indexList
    for ind in refind:
        clustNum[ind] += 1

    results = sorted(results, key=lambda k: k[1])
    nClust = len(clustNum.keys())

    maxClustNum = results[-1][1] if len(results) > 0 else 0
    clusters = [[] for _ in range(maxClustNum)]
    for pair in results:
        clusters[pair[1]-1].append(pair[0])

    clusters = [c for c in clusters if c != []]
  
    score = [0] * max(clustNum.keys())
    for cluster in clusters:
        # Check if all the tags in a clusters are the same.
        tags = set(cluster)
        if len(tags) > 1:
            # This cluster is invalid
            continue
        # Maximize the size of the clusters.
        tag = int(cluster[0])
        score[tag-1] = max(score[tag-1], len(cluster))

    # Compute accuracy
    gamma = [i*0.05 for i in range(0, 21)]
    acc = [0] * len(gamma)
  
    for i, g in enumerate(gamma):
        cnt = 0
        for tag in clustNum.keys():
            if score[tag-1] / clustNum[tag] >= g:
                cnt += 1
        acc[i] = cnt / nClust
    return gamma, acc

def dna_to_seed(dna_str):
    # revert seed from the DNA sequence
    s = ''
    for ch in dna_str:
        s += '{0:02b}'.format(int(dna_to_num[ch]))    
    return int(s, 2)

def seed_to_dna(seed):
    bin_str = '{:032b}'.format(seed)
    dna_str = ''
    for t in range(0, len(bin_str), 2):
        dna_str += num_to_dna[int(bin_str[t:t+2],2)]
    return dna_str

def SSC(indexList, dnaData):
    begin = time.time()
    maxClusterIndex = max(indexList, key=lambda k: k[1])
    tags = [index[1] for index in indexList]
    clusters = [[] for _ in range(maxClusterIndex[1])]
    for i, read in enumerate(dnaData):
        clustInd, tag = indexList[i]
        clusters[tag-1].append(read)

    seeds = []
    seedMap = dict()
    clusterMap = dict()
    for i, cluster in enumerate(clusters):
        # Skip empty clusters
        if len(cluster) == 0:
            continue
        
        # Record the frequencies of seeds
        freq = dict()
        for read in cluster:
            seed = dna_to_seed(read[:16])
            if seed not in freq.keys():
                freq[seed] = 1
            else:
                freq[seed] += 1
        
        # Majority selection for the center seed
        maxSeed = -1
        maxFreq = 0
        for k in freq.keys():
            if freq[k] > maxFreq:
                maxFreq = freq[k]
                maxSeed = k
        if maxSeed == -1:
            continue
        
        # Create a mapping
        clusterInd = i + 1
        if maxSeed not in seeds:
            seeds.append(maxSeed)
            seedMap[maxSeed] = clusterInd # For merging the other clusters
            clusterMap[clusterInd] = clusterInd
        else:
            clusterMap[clusterInd] = seedMap[maxSeed]
            
    for tag in set(tags):
        if tag not in clusterMap.keys():
            clusterMap[tag] = tag # Map to itself
        
    newIndexList = []
    for clustInd, tag in indexList:
        newIndexList.append((clustInd, clusterMap[tag]))
    tags = [index[1] for index in newIndexList]
    print('time spent:', time.time()-begin)
    print('Found %d clusters'%len(clusters))
    print('Generate %d clusters'%len(set(tags)))
    return newIndexList

def SSCRS(indexList, dnaData):
    begin = time.time()
    maxClusterIndex = max(indexList, key=lambda k: k[1])
    tags = [index[1] for index in indexList]
    clusters = [[] for _ in range(maxClusterIndex[1])]
    for i, read in enumerate(dnaData):
        clustInd, tag = indexList[i]
        clusters[tag-1].append(read)
    cnt_corrected = 0
    cnt_detected = 0
    seeds = []
    seedMap = dict()
    clusterMap = dict()
    for i, cluster in enumerate(clusters):
        # Skip empty clusters
        if len(cluster) == 0:
            continue
        
        # Record the frequencies of seeds
        freq = dict()
        for read in cluster:
            if 'N' in read:
                continue
            dna = dna_to_int_array(read)
            detected = False
            try:
                #First is the decoded/repaired message
                #Second is the decoded message and error correction code (because both get repaired in reality)
                #Third is the position of the errors and erasures.
                data_corrected, _, _ = codec.decode(dna)
                detected = True
                
            except:
                detected = False #could not correct the code

            if detected:
                #we will encode the data again to evaluate the correctness of the decoding
                data_again = list(codec.encode(data_corrected)) #list is to convert byte array to int
                if np.count_nonzero(dna != data_again) > max_hamming: #measuring hamming distance between raw input and expected raw input
                    #too many errors to correct in decoding
                    cnt_detected += 1                 
                    seed = dna_to_seed(read[:16])
                    continue
                
                else:
                    cnt_corrected += 1
                    dna_str_corrected = int_array_to_dna(data_again)
                    seed = dna_to_seed(read[:16])
                

            else:
                seed = dna_to_seed(read[:16])
                continue

            if seed not in freq.keys():
                freq[seed] = 1
            else:
                freq[seed] += 1
        
        # Majority selection for the center seed
        maxSeed = -1
        maxFreq = 0
        for k in freq.keys():
            if freq[k] > maxFreq:
                maxFreq = freq[k]
                maxSeed = k
        if maxSeed == -1:
            continue
        
        # Create a mapping
        clusterInd = i + 1
        if maxSeed not in seeds:
            seeds.append(maxSeed)
            seedMap[maxSeed] = clusterInd # For merging the other clusters
            clusterMap[clusterInd] = clusterInd
        else:
            clusterMap[clusterInd] = seedMap[maxSeed]
            
    for tag in set(tags):
        if tag not in clusterMap.keys():
            clusterMap[tag] = tag # Map to itself

    newIndexList = []
    for clustInd, tag in indexList:
        newIndexList.append((clustInd, clusterMap[tag]))
    tags = [index[1] for index in newIndexList]
    print('time spent:', time.time()-begin)
    print('Found %d clusters'%len(clusters))
    print('Generate %d clusters'%len(set(tags)))
    print('Reads corrected: %d'%cnt_corrected)
    return newIndexList

def SSCML(indexList, dnaData, iter):
    begin = time.time()
    maxClusterIndex = max(indexList, key=lambda k: k[1])
    tags = [index[1] for index in indexList]
    clusters = [[] for _ in range(maxClusterIndex[1])]
    for i, read in enumerate(dnaData):
        clustInd, tag = indexList[i]
        clusters[tag-1].append(read)
    cnt_corrected = 0
    cnt_detected = 0
    seeds = []
    seedMap = dict()
    clusterMap = dict()
    for i, cluster in enumerate(clusters):
        # Skip empty clusters
        if len(cluster) == 0:
            continue
        
        # Record the frequencies of seeds
        freq = dict()
        for read in cluster:
            seed = dna_to_seed(read[:16])
            if seed not in freq.keys():
                freq[seed] = 1
            else:
                freq[seed] += 1
        
        # Majority selection for the center seed
        maxSeed = -1
        maxFreq = 0
        for k in freq.keys():
            if freq[k] > maxFreq:
                maxFreq = freq[k]
                maxSeed = k
        if maxSeed == -1:
            # Select random
            maxSeed = cluster[random.randint(0, len(cluster)-1)][:16]
        
        # Create a mapping
        clusterInd = i + 1
        if maxSeed not in seeds:
            seeds.append(maxSeed)
            seedMap[maxSeed] = clusterInd # For merging the other clusters
            clusterMap[clusterInd] = clusterInd
        else:
            clusterMap[clusterInd] = seedMap[maxSeed]
            
    for tag in set(tags):
        if tag not in clusterMap.keys():
            clusterMap[tag] = tag # Map to itself

    # Clustering with Kmeans
    alphabet = np.array(['A', 'T', 'G', 'C'])
    alg = SeqKmeans(n_centroid = nCluster, centroid_length = 16, alphabet=alphabet)
    data = []
    dataClusterInd = []
    for seed in seeds:
        dnaSeed = seed_to_dna(seed)
        data.append(dnaSeed)
        dataClusterInd.append(seedMap[seed])
    
    lcurve = alg.fit(np.array(data), n_iter=iter)
    labels = alg.transform(data)

    print('Number of clusters before KMeans: %d'%len(data))
    print('Number of clusters after KMeans: %d'%len(set(labels)))

    # Create mapping
    MLMap = dict()
    for i, label in enumerate(labels):
        if label not in MLMap.keys():
            MLMap[label] = dataClusterInd[i] # Set a root for merging

        else:
            # Update merging
            clusterMap[dataClusterInd[i]] = MLMap[label]

    newIndexList = []
    for clustInd, tag in indexList:
        newIndexList.append((clustInd, clusterMap[tag]))

    tags = [index[1] for index in newIndexList]
    print('time spent:', time.time()-begin)
    print('Found %d clusters'%len(clusters))
    print('Generate %d clusters'%len(set(tags)))
    return newIndexList

if __name__ == '__main__':
    state = 0b001
    mask = 0b100000000000000000000000011000101
    psnr = lfsr(state, mask)
    pi, pd, ps = 0.001, 0.001, 0.001
    nCluster = 20
    clusterSize = [15, 40]
    extraLen = 16
    rscode = 2
    max_hamming = 2
    codec = RSCodec(rscode)
    num2dna = ['A', 'T', 'G', 'C']
    dna_to_num = {'A':0,'C':1,'G':2,'T':3}
    num_to_dna = {0:'A',1:'C',2:'G',3:'T'}
    pattern = ['AAA', 'TTT', 'GGG', 'CCC']
    psnr = lfsr(state, mask) # PSNR Initialization
    clInds = []
    simerror = np.array([0, 0, 0])
    with open('seeds.txt', 'w') as f:
        f.truncate(0)
        for i in range(nCluster):
            nRepeat = random.randint(clusterSize[0], clusterSize[1])
            txSeed = '{:032b}'.format(next(psnr))
            tx = ''
            for j in range(16):
                tx += num2dna[int(txSeed[j*2:j*2+2], 2)]
            for j in range(extraLen):
                tx += num2dna[random.randint(0, 3)]
            if rscode > 0:
                txInt = dna_to_int_array(tx)
                txCode = codec.encode(txInt)
                tx = int_array_to_dna(txCode)
            while not screen_homopolymers(tx, pattern):
                txSeed = '{:032b}'.format(next(psnr))
                tx = ''
                for j in range(16):
                    tx += num2dna[int(txSeed[j*2:j*2+2], 2)]
                for j in range(extraLen):
                    tx += num2dna[random.randint(0, 3)]
            rxs, e = channel([tx for _ in range(nRepeat)], pi, pd, ps)
            simerror = simerror + np.array(e)
            clInds.extend([i+1 for _ in range(nRepeat)]) 
            for rx in rxs:
                rx = rx[0:16+extraLen]
                f.write(rx + '\n')

    with open('clustInd.txt', 'w') as f:
        f.truncate(0)
        for ind in clInds:
            f.write(str(ind) + '\n')

    # Data preprocessing
    data = np.genfromtxt('seeds.txt' ,delimiter='\n',dtype=str)
    clInds = np.genfromtxt('clustInd.txt' ,delimiter='\n',dtype='uint64')
    nRead = len(data)
    print('Import %d PRNs'%nRead)

    # Random Shuffling
    inds = np.arange(nRead)
    Shuffling = False
    if Shuffling:
        np.random.shuffle(inds)
        data = data[inds]
        clInds = clInds[inds]
    dnaData = [s for s in data]

    with open('refInd.txt', 'w') as f:
        f.truncate(0)
        for ind in clInds:
            f.write(str(ind) + '\n')

    tree = Trie()
    indexList = clust(tree, dnaData)
    gamma, acc = computeAccuracy(indexList)

    indexListSSC = SSC(indexList, dnaData)
    gamma, accSSC = computeAccuracy(indexListSSC)

    indexListSSCRS = SSCRS(indexList, dnaData)
    gamma, accSSCRS = computeAccuracy(indexListSSCRS)

    indexListSSCML = SSCML(indexList, dnaData, iter=20)
    gamma, accSSCML = computeAccuracy(indexListSSCML)

    plt.figure()
    plt.plot(gamma, acc, '-')
    plt.plot(gamma, accSSC, '--^')
    plt.plot(gamma, accSSCRS, '--o')
    plt.plot(gamma, accSSCML, '--*')
    plt.legend(['original', 'ssc', 'sscrs', 'sscml'])
    plt.show()
from clust import tree as tr
from clust import chainer_lcs
dna_to_num = {'A':0,'C':1,'G':2,'T':3}
num_to_dna = {0:'A',1:'C',2:'G',3:'T'}

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

def ssc(indexList, dnaData, read_len):
    maxClusterIndex = max(indexList, key=lambda k: k[1])
    tags = [index[1] for index in indexList]
    clusters = [[] for _ in range(maxClusterIndex[1])]
    for i, pair in enumerate(dnaData):
        _, read = pair
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
            seed = dna_to_seed(read[:read_len])
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
    return newIndexList

def ssc_repeat(indexList, dnaData, read_len, tree_threshold, filter=False, spliter=False):
    maxClusterIndex = max(indexList, key=lambda k: k[1])
    tags = [index[1] for index in indexList]
    clusters = [[] for _ in range(maxClusterIndex[1])]
    for i, pair in enumerate(dnaData):
        _, read = pair
        clustInd, tag = indexList[i]
        clusters[tag-1].append(read)
    seedMap = dict()
    clusterMap = dict()
    readMap = dict()
    tree = tr.Trie()
    for i, cluster in enumerate(clusters):
        # Skip empty clusters
        if len(cluster) == 0:
            continue
        # Record the frequencies of seeds
        freq = dict()
        for read in cluster:
            seed = read[:read_len]
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

        if spliter:
            # inp = [r[:read_len] for r in cluster]
            # out = chainer_lcs.mutual_lcs(inp, read_len)
            pass

        # Create a mapping
        clusterInd = i + 1

        # Clustering with Tree Structure
        seedStr = maxSeed
        align = tree.fuzz_fin(seedStr[:read_len], tree_threshold)
        if sum(align[1]) < tree_threshold:
            if filter:
                # Prefiltering with LCS
                # print('%s Merge To %s Err: %d'%(seedStr[-tree_depth:], readMap[align[0]], sum(align[1])))
                pass
            clusterMap[clusterInd] = align[0] # Merging
        else:
            tree.insert(seedStr[:read_len], clusterInd)
            seedMap[maxSeed] = clusterInd # For merging the other clusters
            clusterMap[clusterInd] = clusterInd
            # if filter:
            #     readMap[clusterInd] = seedStr[:read_len]
            
    for tag in set(tags):
        if tag not in clusterMap.keys():
            clusterMap[tag] = tag # Map to itself

    newIndexList = []
    for clustInd, tag in indexList:
        newIndexList.append((clustInd, clusterMap[tag]))
    tags = [index[1] for index in newIndexList]
    return newIndexList
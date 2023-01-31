import random
from tqdm import tqdm
from multiprocessing.dummy import Pool as ThreadPool
dna_alphabet = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
ind_alphabet = {'A': 0, 'T': 1, 'G': 2, 'C': 3}

params = {
    'r': 20, # Threshold for edit distance
    'q': 3, # Substring's length
    'w': 4, # Anchor's length
    'l': 8, # Hash's length
    'n': 152, # Real segmentation length
    'theta_low': 20,
    'theta_high': 40,
    'local_step': 10,
    'comm_step': 20,
    'core_num' : 1,
    'debug' : False,
}

class Sequence:
    def __init__(self, data, tag, cluster) -> None:
        self.data = data
        self.cluster = cluster
        self.hash = ''
        self.tag = tag
        self.iv = []

def computeIndicator(str, blockLen, q):
    blockIVLen = 4**q
    strLen = len(str)
    iv = [0] * int(strLen/blockLen+1)*blockIVLen
    for i in range(strLen - q):
        index = 0
        carry = 1
        for j in range(q):
            index += carry * ind_alphabet[str[i+j]]
            carry *= 4
        iv[index+int(i/blockLen)*blockIVLen] = 1
    return iv

def editDistance(str1, str2, m, n):
    edit_m = [[0 for _ in range(n+1)] for _ in range(m+1)] # m rows and n cols
    for i in range(m+1):
        edit_m[i][0] = i
    for i in range(n+1):
        edit_m[0][i] = i
    for i in range(m):
        for j in range(n):
            if str2[j] == str1[i]:
                edit_m[i+1][j+1] = edit_m[i][j]
            else:
                edit_m[i+1][j+1] = min([edit_m[i][j+1]+1,
                                        edit_m[i+1][j]+1,
                                        edit_m[i][j]+1, ])
    return edit_m[-1][-1]

def compute_hash(str, anchor, hashLen):
    if(hashLen>len(str)):
        return str
    pos = str.find(anchor)
    # not matched, return the last pos
    if pos < 0:
        return str[-hashLen:]
    
    return str[pos:pos+hashLen]

def compute_bsd(iv1, iv2):
    # compute hamming distance between indicator vectors
    ans = 0
    if len(iv1) < len(iv2):
        iv1, iv2 = iv2, iv1
    for i in range(len(iv2)):
        ans += iv1[i] ^ iv2[i] # xor
    for i in range(len(iv2), len(iv1)):
        ans += iv1[i]
    return ans

def random_anchor(anchorLen):
    anchor = ''
    for _ in range(anchorLen):
        anchor += dna_alphabet[random.randint(0, 3)]
    return anchor

def compute_comm(inData):
    S = []
    blockLen = 25
    print('Preprocessing...')
    def preprocess(seq):
        seq.iv = computeIndicator(seq.data, blockLen, params['q'])
        S.append([seq])

    with ThreadPool() as pool:
        pool.map(preprocess, inData)
        pool.close()
        pool.join()
        
    core = params['core_num']
    hashLen = params['w'] + params['l']
    print('Clustering')
    for _ in tqdm(range(params['comm_step'])):
        #print('comm step %d'%i)
        # initialize a random anchor
        #anchor = random_anchor(params['w'])
        # Partition map
        C = [[] for _ in range(core)]
        for cluster in S:
            C[random.randint(0, core-1)].append(cluster)
        with ThreadPool(params['core_num']) as pool:
            pool.map(compute_local, C)
            pool.close()
            pool.join()
        S = []
        for i in range(core):
            S.extend(C[i])
        S = [c for c in S if len(c) > 0]
        # Assign cluster number
        cnt = 1
        for cluster in S:
            for seq in cluster:
                seq.cluster = cnt
            cnt += 1
    return S

def compute_local(Clusters):
    hashLen = params['w'] + params['l']
    blockLen = 22
    numBlock = params['n'] / blockLen

    for i in range(params['local_step']):
        anchor = random_anchor(params['w'])
        bucket = []
        indicators = []
        merge_map = {}
        # Get representatives and compute hash value
        for j in range(len(Clusters)):
            if len(Clusters[j]) < 1:
                continue
            sample = Clusters[j][random.randint(0, len(Clusters[j])-1)]
            sample.hash = compute_hash(sample.data, anchor, hashLen)
            bucket.append(sample)
            merge_map[sample.cluster] = j

        # Put in bucket and merge clusters and sort the reads based on hashval;
        bucket = sorted(bucket, key=lambda k: k.hash)

        # Compute indicator for samples in buckets

        for sample in bucket:
            indicators.append(sample.iv)

        # Merge
        bucketLen = len(bucket)
        for first in range(bucketLen-1):
            if len(Clusters[merge_map[bucket[first].cluster]]) == 0:
                # is empty, skip
                continue

            # Only compare adjacent elements in sorted bucket
            for second in range(first+1, min(first+20, bucketLen)):
                if len(Clusters[merge_map[bucket[second].cluster]]) == 0:
                    # is empty, skip
                    continue

                ind1 = bucket[first].cluster
                ind2 = bucket[second].cluster
                if ind1 == ind2:
                    continue

                str1 = bucket[first].hash
                str2 = bucket[second].hash
                distance = compute_bsd(indicators[first], indicators[second])
                if distance > params['theta_high']:
                    continue
                elif distance >= params['theta_low']:
                    # compare edit distance
                    if editDistance(str1, str2, len(str1), len(str2)) < params['r']:
                        #  merge
                        
                        for t in range(len(Clusters[merge_map[ind2]])):
                            Clusters[merge_map[ind2]][t].cluster = ind1
                        Clusters[merge_map[ind1]].extend(Clusters[merge_map[ind2]])
                        Clusters[merge_map[ind2]].clear()

                elif distance < params['theta_low']:
                    #  merge
                    for t in range(len(Clusters[merge_map[ind2]])):
                        Clusters[merge_map[ind2]][t].cluster = ind1
                    Clusters[merge_map[ind1]].extend(Clusters[merge_map[ind2]])
                    Clusters[merge_map[ind2]].clear()
    Clusters = [c for c in Clusters if c != []]
        
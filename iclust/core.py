import random
import multiprocessing as mp
from tqdm import tqdm
num_to_dna = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
dna_to_num = {'A': 0, 'T': 1, 'G': 2, 'C': 3}

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

def random_anchor(anchorLen):
    anchor = ''
    for _ in range(anchorLen):
        anchor += num_to_dna[random.randint(0, 3)]
    return anchor

class Sequence:
    def __init__(self, data, tag, cluster) -> None:
        self.data = data
        self.cluster = cluster
        self.hash = ''
        self.tag = tag
        self.iv = []

class Process():
    def __init__(self, params) -> None:
        self.params = params

    def computeIndicator(self, data: Sequence):
        blockLen = self.params['block_len']
        q = self.params['q']
        seq = data.data
        blockIVLen = 4**q
        strLen = len(seq)
        iv = [0] * int(strLen/blockLen+1)*blockIVLen
        for i in range(strLen - q):
            index = 0
            carry = 1
            for j in range(q):
                index += carry * dna_to_num[seq[i+j]]
                carry *= 4
            iv[index+int(i/blockLen)*blockIVLen] = 1
        # out = np.zeros(blockIVLen, dtype=np.uint8)
        # for offset in range(0, len(iv), blockIVLen):
        #     out |= iv[offset:offset+blockIVLen]
        data.iv = int(''.join(map(str, iv)), 2)
        return data

    def compute_comm(self, inData):
        params = self.params
        blockLen = params['block_len']
        print('Preprocessing...')
        # with mp.Pool(processes= params['thread']) as pool:
        #     for result in tqdm(pool.imap_unordered(self.computeIndicator, inData), total=len(inData)):
        #         S.append([result])
        S = [[self.computeIndicator(data)] for data in tqdm(inData)]
        core = params['core_num']
        print('Clustering...')
        for _ in tqdm(range(params['comm_step'])):
            C = [[] for _ in range(core)]
            for cluster in S:
                C[random.randint(0, core-1)].append(cluster)
            for in_clusters in C:
                self.compute_local(in_clusters)
            S = []
            for out_clusters in C:
                S.extend(out_clusters)
            S = [clusters for clusters in S if len(clusters) > 0]
            # Assign cluster number
            cnt = 1
            for cluster in S:
                for seq in cluster:
                    seq.cluster = cnt
                cnt += 1
        return S

    def compute_local(self, in_clusters):
        params = self.params
        hashLen = params['w'] + params['l']
        for i in range(params['local_step']):
            anchor = random_anchor(params['w'])
            bucket = []
            merge_map = {}
            # Get representatives and compute hash value
            for j in range(len(in_clusters)):
                if len(in_clusters[j]) < 1:
                    continue
                sample = in_clusters[j][random.randint(0, len(in_clusters[j])-1)]
                sample.hash = compute_hash(sample.data, anchor, hashLen)
                bucket.append(sample)
                merge_map[sample.cluster] = j

            # Put in bucket and merge in_clusters and sort the reads based on hashval;
            bucket = sorted(bucket, key=lambda k: k.hash)

            # Merge
            bucketLen = len(bucket)
            for first in range(bucketLen-1):
                if len(in_clusters[merge_map[bucket[first].cluster]]) == 0:
                    # is empty, skip
                    continue

                for second in range(first+1, min(first+20+1, bucketLen)):
                    if len(in_clusters[merge_map[bucket[second].cluster]]) == 0:
                        # is empty, skip
                        continue

                    ind1 = bucket[first].cluster
                    ind2 = bucket[second].cluster
                    if ind1 == ind2:
                        continue

                    str1 = bucket[first].hash
                    str2 = bucket[second].hash
                    # comput hamming distance
                    xor_bit = bucket[first].iv ^ bucket[second].iv
                    distance = bin(xor_bit).count('1')
                    if distance > params['theta_high']:
                        continue
                    elif distance >= params['theta_low']:
                        # compare edit distance
                        if editDistance(str1, str2, len(str1), len(str2)) < params['r']:
                            #  merge
                            for t in range(len(in_clusters[merge_map[ind2]])):
                                in_clusters[merge_map[ind2]][t].cluster = ind1
                            in_clusters[merge_map[ind1]].extend(in_clusters[merge_map[ind2]])
                            in_clusters[merge_map[ind2]].clear()

                    elif distance < params['theta_low']:
                        #  merge
                        for t in range(len(in_clusters[merge_map[ind2]])):
                            in_clusters[merge_map[ind2]][t].cluster = ind1
                        in_clusters[merge_map[ind1]].extend(in_clusters[merge_map[ind2]])
                        in_clusters[merge_map[ind2]].clear()
        in_clusters = [c for c in in_clusters if c != []]
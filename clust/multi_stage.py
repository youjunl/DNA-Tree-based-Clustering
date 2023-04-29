from clust import tree
from tqdm import tqdm
from collections import Counter
import random
dna_to_num = {'A':0,'C':1,'G':2,'T':3}
num_to_dna = {0:'A',1:'C',2:'G',3:'T'}

def random_anchor(anchorLen):
    anchor = ''
    for _ in range(anchorLen):
        anchor += num_to_dna[random.randint(0, 3)]
    return anchor

def compute_hash(str, anchor, hashLen):
    if(hashLen>len(str)):
        return -1
    pos = str.find(anchor)
    # not matched, return -1
    if pos < 0 or pos+len(anchor)+hashLen>len(str):
        return -1
    
    return str[pos+len(anchor):pos+len(anchor)+hashLen]

def iter_clust(indexList, dnaData, anchor_len, hash_len, tree_threshold):
    maxClusterIndex = max(indexList, key=lambda k: k[1])[1]
    clust_nums = [index[1] for index in indexList]
    clust_idxes = [[] for _ in range(len(clust_nums))]
    clusters = [[] for _ in range(maxClusterIndex)]
    for i, pair in enumerate(dnaData):
        _, read = pair
        dna_tag, clust_num = indexList[i]
        clusters[clust_num-1].append(read) 
        clust_idxes[clust_num-1].append(dna_tag) # record the order of data
    nextClusterIndex = maxClusterIndex + 1

    newIndexList = []
    # # Spliter
    # for i, cluster in enumerate(tqdm(clusters)):
    #     # Skip empty clusters
    #     if len(cluster) == 0:
    #         continue

    #     anchor = random_anchor(anchor_len)
    #     hash_vals = []
    #     ctree = tree.new_tree(hash_len)
    #     rep = cluster[random.randint(0, len(cluster)-1)]
    #     rep_hash = compute_hash(rep, anchor, hash_len)
        
    #     if rep_hash == -1:
    #         for index in clust_idxes[i]:
    #             newIndexList.append((index, i+1))
    #         continue

    #     # Initialize the tree with representative
    #     tree.insert(ctree, rep_hash, i+1) # i+1 is current tag

    #     # compute hash
    #     for seq in cluster:
    #         hash_val = compute_hash(seq, anchor, hash_len)
    #         if hash_val != -1:
    #             hash_vals.append(compute_hash(seq, anchor, hash_len))
    #         else:
    #             hash_val = seq[-hash_len:]
    #             hash_vals.append(hash_val)

    #     for j, hash_val in enumerate(hash_vals):
    #         align = tree.search(ctree, hash_val, tree_threshold)
    #         if align.label > 0:
    #             newIndexList.append((clust_idxes[i][j], align.label))
                    
    #         else:
    #             tree.insert(ctree, hash_val, nextClusterIndex)
    #             newIndexList.append((clust_idxes[i][j], nextClusterIndex))
    #             nextClusterIndex += 1

    # Merger with anchor
    ctree = tree.new_tree(hash_len)
    anchor = random_anchor(anchor_len)
    clusterMap = dict()
    for i, cluster in enumerate(tqdm(clusters)):
        # Skip empty clusters
        if len(cluster) == 0:
            continue

        rep = cluster[random.randint(0, len(cluster)-1)]
        rep_hash = compute_hash(rep, anchor, hash_len)

        if rep_hash == -1:
            continue

        # Create a mapping
        clusterInd = i + 1

        # Clustering with Suffix
        if len(cluster) > 5:
            tree.insert(ctree, rep_hash, clusterInd)
            clusterMap[clusterInd] = clusterInd
        else:
            align = tree.search(ctree, rep_hash, tree_threshold)
            if align.label > 0:
                clusterMap[clusterInd] = align.label # Merging
            else:
                tree.insert(ctree, rep_hash, clusterInd)
                clusterMap[clusterInd] = clusterInd
            
    for clust_num in set(clust_nums):
        if clust_num not in clusterMap.keys():
            clusterMap[clust_num] = clust_num # Map to itself

    newIndexList = []
    for clustInd, tag in indexList:
        newIndexList.append((clustInd, clusterMap[tag]))

    return newIndexList

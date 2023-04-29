import matplotlib.pyplot as plt
from collections import defaultdict
from tqdm import tqdm
import getopt, sys
import cupy
# inputReads = 'testdata/toClustReads.txt'
# inputTags = 'testdata/toClustTags.txt'
# algResult = 'output_ERR1816980_clust.txt'

# inputReads = 'testdata/toClustSmallReads.txt'
# inputTags = 'testdata/toClustSmallTags.txt'
# algResult = 'output_file_Small_1.txt'

inputReads = 'testdata/toClustRead_p10.txt'
inputTags = 'testdata/toClustTag_p10.txt'
algResult = 'output_file_P10_clust.txt'

def editDistance(A, B):
    N, M = len(A), len(B)
    # Create an array of size NxM
    dp = [[0 for i in range(M + 1)] for j in range(N + 1)]

    # Base Case: When N = 0
    for j in range(M + 1):
        dp[0][j] = j
    # Base Case: When M = 0
    for i in range(N + 1):
        dp[i][0] = i
    # Transitions
    for i in range(1, N + 1):
        for j in range(1, M + 1):
            if A[i - 1] == B[j - 1]:
                dp[i][j] = dp[i-1][j-1]
            else:
                dp[i][j] = 1 + min(
                    dp[i-1][j], # Insertion
                    dp[i][j-1], # Deletion
                    dp[i-1][j-1] # Replacement
                )

    return dp[N][M]

def compute(fileIn, tags, clustNum, gamma):
    # Total number of clusters in the labeled data
    nClust = 0
    for k in clustNum.keys():
        if clustNum[k] > 1:
            nClust += 1

    # Read clustering results
    results = []

    with open(fileIn, 'r') as f:
        for i, text in enumerate(f.readlines()):
            # The second element of each line is the output clustering index of the algorithms.
            ind, cluster = map(int, text.strip().split(','))
            # Save these pairs in a list.
            if tags[ind-1] != -1:
                results.append((tags[ind-1], cluster, ind))
            
    reads = []
    with open(inputReads, 'r') as f:
        for i, text in enumerate(f.readlines()):
            # The second element of each line is the output clustering index of the algorithms.
            ind, read = text.strip().split(' ')
            # Save these pairs in a list.
            reads.append(read) 

    # Sort the result according to the clustering index and get the max index number of clusters that algorithms output.
    sorted_results = sorted(results, key=lambda k: k[1])
    maxClustNum = sorted_results[-1][1] if len(sorted_results) > 0 else 0

    # Cluster the tags according to the clustering indexes.
    clusters = [[] for _ in range(maxClustNum)]
    for i in range(len(results)):
        clusters[results[i][1]-1].append((results[i][0], reads[i], results[i][2]))

    ref_cluster_map = defaultdict(list)
    for i in range(len(tags)):
        ref_cluster_map[tags[i]].append(i+1)

    # Remove empty clusters
    clusters = [c for c in clusters if c != []]
    # clusters = [c for c in clusters if len(c)>1]
    # A list that stores the maximum size of correct clustering for different tags.
    score = [0] * max(clustNum.keys())
    for cluster in clusters:
        # Check if all the tags in a clusters are the same.
        tags = [p[0] for p in cluster]
        tags_key = set(tags)
        if len(tags_key) > 1:
            # This cluster is invalid
            print('Wrong clusters...')
            print(tags_key)
            for p in cluster:
                print(p[1][:21], p[0])
            continue
        # Maximize the size of the clusters.
        tag = cluster[0][0]


        # if len(cluster) > 15:
        #     inds = [p[2] for p in cluster]
        #     if len(cluster) != len(ref_cluster_map[tag]):
        #         print(cluster[0][1][:21])
        #         print('Not clusted...')
        #         for ind in ref_cluster_map[tag]:
        #             if ind not in inds:
        #                 dist = editDistance(cluster[0][1][:21], reads[ind-1][:21])
        #                 print('Seq:%s\t Ind:%d\t Now:%s\t Dist:%d'%(reads[ind-1][:21], ind, results[ind-1][1],dist))
        #     print(inds)

        score[tag-1] = max(score[tag-1], len(cluster))

    # Compute accuracy under different gammas.
    acc = [0 for _ in gamma]
    for i, g in enumerate(gamma):
        cnt = 0
        for tag in clustNum.keys():
            if clustNum[tag] > 1 and score[tag-1] / clustNum[tag] >= g:
                cnt += 1
        acc[i] = cnt / nClust
    print('Num cluster: %d. Inp cluster: %d'%(len(clustNum.keys()), len(clusters)))
    return acc

if __name__ == '__main__':

    labeled = inputTags
    indexes = algResult

    # Set range of gamma
    step_num = 20
    gamma = [i/step_num for i in range(0, step_num+1)]
    acc = [0] * len(gamma)

    # Count frequencies of tags in the labeled data.
    print('Counting tags in the labeled data...')
    clustNum = defaultdict(int)
    lines = open(labeled, 'r').readlines()
    tags = [-1 for _ in range(int(lines[-1].split(',')[0]))]
    for text in lines:
        ind, tag = map(int, text.strip().split(','))
        clustNum[tag] += 1
        tags[ind-1] = tag
        
    # Compute accuracy for each input clustering indexes file
    print('Computing accuracy...')
    acc = compute(indexes, tags, clustNum, gamma)
    print(acc)

    

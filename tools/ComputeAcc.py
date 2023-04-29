import matplotlib.pyplot as plt
from collections import defaultdict
from tqdm import tqdm
import getopt, sys

text = '''
#############################################################################################
#Implementation of accuracy in paper:                                                       #
#Rashtchian, Cyrus, et al. "Clustering billions of reads for DNA data storage." NIPS 2017.  #
#############################################################################################
'''

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
                results.append((tags[ind-1], cluster))
            

    # Sort the result according to the clustering index and get the max index number of clusters that algorithms output.
    results = sorted(results, key=lambda k: k[1])
    maxClustNum = results[-1][1] if len(results) > 0 else 0

    # Cluster the tags according to the clustering indexes.
    clusters = [[] for _ in range(maxClustNum)]
    for i in range(len(results)):
        clusters[results[i][1]-1].append(results[i][0])

    # Remove empty clusters
    clusters = [c for c in clusters if c != []]
    # clusters = [c for c in clusters if len(c)>1]
    # A list that stores the maximum size of correct clustering for different tags.
    score = [0] * max(clustNum.keys())
    for cluster in clusters:
        # Check if all the tags in a clusters are the same.
        tags = set(cluster)
        if len(tags) > 1:
            # This cluster is invalid
            continue
        
        # Maximize the size of the clusters.
        tag = cluster[0]
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
    print(text)
    # Inputs Management
    helpInfo = '[Usage]\nComputeAcc.py <labeled data> <cluster indexes file> <cluster indexes file 2> ... <output file>'
    try:
        opts, args = getopt.getopt(sys.argv[1:],"h",[])
        for opt, arg in opts:
            print(opt,arg)
            if opt in ['-h']:
                print(helpInfo)
                sys.exit()
    
    except getopt.GetoptError:
        print('Error!')
        print(helpInfo)
        sys.exit(2)

    if len(args)<3:
        print('No enough inputs, labeled data, cluster indexes and output filename are expected')
        print(helpInfo)
        sys.exit(2)

    labeled = args[0]
    indexes = args[1:-1]
    outfile = args[-1]

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
    outAcc = [[] for _ in indexes]
    for i, infile in enumerate(indexes):
        acc = compute(infile, tags, clustNum, gamma)
        outAcc[i] = acc

    with open(outfile, 'w') as f:
        f.truncate(0)
        
        # header
        f.write('Gamma,')
        for g in gamma:
            f.write('%f,'%g)
        f.write('\n')

        # results
        for i, accData in enumerate(outAcc):
            f.write('%s,'%indexes[i])
            for acc in accData:
                f.write('%.4f,'%acc)
            f.write('\n')
    
    print('Finished!')

import matplotlib.pyplot as plt
from collections import defaultdict
from tqdm import tqdm
import numpy as np
import getopt, sys

text = '''
#############################################################################################
#Cluster size comparison                                                                    #
#############################################################################################
'''
def plot_hist(ax, cnt_dist, label):
    max_num = max(cnt_dist.keys())
    x = []
    for i in range(1, max_num+1):
        x.append(cnt_dist[i])
    bins = list(range(1, max_num+1))
    ax.hist(x, bins, histtype="stepfilled", alpha=0.6, density=True, label=label)
    
def compute(fileIn, tags, clustNum, gamma):
    # Total number of clusters in the labeled data
    nClust = len(clustNum.keys())

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
    
    # A list that stores the maximum size of correct clustering for different tags.
    new_clustCnt = defaultdict(int)
    for cluster in clusters:
        # Check if all the tags in a clusters are the same.
        tags = set(cluster)
        if len(tags) > 1:
            # This cluster is invalid
            continue

        new_clustCnt[len(cluster)] += 1

    return new_clustCnt

if __name__ == '__main__':
    print(text)
    # Inputs Management
    helpInfo = '[Usage]\nComputeNum.py <labeled data> <cluster indexes file> <cluster indexes file 2> ... <output file>'
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
    gamma = [i*0.1 for i in range(0, 11)]
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
        
    clustCnt = defaultdict(int)
    for k in clustNum.keys():
        clustCnt[clustNum[k]] += 1
    # Compute accuracy for each input clustering indexes file

    fig, ax = plt.subplots()
    plot_hist(ax, clustCnt, 'Ground Truth')
    for i, infile in enumerate(indexes):
        new_clustCnt = compute(infile, tags, clustNum, gamma)
        plot_hist(ax, new_clustCnt, infile)

    plt.legend()
    plt.show()

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
    ax.bar(bins, x, alpha=0.6, label=label)
    
def compute(fileIn):
    # A list that stores the maximum size of correct clustering for different tags.
    clust_size = defaultdict(int)
    with open(fileIn, 'r') as f:
        for line in f.readlines():
            _, tag, _ = line.strip().split('\t')
            tag = int(tag)
            clust_size[tag] += 1
    clust_num_dist = defaultdict(int)
    for k in clust_size.keys():
        clust_num_dist[clust_size[k]] += 1
    return clust_num_dist

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

    # Count frequencies of tags in the labeled data.
    print('Counting tags in the labeled data...')
    clustNum = defaultdict(int)
    lines = open(labeled, 'r').readlines()
    tags = [-1 for _ in range(int(lines[-1].split(',')[0]))]
    for text in lines:
        ind, tag = map(int, text.strip().split(','))
        clustNum[tag] += 1
        tags[ind-1] = tag
        
    clust_num_dist = defaultdict(int)
    for k in clustNum.keys():
        clust_num_dist[clustNum[k]] += 1
    # Compute accuracy for each input clustering indexes file

    with open(outfile, 'w') as fout:
        fout.truncate(0)
        fout.write(',')
        for i in range(1, 200 + 1):
            fout.write('%d,'%i)
        fout.write('labeled,')
        max_num = max(clust_num_dist.keys())
        for i in range(1, max_num + 1):
            fout.write('%d,'%clust_num_dist[i])
        fout.write('\n')

        for i, infile in enumerate(indexes):
            clust_num_dist = compute(infile)
            max_num = max(clust_num_dist.keys())
            fout.write('%s,'%infile)
            for i in range(1, max_num + 1):
                fout.write('%d,'%clust_num_dist[i])
            fout.write('\n')

import matplotlib.pyplot as plt
from collections import defaultdict
from tqdm import tqdm
import getopt, sys

text = '''
#############################################################################################
#Implementation of accuracy (or called purity) in these papers:                             #
#Clover: tree structure-based efficient DNA clustering for DNA-based data storage           #
#MeShClust: an intelligent tool for clustering DNA sequences                                #
#############################################################################################
'''

def compute(fileIn, tags):
    # Read clustering results
    results = []
    with open(fileIn, 'r') as f:
        for i, text in enumerate(f.readlines()):
            text = text.strip()
            content = text.split(',')
            # The second element of each line is the output clustering index of the algorithms.
            cluster = int(content[1])
            # Save these pairs in a list.
            results.append((tags[i], cluster))

    # Sort the result according to the clustering index and get the max index number of clusters that algorithms output.
    results = sorted(results, key=lambda k: k[1])
    maxClustNum = results[-1][1] if len(results) > 0 else 0

    # Cluster the tags according to the clustering indexes.
    clusters = [[] for _ in range(maxClustNum)]
    for i in range(len(results)):
        clusters[results[i][1]-1].append(results[i][0])
    clusters = [c for c in clusters if c != []]
    
    # Count the maximum frequencies in every clasters
    cnt = 0
    for cluster in clusters:
        # A dictionary for storing the number of tags
        tagNums = defaultdict(int)
        tags = set(cluster)
        for tag in cluster:
            tagNums[tag] += 1
        
        maxNum = 0
        for tag in tags:
            if tagNums[tag] > maxNum:
                maxNum = tagNums[tag]
        cnt += maxNum
    acc = cnt / len(results)
    return acc

if __name__ == '__main__':
    print(text)
    # Inputs Management
    helpInfo = '[Usage]\nComputePur.py <labeled data> <cluster indexes file> <luster indexes file 2> ... <output file>'
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

    tags = []
    with open(labeled, 'r') as f:
        for text in tqdm(f.readlines()):
            ind, tag = map(int, text.strip().split(','))
            tags.append(tag)

    # Compute accuracy for each input clustering indexes file
    print('Computing accuracy...')
    outAcc = [0 for _ in indexes]
    for i, infile in enumerate(indexes):
        acc = compute(infile, tags)
        outAcc[i] = acc

    with open(outfile, 'w') as f:
        f.truncate(0)
        for i, acc in enumerate(outAcc):
            f.write(indexes[i])
            f.write('\n')
            f.write('%.4f, '%acc)
            f.write('\n')
    
    print('Finished!')
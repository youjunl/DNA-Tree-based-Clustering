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

def compute(fileIn, clustNum, gamma):
    # Total number of clusters in the labeled data
    nClust = len(clustNum.keys())

    # Read clustering results
    results = []
    with open(fileIn, 'r') as f:
        for text in f.readlines():
            text = text.strip()
            content = text.split(',')
            # The first element of each line, which is the original tag in the labeled data.
            tag = int(content[0])
            # The second element of each line, which is the output clustering index of the algorithms.
            cluster = int(content[1])
            # Save these pairs in a list.
            results.append((tag, cluster))

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
            if score[tag-1] / clustNum[tag] >= g:
                cnt += 1
        acc[i] = cnt / nClust
    return acc

if __name__ == '__main__':
    print(text)
    # Inputs Management
    helpInfo = '[Usage]\nComputeAcc.py <labeled data> <cluster indexes file> <luster indexes file 2> ... <output file>'
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

    if len(args<3):
        print('No enough inputs, labeled data, cluster indexes and output filename are expected')
        print(helpInfo)
        sys.exit(2)

    labeled = args[0]
    indexes = args[1:-1]
    outfile = args[-1]

    # Set range of gamma
    gamma = [i*0.02 for i in range(20, 51)]
    acc = [0] * len(gamma)

    # Count frequencies of tags in the labeled data.
    print('Counting tags in the labeled data...')
    clustNum = defaultdict(int)
    with open(labeled, 'r') as f:
        for text in tqdm(f.readlines()):
            origin = int(text.split(' ')[0])
            clustNum[origin] += 1

    # Compute accuracy for each input clustering indexes file
    print('Computing accuracy...')
    outAcc = [0 for _ in indexes]
    for i, infile in enumerate(indexes):
        acc = compute(infile, clustNum, gamma)
        outAcc[i] = acc

    with open(outfile, 'w') as f:
        f.truncate(0)
        f.write('Gamma:\n')
        for g in gamma:
            f.write('%.4f, '%g)
        f.write('\n')
        for acc in outAcc:
            f.write(indexes[i])
            f.write('\n')
            for g in acc:
                f.write('%.4f, '%g)
            f.write('\n')
    
    print('Finished!')
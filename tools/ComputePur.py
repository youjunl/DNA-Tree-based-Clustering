from sklearn.metrics.cluster import normalized_mutual_info_score
from collections import defaultdict
from tqdm import tqdm
import getopt, sys

text = '''
#############################################################################################
#Implementation of identity (or called purity) in these papers:                             #
#Clover: tree structure-based efficient DNA clustering for DNA-based data storage           #
#MeShClust: an intelligent tool for clustering DNA sequences                                #
#############################################################################################
'''

def compute(fileIn, tags):
    # Read clustering results
    results = []
    # Number of alignment results
    truth = []
    pred = []
    with open(fileIn, 'r') as f:
        for i, text in enumerate(f.readlines()):
            # The second element of each line is the output clustering index of the algorithms.
            ind, cluster_ind = map(int, text.strip().split(','))
            # Save these pairs in a list.
            if ind > len(tags):
                continue
            if tags[ind-1] != 0:
                truth_label = tags[ind-1]
                results.append((truth_label, cluster_ind))
                truth.append(truth_label)
                pred.append(cluster_ind)

    # Calculate NMI
    nmi = normalized_mutual_info_score(truth, pred)
    print(fileIn)
    print('NMI: %f'%nmi)
    # Sort the result according to the clustering index and get the max index number of clusters that algorithms output.
    results = sorted(results, key=lambda k: k[1])
    maxClustNum = results[-1][1] if len(results) > 0 else 0

    # Cluster the tags according to the clustering indexes.
    clusters = [[] for _ in range(maxClustNum)]
    for i in range(len(results)):
        clusters[results[i][1]-1].append(results[i][0])
    clusters = [c for c in clusters if c != []]
    # Remove the cluster that size = 1
    clusters = [c for c in clusters if len(c) > 1]
    # Count the maximum frequencies in every clasters
    cnt = 0
    seq_cnt = 0
    for cluster in clusters:
        seq_cnt += len(cluster)
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
    acc = cnt / seq_cnt
    print('Purity: %f'%acc)
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

    lines = open(labeled, 'r').readlines()
    tags = [0 for _ in range(int(lines[-1].split(',')[0]))]
    for text in tqdm(lines):
        ind, tag = map(int, text.strip().split(','))
        tags[ind-1] = tag

    # Compute purity for each input clustering indexes file
    print('Computing purity...')
    outAcc = [0 for _ in indexes]
    for i, infile in enumerate(indexes):
        acc = compute(infile, tags)
        outAcc[i] = acc

    with open(outfile, 'w') as f:
        f.truncate(0)
        for i, acc in enumerate(outAcc):
            f.write('%s,'%indexes[i])
            f.write('%f\n'%acc)
    
    print('Finished!')

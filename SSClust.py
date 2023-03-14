import numpy as np
from reedsolo import RSCodec
from collections import defaultdict
from tqdm import tqdm
import getopt, sys

text = '''
#############################################################################################
#Implementation of second stage clustering:                                                 #
#############################################################################################
'''

dna_to_num = {'A':0,'C':1,'G':2,'T':3}
num_to_dna = {0:'A',1:'C',2:'G',3:'T'}
max_hamming = 2
codec = RSCodec(2) # 2Bytes

def dna_to_int_array(dna_str):
    #convert a string like ACTCA to an array of ints like [10, 2, 4]
    s = ''
    for ch in dna_str:
        s += '{0:02b}'.format(int(dna_to_num[ch]))
    
    data = [int(s[t:t+8],2) for t in range(0,len(s), 8)]
    return data

def byte_to_bin(s):
    #convert byte data (\x01 \x02) to DNA data: ACTC
    bin_data = ''
    for b in s:
        bin_data += '{0:08b}'.format(b)
    return bin_data

def bin_to_dna(bin_str):
    dna_str = ''
    for t in range(0, len(bin_str), 2):
        dna_str += num_to_dna[int(bin_str[t:t+2],2)]
    return dna_str

def int_array_to_dna(data):
    s = ''
    bin_data = ''
    for num in data:
        bin_data += '{0:08b}'.format(num)
    for i in range(0, len(bin_data), 2):
        s += num_to_dna[int(bin_data[i:i+2], 2)]
    return s
    

def dna_to_seed(dna_str):
    # revert seed from the DNA sequence
    s = ''
    for ch in dna_str:
        s += '{0:02b}'.format(int(dna_to_num[ch]))    
    return int(s, 2)

if __name__ == '__main__':
    print(text)
    # Inputs Management
    helpInfo = '[Usage]\SSClustering.py <reads data> <cluster indexes file> <output file>'
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

    reads = args[0] # Reads
    tags = args[1] # Cluster indexes
    outfile = args[2]
    
    # Processing original reads
    lines = open(reads, 'r').readlines()
    reads = ['' for _ in range(len(lines))]
    for i, text in enumerate(lines):
        tag, read = text.strip().split(' ')
        reads[i] = read

    lines = open(tags, 'r').readlines()
    tags = [0 for _ in range(len(lines))]
    inds = [0 for _ in range(len(lines))]
    for i, text in enumerate(lines):
        ind, tag = text.strip().split(',')
        tags[i] = int(tag)
        inds[i] = int(ind)
        
    maxClusterIndex = max(tags)
    numCluster = len(set(tags))
    print('Detect number of cluster: %d'%numCluster)
    
    clusters = [[] for _ in range(maxClusterIndex)]
    for i, read in enumerate(reads):
        tag = tags[i]
        clusters[tag-1].append(read)

    # Remove empty clusters
    # clusters = [c for c in clusters if len(c) != 0]
    
    cnt_corrected = 0
    cnt_detected = 0
    seeds = []
    seedMap = dict()
    clusterMap = dict()
    for i, cluster in enumerate(clusters):
        # Skip empty clusters
        if len(cluster) == 0:
            continue
        
        # Record the frequencies of seeds
        freq = dict()
        for read in cluster:
            if 'N' in read:
                continue

            seed = dna_to_seed(read[:16])
            if seed not in freq.keys():
                freq[seed] = 1
            else:
                freq[seed] += 1
        
        # Majority selection for the center seed
        maxSeed = -1
        maxFreq = 0
        for k in freq.keys():
            if freq[k] > maxFreq:
                maxFreq = freq[k]
                maxSeed = k
        if maxSeed == -1:
            continue
        
        # Create a mapping
        clusterInd = i + 1
        if maxSeed not in seeds:
            seeds.append(maxSeed)
            seedMap[maxSeed] = clusterInd # For merging the other clusters
            clusterMap[clusterInd] = clusterInd
        else:
            clusterMap[clusterInd] = seedMap[maxSeed]
            
    for tag in tags:
        if tag not in clusterMap.keys():
            clusterMap[tag] = tag # Map to itself
        
    print('Found %d clusters'%len(clusters))
    print('Generate %d clusters'%len(seeds))
    
    # Process original clustering result
    f = open(outfile, 'w')
    f.truncate(0)
    for i in range(len(tags)):
        f.write('%d,%d\n'%(i+1, clusterMap[tags[i]]))
    f.close()

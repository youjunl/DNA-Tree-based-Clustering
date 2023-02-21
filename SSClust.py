import numpy as np
from reedsolo import RSCodec
from collections import defaultdict
from tqdm import tqdm
import getopt, sys

text = '''
#############################################################################################
#Implementation of second stage clustering (for DNA fountain):                              #
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

if __name__ == '__main__':
    print(text)
    # Inputs Management
    helpInfo = '[Usage]\nSSClust.py <reads data> <cluster indexes file> <output file>'
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
    indexes = args[1]
    outfile = args[2]
    
    # Processing original reads
    lines = open(labeled, 'r').readlines()
    reads = ['' for _ in range(len(lines))]
    for text in lines:
        ind, read = text.strip().split(' ')
        reads[int(ind)-1] = read

    # Processing clustering
    results = []
    lines = open(indexes, 'r').readlines()
    tags = [[] for _ in range(len(lines))]
    for text in lines:
        ind, cluster = map(int, text.strip().split(','))
        results.append((ind, cluster))

    # Sort the result according to the clustering index and get the max index number of clusters that algorithms output.
    results = sorted(results, key=lambda k: k[1])
    maxClustNum = results[-1][1] if len(results) > 0 else 0   

    # Cluster the reads according to the clustering indexes.
    clusters = [[] for _ in range(maxClustNum)]
    for i in range(len(results)):
        clusters[results[i][1]-1].append((results[i][0], reads[results[i][0]-1]))

    # Remove empty clusters
    clusters = [c for c in clusters if c != []]
    
    cnt_corrected = 0
    cnt_detected = 0
    
    for cluster in clusters:
        for pair in cluster:
            ind, read = pair
            if 'N' in read:
                continue
            dna = dna_to_int_array(read)
            try:
                #First is the decoded/repaired message
                #Second is the decoded message and error correction code (because both get repaired in reality)
                #Third is the position of the errors and erasures.
                data_corrected, _, _ = codec.decode(dna)

            except:
                continue #could not correct the code
            
            #we will encode the data again to evaluate the correctness of the decoding
            data_again = list(codec.encode(data_corrected)) #list is to convert byte array to int
            if np.count_nonzero(dna != data_again) > max_hamming: #measuring hamming distance between raw input and expected raw input
                #too many errors to correct in decoding
                cnt_detected += 1                 
                continue
            
            cnt_corrected += 1
            dna_corrected = bin_to_dna(byte_to_bin(data_corrected))

            # Get seed from the sample and compare with the tree

    
    
    print('%d out of %d are corrected'%(cnt_corrected, len(lines)))
    print('%d out of %d are detected'%(cnt_detected, len(lines)))
    print('Finished!')

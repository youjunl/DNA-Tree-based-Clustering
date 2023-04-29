"""Main program



"""
import time
import math
import random
from multiprocessing import Process, Manager, Queue, Pool
from tqdm import tqdm
from clust import load_config as lc
from clust.multi_stage import *
from clust import tree
from collections import defaultdict, Counter
import numpy as np
import mmap

def mapcount(filename):
    f = open(filename, "r+")
    buf = mmap.mmap(f.fileno(), 0)
    lines = 0
    readline = buf.readline
    while readline():
        lines += 1
    return lines

def one_hot_encoding(seq, dict_alphabet, max_seq_length, smooth=1.):
    out = np.zeros((len(dict_alphabet), max_seq_length), dtype=np.float32)
    if smooth < 1:
        out[:] = (1 - smooth) / (len(dict_alphabet) - 1)
    l = len(seq)
    for i, c in enumerate(seq):
        out[dict_alphabet[c], i] = smooth
    out[:, l:] = 0
    return out

def one_hot_decoding(inp, dict_alphabet):
    tmp = np.argmax(inp, axis=0)
    out = ''.join(dict_alphabet[i] for i in tmp)
    return out

class SingleProcess():

    def __init__(self, infile, config_dict, indexBegin = 0):
        self.config_dict = config_dict
        self.h_index = self.config_dict["h_index_nums"]
        self.read_len = self.config_dict["end_tree_len"]
        self.infile = infile
        self.tree = tree.new_tree(self.read_len)
        self.sub_tree_depth = self.read_len
        self.sub_tree = tree.new_tree(self.sub_tree_depth)
        self.indexList = []
        self.clust_num = indexBegin + 0
        self.branch_num = 0
        self.core_set = {}
        self.alphabet = {'A':0, 'T':1, 'G':2, 'C':3}
        self.reversed_alphabet = {0:'A', 1:'T', 2:'G', 3:'C'}
        self.merge_map = {}

    def cluster(self, dna_tag, dna_str):
        dna_str = dna_str[self.h_index:self.h_index+self.read_len]
        if len(dna_str) == self.read_len:
            align_result = tree.quick_search(self.tree, dna_str, self.config_dict['tree_threshold'], 4)
            label = align_result.label
            # If the match is successful, it is recorded.
            if label > 0:
                self.indexList.append((dna_tag, label))
            else:
                # Check umatched sequence in subtree
                sub_dna_str = dna_str[:self.sub_tree_depth]
                sub_align_result = tree.search(self.sub_tree, sub_dna_str, self.config_dict['tree_threshold'])
                if sub_align_result.label > 0:
                    new_label = sub_align_result.label

                else:
                    self.clust_num += 1
                    new_label = self.clust_num
                    tree.insert(self.sub_tree, sub_dna_str, new_label)
                tree.insert(self.tree, dna_str, new_label)
                self.indexList.append((dna_tag, new_label))
        else:
            self.clust_num += 1
            self.indexList.append((dna_tag, self.clust_num))

    def cluster_with_index(self, dna_tag, dna_str):
        if self.h_index == 0:
            dna_str = dna_str[:self.read_len]
        else:
            dna_str = dna_str[self.h_index:self.h_index+self.read_len]

        align_result = tree.search(self.tree, dna_str, self.config_dict['tree_threshold'])

        # If the match is successful, it is recorded.
        if align_result[1] < self.config_dict['tree_threshold']:
            self.indexList.append((dna_tag, align_result[0]))

        else:
            self.clust_num += 1
            self.tree.insert(dna_str[:self.read_len], self.clust_num)
            self.indexList.append((dna_tag, self.clust_num))

    def run(self):
        # Get number of reads
        numReads = mapcount(self.infile)

        print('Number of reads: %d'%numReads)
        with tqdm(total=numReads) as pbar:
            with open(self.infile, 'r') as f:
                m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
                line = m.readline()
                while line != b'':     
                    pbar.update(1)
                    tag, read = line.decode('utf8').strip().replace('N', '').split(' ')
                    self.cluster(int(tag), read)
                    line = m.readline()
        return self.indexList
            
def chunks(arr, m):
    n = int(math.ceil(len(arr) / float(m)))
    return [arr[i:i + n] for i in range(0, len(arr), n)], [i for i in range(0, len(arr), n)]

def clust(data, config_dict):
    p = SingleProcess(data, config_dict)
    return p.run()

def clustMP(infile, config_dict, indexList, indexBegin):
    p = SingleProcess(infile, config_dict, indexBegin)
    p.run()
    for i, index in enumerate(p.indexList):
        indexList[indexBegin + i] = index

if __name__ == '__main__':
    config_dict = lc.out_put_config()
    print(config_dict['tag'])

    # Memory for output
    indexList = []
    print('Start clustering...')
    st = time.time()
    # Start clustering
    indexList = clust(config_dict['input_path'], config_dict)

    print('Time: %f'%(time.time()-st))
    if 'output_file' in config_dict:
        output_file = open(config_dict['output_file'], 'w')
        output_file.truncate(0)
        for pair in indexList:
            output_file.write(str(pair[0]) + ',' + str(pair[1]) + '\n')
        output_file.close()

"""Main program
"""
import time
import math
from tqdm import tqdm
from clust import load_config as lc
from clust.multi_stage import *
from clust import tree
import mmap
def mapcount(filename):
    f = open(filename, "r+")
    buf = mmap.mmap(f.fileno(), 0)
    lines = 0
    readline = buf.readline
    while readline():
        lines += 1
    return lines

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
        self.clust_num = indexBegin
        self.merge_end = 0
        self.read_num = 0
        self.depth_limit = self.config_dict['depth_limit']
        self.main_tree_th = self.config_dict['tree_threshold']
        self.sub_tree_th = self.config_dict['sub_tree_threshold']
        # self.config = tree.Config()
        # self.config.tree_depth = self.read_len
        # self.config.index_begin = indexBegin
        # self.config.main_tree_threshold = self.config_dict['tree_threshold']
        # self.config.sub_tree_threshold = self.config_dict['sub_tree_threshold']
        # self.config.depth_limit = self.config_dict['depth_limit']
        # self.alg = tree.Process(self.config)

    def cluster(self, dna_tag, dna_str):
        # label = self.alg.cluster(dna_str) # C++ interface
        # self.indexList.append((dna_tag, label))
        dna_str = dna_str[self.h_index:self.h_index+self.read_len]
        if len(dna_str) == self.read_len:
            align_result = tree.quick_search(self.tree, dna_str, self.main_tree_th, self.depth_limit)
            label = align_result.label

            if label > 0:
                self.indexList.append((dna_tag, label))
            else:
                # Check umatched sequence in subtree
                if dna_tag < self.merge_end:
                    sub_dna_str = dna_str[:self.sub_tree_depth]
                    sub_align_result = tree.search(self.sub_tree, sub_dna_str, self.sub_tree_th)
                    if sub_align_result.label > 0:
                        new_label = sub_align_result.label

                    else:
                        self.clust_num += 1
                        new_label = self.clust_num
                        tree.insert(self.sub_tree, sub_dna_str, new_label)
                else:
                    self.clust_num += 1
                    new_label = self.clust_num                    
                tree.insert(self.tree, dna_str, new_label)
                self.indexList.append((dna_tag, new_label))
        else:
            self.clust_num += 1
            self.indexList.append((dna_tag, self.clust_num))

    def cluster_with_index(self, dna_tag, dna_str):
        dna_str = dna_str[self.h_index:self.h_index+self.read_len]
        if len(dna_str) == self.read_len:
            align_result = tree.quick_search(self.tree, dna_str, self.main_tree_th, self.depth_limit)
            label = align_result.label

            if label > 0:
                self.indexList.append((dna_tag, label))
            else:
                # Check umatched sequence in subtree
                if dna_tag < self.merge_end:
                    sub_dna_str = dna_str[:self.sub_tree_depth]
                    sub_align_result = tree.search(self.sub_tree, sub_dna_str, self.sub_tree_th)
                    if sub_align_result.label > 0:
                        new_label = sub_align_result.label

                    else:
                        self.clust_num += 1
                        new_label = self.clust_num
                        tree.insert(self.sub_tree, sub_dna_str, new_label)
                else:
                    self.clust_num += 1
                    new_label = self.clust_num                    
                #tree.insert(self.tree, dna_str, new_label)
                self.indexList.append((dna_tag, new_label))
        else:
            self.clust_num += 1
            self.indexList.append((dna_tag, self.clust_num))

    def run(self):
        # Get number of reads
        numReads = mapcount(self.infile)
        self.read_num = numReads
        self.merge_end = int(self.config_dict["frac"] * numReads)
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

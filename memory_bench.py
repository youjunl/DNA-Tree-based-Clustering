from clust import tree
from tqdm import tqdm
from clust.multi_stage import *
from clust import tree
import psutil
import os
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

    def __init__(self, infile, config_dict, main_tree, sub_tree, indexBegin = 0):
        self.config_dict = config_dict
        self.h_index = self.config_dict["h_index_nums"]
        self.read_len = self.config_dict["end_tree_len"]
        self.infile = infile
        self.tree = main_tree
        self.sub_tree_depth = self.read_len
        self.sub_tree = sub_tree
        self.indexList = []
        self.clust_num = indexBegin
        
        self.config = tree.Config()
        self.config.tree_depth = self.read_len
        self.config.index_begin = indexBegin
        self.config.main_tree_threshold = self.config_dict['tree_threshold']
        self.config.sub_tree_threshold = self.config_dict['sub_tree_threshold']
        self.config.depth_limit = 3
        self.alg = tree.Process(self.config)

    def cluster(self, dna_tag, dna_str):
        # label = self.alg.cluster(dna_str)
        # self.indexList.append((dna_tag, label))
        dna_str = dna_str[self.h_index:self.h_index+self.read_len]
        if len(dna_str) == self.read_len:
            align_result = tree.quick_search(self.tree, dna_str, self.config_dict['tree_threshold'], 3)
            label = align_result.label

            if label > 0:
                self.indexList.append((dna_tag, label))
            else:
                # Check umatched sequence in subtree
                if dna_tag > self.config_dict["train_start"] and dna_tag < self.config_dict["train_end"]:
                    sub_dna_str = dna_str[:self.sub_tree_depth]
                    sub_align_result = tree.search(self.sub_tree, sub_dna_str, self.config_dict['sub_tree_threshold'])
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

    def run(self, check_points):
        # Get number of reads
        numReads = mapcount(self.infile)
        usage = [0 for _ in check_points]
        print('Number of reads: %d'%numReads)
        with tqdm(total=numReads) as pbar:
            with open(self.infile, 'r') as f:
                m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
                line = m.readline()
                count = 1
                while line != b'':
                    pbar.update(1)
                    tag, read = line.decode('utf8').strip().replace('N', '').split(' ')
                    self.cluster(int(tag), read)
                    line = m.readline()
                    if count%1000==0:
                        # Output the memory usage
                        ind = count//1000
                        usage[ind] = psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024

                    count += 1

        return usage

if __name__ == '__main__':
    # DLS configurations:
    config_dict={
        "end_tree_len" : 21,
        "tree_threshold" : 4,
        "sub_tree_threshold" : 2,
        "h_index_nums" : 0,
        "train_start": 0,
        "train_end": 0
    }

    # # DLSM configurations:
    # config_dict={
    #     "end_tree_len" : 13,
    #     "tree_threshold" : 3,
    #     "sub_tree_threshold" : 2,
    #     "h_index_nums" : 0,
    #     "train_start": 0,
    #     "train_end": 1000000
    # }

    infile = 'testdata/toClustReads_1000000.txt'
    outfile = 'mem_usage_dls_21.csv'
    check_points = [i for i in range(0, 1000000+1000, 1000)]
    main_tree = tree.new_tree(config_dict['end_tree_len'])
    sub_tree = tree.new_tree(config_dict['end_tree_len'])
    p = SingleProcess(infile, config_dict, main_tree, sub_tree)
    usage = p.run(check_points)
    with open(outfile, 'w') as f:
        f.truncate(0)
        for pts in check_points:
            f.write('%d,'%pts)
        f.write('\n')
        for pts in usage:
            f.write('%.4f,'%pts)
        f.write('\n')        
    


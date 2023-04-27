"""Main program



"""
import time
import math
import random
from multiprocessing import Process, Manager, Queue, Pool
from tqdm import tqdm
from clust import load_config as lc
from clust import pytree as tr
from clust.multi_stage import ssc, ssc_repeat
from clust import tree

class SingleProcess():

    def __init__(self, data, config_dict, indexBegin = 0):
        self.config_dict = config_dict
        self.h_index = self.config_dict["h_index_nums"]
        self.read_len = self.config_dict["end_tree_len"]
        self.data = data
        self.tree = tree.new_tree(self.read_len)
        self.indexList = []
        self.clust_num = indexBegin + 0
        self.branch_num = 0
    def cluster(self, dna_tag, dna_str):
        if self.h_index == 0:
            dna_str = dna_str[:self.read_len]
        else:
            dna_str = dna_str[self.h_index:self.h_index+self.read_len]

        if len(dna_str) == self.read_len:
            align_result = tree.quick_search(self.tree, dna_str, self.config_dict['tree_threshold'], 2)
            label, distance = align_result.label, align_result.distance
            # align_result_slow = tree.search(self.tree, dna_str, self.config_dict['tree_threshold'])
            # print(label, distance, align_result_slow.label, align_result_slow.distance)
            # If the match is successful, it is recorded.
            if label > 0 and distance < self.config_dict['tree_threshold']:
                self.indexList.append((dna_tag, label))

            else:
                self.clust_num += 1
                self.branch_num += 1
                tree.insert(self.tree, dna_str, self.clust_num)
                self.indexList.append((dna_tag, self.clust_num))
        else:
            return

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
            self.branch_num += 1
            self.tree.insert(dna_str[:self.read_len], self.clust_num)
            self.indexList.append((dna_tag, self.clust_num))

    def train(self, train_num):
        print('Training...')
        samples = random.sample(self.data, min(train_num, len(self.data)))

        for _, seq in tqdm(samples):
            result = tree.search(self.tree, seq[:self.read_len], self.config_dict['tree_threshold'])
            if result.label == -1:
                self.clust_num += 1
                tree.insert(self.tree, seq[:self.read_len], self.clust_num)
                self.branch_num += 1
            elif result.label > 0 and result.distance <= self.config_dict['tree_threshold']:
                tree.insert(self.tree, seq[:self.read_len], result.label) # Merging

        print('Test %d reads, %d clusters have been created'%(train_num, self.branch_num))

    def run(self):
        if self.config_dict['use_index']:
            with open(self.config_dict['index_file'], 'r') as f:
                for line in f.readlines():
                    seq = line.strip()
                    self.clust_num += 1
                    self.tree.insert(seq, self.clust_num)
                self.read_len = len(seq)

            for tag, read in tqdm(self.data):
                self.cluster_with_index(tag, read)

        else:
            # Train clustering model
            self.train(self.config_dict['train_num'])
            # Clustering
            for tag, read in tqdm(self.data):
                self.cluster(tag, read)
            
def chunks(arr, m):
    n = int(math.ceil(len(arr) / float(m)))
    return [arr[i:i + n] for i in range(0, len(arr), n)], [i for i in range(0, len(arr), n)]

def clust(data, config_dict):
    p = SingleProcess(data, config_dict)
    p.run()
    return p.indexList

def clustMP(data, config_dict, indexList, indexBegin):
    p = SingleProcess(data, config_dict, indexBegin)
    p.run()
    for i, index in enumerate(p.indexList):
        indexList[indexBegin + i] = index

if __name__ == '__main__':
    config_dict = lc.out_put_config()
    print(config_dict['tag'])

    # Open file
    f = open(config_dict['input_path'], "r")

    # Read file into the memory
    lines = f.readlines()
    data = []
    for line in lines:
        tag, read = line.strip().replace('N', '').split(' ')
        data.append((int(tag), read))
    # Memory for output
    indexList = []
    st = time.time()
    # Start clustering
    if not config_dict["multi_stage"]:
        indexList = clust(data, config_dict)

    else: # Multi processor mode is only valid for multi stage clustering
        nProcess = config_dict['processes_nums']
        if nProcess < 2:
            indexList = clust(data, config_dict)          
        else:
            splitData, indexBegin = chunks(data, nProcess)
            Processes = []
            manager = Manager()
            indexListMP = manager.list([[] for _ in range(len(data))])
            for i in range(nProcess):
                p = Process(target=clustMP, args=(splitData[i], config_dict, indexListMP, indexBegin[i]))
                p.start()
                Processes.append(p)

            for p in Processes:
                p.join()
                indexList = indexListMP
        # Multi stage clustering
        start_tree_depth = config_dict["end_tree_len"]
        tree_depth = start_tree_depth
        stageNum = config_dict['extra_stage_num']
        tree_threshold = config_dict["tree_threshold"]
        spliterFlag = config_dict["spliter"]
        filterFlag = config_dict["filter"]
        for i in range(stageNum):
            indexList = ssc_repeat(indexList, data, tree_depth, tree_threshold, filter=filterFlag, spliter=spliterFlag)
            tree_threshold += 2
            print('st %d'%(i+2))

    print('Time: %f'%(time.time()-st))
    if 'output_file' in config_dict:
        output_file = open(config_dict['output_file'], 'w')
        output_file.truncate(0)
        for pair in indexList:
            output_file.write(str(pair[0]) + ',' + str(pair[1]) + '\n')
        output_file.close()

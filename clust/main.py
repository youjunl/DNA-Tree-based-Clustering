"""Main program



"""
import time
import math
from multiprocessing import Process, Manager, Queue, Pool
from tqdm import tqdm
from clust import load_config as lc
from clust import tree as tr
from clust.multi_stage import ssc, ssc_repeat


class SingleProcess():

    def __init__(self, data, config_dict, indexBegin = 0):
        self.config_dict = config_dict
        self.h_index = self.config_dict["h_index_nums"]
        self.read_len = self.config_dict["end_tree_len"]
        self.data = data
        self.tree = tr.Trie()
        self.indexList = []
        self.test_num = 0
        self.clust_num = indexBegin + 0
        self.insertion = 0
        self.deletion = 0
        self.substitution = 0
        self.chCnt = 0

    def cluster(self, dna_tag, dna_str):
        self.chCnt += len(dna_str)
        if self.h_index == 0:
            dna_str = dna_str[:self.read_len]
        else:
            dna_str = dna_str[self.h_index:self.h_index+self.read_len]

        align_result = self.tree.fuzz_fin(
            dna_str, self.config_dict['tree_threshold'])

        # If the match is successful, it is recorded.
        if sum(align_result[1]) < self.config_dict['tree_threshold']:
            self.indexList.append((dna_tag, align_result[0]))
            self.insertion += align_result[1][0]
            self.deletion += align_result[1][1]
            self.substitution += align_result[1][2]

        else:
            self.clust_num += 1
            self.tree.insert(dna_str, self.clust_num)
            self.indexList.append((dna_tag, self.clust_num))

    # Process flow
    def run(self):
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
            tree_depth -= 2
            print('st %d'%(i+2))
    print('Time: %f'%(time.time()-st))
    if 'output_file' in config_dict:
        output_file = open(config_dict['output_file'], 'w')
        output_file.truncate(0)
        for result in indexList:
            output_file.write(str(result[0]) + ',' + str(result[1]) + '\n')
        output_file.close()

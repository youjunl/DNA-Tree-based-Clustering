"""Main program



"""
import os
import time
from collections import  Counter
from multiprocessing import Process, Queue
from tqdm import tqdm
from clust import load_config as lc
from clust import tree as tr

class MyProcess():

    def __init__(self, data):
        self.config_dict=lc.out_put_config()
        self.data = data
        self.tree = tr.Trie()
        self.Cluster_size_threshold = self.config_dict['Cluster_size_threshold']
        self.ref_list={}    
        self.ref_dict={}     
        self.ref_error_dict={}
        self.num_dict={}
        self.tag_dict={}
        self.index_list=[]   
        self.now_clust_threshold = self.config_dict['now_clust_threshold']
        self.read_len= self.config_dict['read_len']  
        self.dntree_nums = self.config_dict['end_tree_len']
        self.fuzz_list = [self.config_dict['thd_tree_loc'],self.config_dict['four_tree_loc'],self.config_dict['other_tree_len']] 
        self.loc_nums = self.config_dict['Vertical_drift'] 
        self.tag_nums = self.config_dict['tag_nums']     
        self.align_swicth= self.config_dict['align_fuc'] 
        self.fuzz_tree_nums= self.config_dict['Horizontal_drift']   
        self.h_index=self.config_dict['h_index_nums']
        self.e_index=self.config_dict['e_index_nums']
        self.test_num = 0
        self.file_format = "txt"
        if 'input_path' in self.config_dict:
            if self.config_dict['input_path'][-1] == "a" :
                self.file_format = "fasta"
            elif self.config_dict['input_path'][-1] == "q" :
                self.file_format = "fastq"
        self.clust_num = 0


    def cluster(self, dna_tag, dna_str):
        self.test_num += 1
        if self.h_index == 0 :
            dna_str = dna_str[:self.dntree_nums]
        else:
            dna_str = dna_str[self.h_index:self.h_index+self.dntree_nums]

        if "N" in dna_str:
            pass 
        else:
            align_result = self.tree.match(dna_str,self.config_dict['tree_threshold']) 

            if align_result > 0:  #If the match is successful, it is recorded.
                self.index_list.append((dna_tag, align_result))
                
            else:
                self.clust_num += 1
                self.tree.insert(dna_str,self.clust_num)
                self.index_list.append((dna_tag, self.clust_num))

    #Process flow
    def run(self):
        for line in tqdm(self.data):
            if line == [] or line == "" :
                break
            tag, read = line.strip().split(' ')
            self.cluster(tag, read)

def main():
    pass

if __name__ == '__main__':
    config_dict = lc.out_put_config()
    f=open(config_dict['input_path'],"r")
    lines = f.readlines()
    st = time.time()
    p = MyProcess(lines)
    p.run()
    print("Time:",time.time()-st)


    if 'output_file' in config_dict:
        output_file = open(config_dict['output_file'], 'w')
        output_file.truncate(0)
        for result in p.index_list:
            output_file.write(str(result[0]) + ',' + str(result[1]) + '\n')
        output_file.close()

    main()
"""Main program

This module is the main code of Clover, which contains the algorithm flow,
multi-process implementation and data statistics of Clover.
"""
import os
import time
from collections import  Counter
from multiprocessing import Process, Queue
from tqdm import tqdm
from clust import load_config as lc
from clust import tree as tr

class MyProcess(Process):

    def __init__(self, name, data, q_output):
        Process.__init__(self)
        self.config_dict=lc.out_put_config()
        self.name = name
        self.data = data
        self.q_output = q_output
        self.q_output_temp = []
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



    def cluster(self,read):
        self.test_num=self.test_num+1
        line_=read.split()

        dna_num=self.test_num
        if self.config_dict['Virtual_mode'] == False or self.config_dict['mmr_mode'] == True:
            dna_index=line_[0]
            dna_str = line_[1]
            dna_tag=line_[0]
        else:
            dna_str = line_[1]
            dna_tag = line_[0]
        if self.h_index == 0 :
            dna_str = dna_str[:self.dntree_nums]
        else:
            dna_str = dna_str[self.h_index:self.h_index+self.dntree_nums]

        if "N" in dna_str:#dna_str_num < self.config_dict['read_len_min'] or "N" in dna_str:
            pass 
        else:
            align_result = self.tree.fuzz_fin(dna_str,self.config_dict['tree_threshold']) 

            if align_result[1]<self.fuzz_tree_nums:  #If the match is successful, it is recorded.
                self.index_list.append((dna_index,self.ref_dict[align_result[0]][0]))
                
            else:
                if self.config_dict['align_fuc'] == True:
                    self.ref_list[dna_num]=dna_str
                self.ref_dict[dna_num]=[dna_tag]
                self.tree.insert(dna_str,dna_num)
                self.index_list.append((dna_index,self.ref_dict[dna_num][0]))
    #Process flow
    def run(self):

        self.num_dict[self.name+"sum_read_num"]=0
        self.num_dict[self.name+"error_num"]=0
        self.num_dict[self.name+"sum_cluster_num"]=0
        self.num_dict[self.name+"sum_tag"]=[]
        len_=len(self.name) 
    
        if self.config_dict['fast_mode'] == True :
            
            for line in tqdm(self.data):
                self.cluster(line)
            if 'output_file' in self.config_dict :

                self.num_dict[self.name+"index_list"]=self.index_list
            if self.config_dict['Virtual_mode'] == True :
                tag_sum=0
                tag_error=0
                nums_=0
                nums_sum=0
                if self.config_dict['Virtual_mode'] == True:
                    for key in self.ref_dict:
                        tag_len=len(self.ref_dict[key])
                        
                        if tag_len > self.Cluster_size_threshold :
                            nums_sum=nums_sum+1 
                            tag_c=Counter(self.ref_dict[key]).most_common(1)[0][1]
                            tag_error=tag_error+tag_len-tag_c
                            tag_sum=tag_sum+tag_len   
                            tag=self.ref_dict[key][0]
                            if tag in self.tag_dict:
                                pass
                            else:
                                nums_=nums_+1
                                self.tag_dict[tag]=1
                        self.ref_dict[key]=[]
                    
                    self.num_dict[self.name+"sum_read_num"]+=tag_sum
                    self.num_dict[self.name+"error_num"]+=tag_error
                    self.num_dict[self.name+"sum_cluster_num"]=nums_sum
                    for key in self.tag_dict:
                        self.num_dict[self.name+"sum_tag"].append(key)
        elif self.config_dict['mmr_mode'] == True:

            output_nums = 0
            data_file=open(self.config_dict['input_path'],"r")
            nums_line = 0
            if 'output_file' in self.config_dict :
                output_file = open(self.name+"_"+self.config_dict['output_file'],'a')
                output_file.write("[")
            while True :
                output_nums += 1
                line = data_file.readline()
                if line == [] or line == "" :
                    break
                if self.file_format ==  "txt" :
                    dna_reads = line
                elif self.file_format ==  "fasta" :
                    if nums_line%2 == 0 :
                        read_name = line[1:]
                        nums_line += 1
                        continue
                    if nums_line%2 == 1:
                        dna_reads = read_name + " " + line
                        nums_line += 1
                elif self.file_format ==  "fastq" :
                    if nums_line%4 == 0 :
                        read_name = line[1:]
                        nums_line += 1
                        continue
                    elif nums_line%4 == 1 :
                        dna_reads = read_name + " " + line
                        nums_line += 1
                    else:
                        continue
                if int(self.config_dict['processes_nums']) == 0:
                    self.cluster(dna_reads)
                elif line.split()[1][:len_] == self.name :
                    self.cluster(dna_reads)
                if 'output_file' in self.config_dict and output_nums%100000 == 0:
                    for i in self.index_list:
                        output_file.write(str(i))
                        output_file.write(',')
                        self.index_list=[]
            if 'output_file' in self.config_dict :
                for i in self.index_list:
                    output_file.write(str(i))
                    output_file.write(',')
                output_file.close()
                with open(self.name+"_"+self.config_dict['output_file'], 'rb+') as output_file:
                    output_file.seek(-1, os.SEEK_END)
                    output_file.truncate()
                output_file = open(self.name+"_"+self.config_dict['output_file'],'a')
                output_file.write("]")    
                output_file.close()
        else :
            if self.config_dict['Virtual_mode'] == True :
                f=open(self.config_dict['input_path'],"r")
                while True :
                    line = f.readline()
                    if line == [] or line == "" :
                        break
                    if int(self.config_dict['processes_nums']) == 0:
                        self.cluster(line)
                    elif line.split()[1][:len_] == self.name :
                        self.cluster(line)

                if 'output_file' in self.config_dict :
                    self.num_dict[self.name+"index_list"]=self.index_list
                tag_sum=0
                tag_error=0
                nums_=0
                nums_sum=0
                for key in self.ref_dict:
                    tag_len=len(self.ref_dict[key])
                    
                    if tag_len > self.Cluster_size_threshold :
                        nums_sum=nums_sum+1 
                        tag_c=Counter(self.ref_dict[key]).most_common(1)[0][1]
                        tag_error=tag_error+tag_len-tag_c
                        tag_sum=tag_sum+tag_len   
                        tag=self.ref_dict[key][0]
                        if tag in self.tag_dict:
                            pass
                        else:
                            nums_=nums_+1
                            self.tag_dict[tag]=1
                    self.ref_dict[key]=[]
                
                self.num_dict[self.name+"sum_read_num"]+=tag_sum
                self.num_dict[self.name+"error_num"]+=tag_error
                self.num_dict[self.name+"sum_cluster_num"]=nums_sum
                for key in self.tag_dict:
                    self.num_dict[self.name+"sum_tag"].append(key)
            else:
                f=open(self.config_dict['input_path'],"r")
                nums_line = 0
                while True :
                    line = f.readline()

                    if line == [] or line == "" :
                        break
                    if self.file_format ==  "txt" :
                        dna_reads = line
                    elif self.file_format ==  "fasta" :
                        if nums_line%2 == 0 :
                            read_name = line[1:]
                            nums_line += 1
                            continue
                        if nums_line%2 == 1:
                            dna_reads = read_name + " " + line
                            nums_line += 1
                    elif self.file_format ==  "fastq" :
                        if nums_line%4 == 0 :
                            read_name = line[1:]
                            nums_line += 1
                            continue
                        elif nums_line%4 == 1 :
                            dna_reads = read_name + " " + line
                            nums_line += 1
                        else:
                            continue
                    if int(self.config_dict['processes_nums']) == 0:
                        self.cluster(dna_reads)
                    elif line.split()[1][:len_] == self.name :
                        self.cluster(dna_reads)
                    if 'output_file' in self.config_dict :
                        self.num_dict[self.name+"index_list"]=self.index_list
        self.q_output.put(self.num_dict)


def all_permutations(items, length):
    res = []
    def track_back(tmp_permutation):
        if len(tmp_permutation) == length:
            res.append(tmp_permutation[:])
        else:
            for i in items:
                tmp_permutation.append(i)
                track_back(tmp_permutation)
                tmp_permutation.pop()
    track_back([])
    return ["".join(i) for i in res]

def main():
    pass

if __name__ == '__main__':

    #******************************************************************************************
    config_dict=lc.out_put_config()

    tag_nums=config_dict['tag_nums']
    PROCESS_INDEX= int(config_dict['processes_nums'])
    
    #******************************************************************************************
    
    print(config_dict['tag'])

    if PROCESS_INDEX == 0 and config_dict['fast_mode'] == True :
        N_PROCESS=1
        process_names=['all']
        print("Generating data")
        data_dict = {}
        data_dict['all']=[]
        f= open(config_dict['input_path'],'r').readlines()
        for i in tqdm(f):
            if "N" in i :
                pass
            else:
                data_dict['all'].append(i)
    elif PROCESS_INDEX == 0 and config_dict['fast_mode'] == False :
        N_PROCESS=1
        process_names=['all']
        data_dict = {}
        data_dict['all']=[]
    elif PROCESS_INDEX > 0 :
        N_PROCESS = 4**PROCESS_INDEX
        process_names = all_permutations(["A","T","G","C"],PROCESS_INDEX)


        
        st_pre = time.time()
        data_dict = {}.fromkeys(process_names)
        for i in data_dict:
            data_dict[i] = []


        if config_dict['fast_mode'] == True:
            print("Pre-process the data")
            f= open(config_dict['input_path'],'r').readlines()
            for i in tqdm(f):
                if "N" in i or "*" in i :
                    pass
                else:
                    data_dict[i.split()[1][:PROCESS_INDEX]].append(i)
        else:
            data_dict[i] = []

    q_output = Queue(N_PROCESS*2)

    st = time.time()
    process_dict = {}.fromkeys(process_names)

    for i in process_dict:
        process_dict[i] = MyProcess(i, data_dict[i], q_output)

    for i in process_dict:
        process_dict[i].start()

    st = time.time()
    count_dict={}
    i=0
    while True :
        if i == N_PROCESS :
            break
        count_dict.update(q_output.get())
        i=i+1
    print("Time:",time.time()-st)

    for i in process_dict:
        process_dict[i].join()    
    
    new_count_dict={}
    new_count_dict["sum_read_num"]=0
    new_count_dict["error_num"]=0
    new_count_dict["sum_cluster_num"]=0
    new_count_dict["sum_tag"]={}
    new_count_dict["index_list"]=[]
    for key in process_dict :
        new_count_dict["sum_read_num"]+=count_dict[key+"sum_read_num"]
        new_count_dict["error_num"]+=count_dict[key+"error_num"]
        new_count_dict["sum_cluster_num"]+=count_dict[key+"sum_cluster_num"]
        tag_list=count_dict[key+"sum_tag"]
        
        for i in tag_list:
            if i in new_count_dict["sum_tag"]:
                pass
            else:
                new_count_dict["sum_tag"][i]=1
    if 'output_file' in config_dict and config_dict["mmr_mode"] is not True:
        for key in process_dict :
            new_count_dict["index_list"]+=count_dict[key+"index_list"]
        output_file = open(config_dict['output_file'],'w')
        for result in new_count_dict["index_list"]:
            output_file.write(result[0] + ',' + result[1] + '\n')
        output_file.close()

    if config_dict['Virtual_mode'] == True :
        print("Number of reads processed:",new_count_dict["sum_read_num"])
        print("Accuracy：",(new_count_dict["sum_read_num"]-new_count_dict["error_num"])/new_count_dict["sum_read_num"])
        if config_dict['tag_mode'] == True :
            print("Coverage：",len(new_count_dict["sum_tag"])/int(tag_nums))
            print("Redundancy Rate：",(new_count_dict["sum_cluster_num"]-int(tag_nums))/int(tag_nums))
    elif config_dict['Statistical_model'] == True :
        pass

    main()
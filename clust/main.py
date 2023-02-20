"""Main program



"""
import time
from tqdm import tqdm
from clust import load_config as lc
from clust import tree as tr

class MyProcess():

    def __init__(self, data):
        self.config_dict=lc.out_put_config()
        self.h_index = self.config_dict["h_index_nums"]
        self.dntree_nums = self.config_dict["end_tree_len"]
        self.data = data
        self.tree = tr.Trie()
        self.index_list=[]
        self.test_num = 0
        self.clust_num = 0
        self.insertion = 0
        self.deletion = 0
        self.substitution = 0
        self.chCnt = 0

    def cluster(self, dna_tag, dna_str):
        self.test_num += 1
        dna_str = dna_str.strip('N')
        self.chCnt += len(dna_str)
        if self.h_index == 0 :
            dna_str = dna_str[:self.dntree_nums]
        else:
            dna_str = dna_str[self.h_index:self.h_index+self.dntree_nums]

        align_result = self.tree.fuzz_fin(dna_str,self.config_dict['tree_threshold']) 
        
        if sum(align_result[1]) < self.config_dict['tree_threshold']:  #If the match is successful, it is recorded.
            self.index_list.append((dna_tag, align_result[0]))
            self.insertion += align_result[1][0]
            self.deletion += align_result[1][1]
            self.substitution += align_result[1][2]

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
    print(config_dict['tag'])
    
    # Open file
    f=open(config_dict['input_path'],"r")

    # Read file into the memory
    lines = f.readlines()

    # Timer
    st = time.time()
    
    # Start clustering
    p = MyProcess(lines)
    p.run()
    print("Time:",time.time()-st)
    print("Estimated Channel:")
    print("I: %d, %.4f%%"%(p.insertion, p.insertion/p.test_num/p.dntree_nums))
    print("D: %d, %.4f%%"%(p.deletion, p.deletion/p.test_num/p.dntree_nums))
    print("S: %d, %.4f%%"%(p.substitution, p.substitution/p.test_num/p.dntree_nums))
    if 'output_file' in config_dict:
        output_file = open(config_dict['output_file'], 'w')
        output_file.truncate(0)
        for result in p.index_list:
            output_file.write(str(result[0]) + ',' + str(result[1]) + '\n')
        output_file.close()

    main()
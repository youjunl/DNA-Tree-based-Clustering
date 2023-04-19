"""Module for input parameters.

This module is responsible for importing the input parameters to Clover.

"""

import json
import getopt
import sys

config_dict={
    "end_tree_len" : 21,
    "tree_threshold" : 5,
    "h_index_nums" : 0,
    "processes_nums" : 1,
    "multi_stage": False,
    "reduce_size": 2,
    "spliter": False,
    "filter": False,
    "use_index": False,
    "index_file": ''
}

opt,args = getopt.getopt(sys.argv[1:],'-I:-S:-H:-D:-P:-O:-R:-F:-h',['help','multi-stage','enable-spliter'])

#Read input info
def load_json(path):
    lines = []
    with open(path) as f :
        for row in f.readlines():
            if row.strip().startswith('//'):
                continue
            lines.append(row)
    return json.loads("\n".join(lines))

#Generate a sequence of vertical drifts
def generate_vertical_drifts_list(x):
    list=[]
    i=0
    k=-x-1
    while True:
        k=k+1
        list.append(k)
        i+=1
        if k == x :
            break
    return list

#Write the input to config.json
helpInfo = '''
DNA Tree Based Clustering
#################################
Command: python -m clust.main -I [input file] -O [output_file] ... [options]
-I file to be clustered
-O output result

Options:
-D tree depth
-H tree threshold
-P number of processors
-S number of extra stages of clustering
-R reduce size
-F index file path
--multi-stage enable multistage clustering 
'''
def out_put_config():
    opt,args = getopt.getopt(sys.argv[1:],'-I:-S:-H:-D:-P:-O:-R:-F:-h',['help','multi-stage','enable-spliter'])

    for opt_name,opt_value in opt :
        if '-h' in opt_name or '--help' in opt_name:
            print(helpInfo)
        if '-I' in opt_name :
            config_dict['input_path'] = opt_value
        if '-H' in opt_name :
            config_dict['tree_threshold'] = int(opt_value)
        if '-P' in opt_name :
            config_dict['processes_nums'] = int(opt_value)
        if '-O' in opt_name :
            config_dict['output_file'] = opt_value
        if '-D' in opt_name :
            config_dict['end_tree_len'] = int(opt_value)
        if '-S' in opt_name:
            config_dict['multi_stage'] = True
            config_dict['extra_stage_num'] = int(opt_value)
        if '-R' in opt_name:
            config_dict['reduce_size'] = int(opt_value)
        if '-F' in opt_name:
            config_dict['index_file'] = opt_value
            config_dict['use_index'] = True
        if 'multi-stage' in opt_name:
            config_dict['multi_stage'] = True
            if 'extra_stage_num' not in config_dict.keys():
                config_dict['extra_stage_num'] = 1
        if 'enable-spliter' in opt_name:
            config_dict['spliter'] = True
        if 'enable-filter' in opt_name:
            config_dict['filter'] = True
    config_dict['tag']="""
    
 _____  _      _   _  _____  _____  _____ ______  _____  _   _  _____      ______  _____ ___  ___ _____ 
/  __ \| |    | | | |/  ___||_   _||  ___|| ___ \|_   _|| \ | ||  __ \     |  _  \|  ___||  \/  ||  _  |
| /  \/| |    | | | |\ `--.   | |  | |__  | |_/ /  | |  |  \| || |  \/     | | | || |__  | .  . || | | |
| |    | |    | | | | `--. \  | |  |  __| |    /   | |  | . ` || | __      | | | ||  __| | |\/| || | | |
| \__/\| |____| |_| |/\__/ /  | |  | |___ | |\ \  _| |_ | |\  || |_\ \     | |/ / | |___ | |  | |\ \_/ /
 \____/\_____/ \___/ \____/   \_/  \____/ \_| \_| \___/ \_| \_/ \____/     |___/  \____/ \_|  |_/ \___/ 
                                                                                                    
                                                                                                    

    """
    return config_dict


"""Module for input parameters.

This module is responsible for importing the input parameters to Clover.

"""

import json
import getopt
import sys

config_dict={
    "end_tree_len" : 20,
    "tree_threshold" : 4,
    "sub_tree_threshold" : 1,
    "depth_limit" : 4,
    "h_index_nums" : 0,
    "use_index": False,
    "index_file": '',
    "frac": 0,
}

#Read input info
def load_json(path):
    lines = []
    with open(path) as f :
        for row in f.readlines():
            if row.strip().startswith('//'):
                continue
            lines.append(row)
    return json.loads("\n".join(lines))

#Write the input to config.json
helpInfo = '''
DNA Tree Based Clustering
Latest Update: 2023-5-26
#################################
Command: python -m clust.main -I [input file] -O [output_file] ... [options]
-I file to be clustered
-O output result

Options:
-D tree depth
-H tree threshold
-S secondary tree threshold
-F index file path
-T fraction of reads that enables merging
-A starting position
'''
def out_put_config():
    opt,args = getopt.getopt(sys.argv[1:],'-I:-S:-H:-D:-P:-O:-A:-R:-F:-T:-h',['help'])

    for opt_name,opt_value in opt :
        if '-h' in opt_name or '--help' in opt_name:
            print(helpInfo)
        if '-I' in opt_name :
            config_dict['input_path'] = opt_value
        if '-H' in opt_name :
            config_dict['tree_threshold'] = int(opt_value)
        if '-O' in opt_name :
            config_dict['output_file'] = opt_value
        if '-D' in opt_name :
            config_dict['end_tree_len'] = int(opt_value)
        if '-S' in opt_name:
            config_dict['sub_tree_threshold'] = int(opt_value)
        if '-A' in opt_name:
            config_dict['h_index_nums'] = int(opt_value)
        if '-F' in opt_name:
            config_dict['index_file'] = opt_value
            config_dict['use_index'] = True
        if '-T' in opt_name:
            config_dict['frac'] = float(opt_value)
        if 'multi-stage' in opt_name:
            config_dict['multi_stage'] = True
            if 'extra_stage_num' not in config_dict.keys():
                config_dict['extra_stage_num'] = 1
    config_dict['tag']="""
    
    Depth-Limited Searching for DNA Clustering                                                             

    """
    return config_dict


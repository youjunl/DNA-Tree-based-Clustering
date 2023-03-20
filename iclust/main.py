from iclust import core
from tqdm import tqdm
import getopt, sys
import os
import time
#Write the input to config.json
helpInfo = '''
#############################################################################################
#Implementation of iterative clustering in paper:                                           #
#Rashtchian, Cyrus, et al. "Clustering billions of reads for DNA data storage." NIPS 2017.  #
#############################################################################################
Command: python -m iclust.main -I [input file] -O [output_file] ... [options]
-I file to be clustered
-O output result

Options:
-R Threshold for edit distance
-Q k-mers
-W Anchor's length
-L Hash length
-N Real segmentation length
-TL Higher threshold for hamming distance
-TH Lower threshold for hamming distance
-SL Number of local steps
-SC Number of communication steps
-C Number of cores
'''

params = {
    'r': 10, # Threshold for edit distance
    'q': 3, # Substring's length
    'w': 4, # Anchor's length
    'l': 6, # Hash's length
    'n': 150, # Real segmentation length
    'theta_low': 10,
    'theta_high': 20,
    'local_step': 15,
    'comm_step': 24,
    'core_num' : 1,
    'block_len' : 25,
    'thread': 1,
}

if __name__ == '__main__':
    opts,args = getopt.getopt(sys.argv[1:],'-I:-O:-R:-Q:-W:-L:-N:-TL:-TH:-SL:-SC:-C:P:-h',['help'])

    for opt_name,opt_value in opts:
        if '-h' in opt_name or '--help' in opt_name:
            print(helpInfo)
        if '-I' in opt_name :
            infile = opt_value
        if '-O' in opt_name :
            outfile = opt_value
        if '-R' in opt_name :
            params['r'] = int(opt_value)
        if '-Q' in opt_name :
            params['q'] = int(opt_value)
        if '-W' in opt_name :
            params['w'] = int(opt_value)
        if '-L' in opt_name :
            params['l'] = int(opt_value)
        if '-N' in opt_name :
            params['n'] = int(opt_value)
        if '-TL' in opt_name :
            params['theta_low'] = int(opt_value)
        if '-TH' in opt_name :
            params['theta_high'] = int(opt_value)
        if '-SL' in opt_name :
            params['local_step'] = int(opt_value)
        if '-SC' in opt_name :
            params['comm_step'] = int(opt_value)
        if '-P' in opt_name :
            params['thread'] = int(opt_value)            
    try:
        print('input: %s'%infile)
        print('output: %s'%outfile)
        print('-------params--------')
        for opt_name, opt_value in opts:
            print(opt_name, opt_value)
    except:
        print('Wrong command')
        print(helpInfo)
        sys.exit(2)

    print('Reading file...')

    inData = []
    cnt = 0    
    with open(infile, 'r') as f:
        for text in tqdm(f.readlines()):
            tag, seq = text.strip().replace('N', '').split(' ')
            data = core.Sequence(seq[:params['n']], int(tag), cnt)
            inData.append(data)
            cnt += 1

    # Run clustering
    begin = time.time()
    process = core.Process(params)
    algout = process.compute_comm(inData)

    # Sort ans
    algout = sorted(algout, key=lambda k: (k[0].tag, k[0].cluster))
    print('Output %d clusters.'%len(algout))
    with open(outfile, 'w') as f:
        f.truncate(0)
        for clust in algout:
            for seq in clust:
                f.write('%d,%d\n'%(seq.tag, seq.cluster))
    
    print('Finished. Time: %d seconds'%(time.time()-begin))
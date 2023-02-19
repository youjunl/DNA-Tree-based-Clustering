import getopt, sys
import numpy as np
text = '''
#############################################################################################
#Shuffle Data                                                                               #
#############################################################################################
'''

if __name__ == '__main__':
    print(text)
    # Inputs Management
    helpInfo = '[Usage]\DataGen.py <labeled data> <output reads> <tags>'
    try:
        opts, args = getopt.getopt(sys.argv[1:],"h",[])
        for opt, arg in opts:
            print(opt,arg)
            if opt in ['-h']:
                print(helpInfo)
                sys.exit()
    
    except getopt.GetoptError:
        print('Error!')
        print(helpInfo)
        sys.exit(2)

    if len(args)<3:
        print('No enough inputs, labeled data, output filenames are expected')
        print(helpInfo)
        sys.exit(2)

    labeled = args[0]
    readfile = args[1]
    tagfile = args[2]

    f = open(labeled, 'r')
    reads = []
    tags = []
    print('Reading files...')
    for line in f.readlines():
        tag, read = line.strip().split(' ')
        reads.append(read)
        tags.append(tag)
    f.close()
    print('Shuffleing results')
    nRead = len(reads)
    inds = np.arange(nRead)
    np.random.shuffle(inds)
    reads = np.array(reads)[inds]
    tags = np.array(tags)[inds]
    print('Write data')
    with open(readfile, 'w') as f:
        f.truncate(0)
        for i, read in enumerate(reads):
            f.write(str(i+1) + ' ' + read + '\n')
    with open(tagfile, 'w') as f:
        f.truncate(0)
        for i, tag in enumerate(tags):
            f.write(str(i+1) + ',' + tag + '\n')
    print('Finished!')
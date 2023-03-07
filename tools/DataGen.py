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
    print('Shuffling results')
    nRead = len(reads)
    inds = np.arange(nRead)
    np.random.shuffle(inds)
    print('Write data')

    # Open target file
    fread = open(readfile, 'w')
    ftag = open(tagfile, 'w')

    # clear
    fread.truncate(0)
    ftag.truncate(0)

    for ind in inds:
        fread.write(str(ind+1) + ' ' + str(reads[ind]) + '\n')
        ftag.write(str(ind+1) + ',' + str(tags[ind]) + '\n')
    fread.close()
    ftag.close()
    print('Finished!')
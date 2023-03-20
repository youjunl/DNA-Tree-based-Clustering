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
    keep_order = False
    try:
        opts, _ = getopt.getopt(sys.argv[1:],"-h:-I:-R:-T:",['keep-order', 'help'])
        for opt, arg in opts:
            print(opt,arg)
            if opt in ['-h']:
                print(helpInfo)
                sys.exit()
            if '-I' in opt:
                labeled = arg
            if '-R' in opt:
                readfile = arg
            if '-T' in opt:
                tagfile = arg
            if 'keep-order' in opt:
                keep_order = True

    except getopt.GetoptError:
        print('Error!')
        print(helpInfo)
        sys.exit(2)

    if len(opts)<3:
        print('No enough inputs, labeled data, output filenames are expected')
        print(helpInfo)
        sys.exit(2)

    f = open(labeled, 'r')
    reads = []
    tags = []
    print('Reading files...')
    for line in f.readlines():
        tag, read = line.strip().split(' ')
        reads.append(read)
        tags.append(tag)
    f.close()
    nRead = len(reads)
    inds = np.arange(nRead)
    if not keep_order:
        print('Shuffling results')
        np.random.shuffle(inds)

    print('Writing data...')

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
    print('Finished!\nReads has been saved to %s\nLabels has been saved to %s'%(readfile, tagfile))
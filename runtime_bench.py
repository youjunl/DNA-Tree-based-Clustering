import subprocess
import time
infiles = ['']

outfile_prefix = 'testdata/toClustReads_'

sample_num = [10000, 100000, 1000000, 10000000]
outfile_names = []
for num in sample_num:
    outfile_name = outfile_prefix+str(num)+'.txt'
    outfile_names.append(outfile_name)

    command = "python -m clust.main -I {}".format(outfile_name)
    begin = time.time()
    subprocess.call (command,shell=True)
    end = time.time()
    print(end-begin)
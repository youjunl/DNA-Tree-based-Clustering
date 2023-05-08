import subprocess
import time
infiles = ['']

outfile_prefix = 'testdata/toClustReads_'

sample_num = [100000, 1000000, 10000000]
outfile_names = []
depth = [16, 18, 20]
runtime_dls = [[0 for _ in sample_num] for _ in depth]
for i, d in enumerate(depth):
    for j, num in enumerate(sample_num):
        outfile_name = outfile_prefix+str(num)+'.txt'
        outfile_names.append(outfile_name)

        command = "python -m clust.main -I {} -D {}".format(outfile_name, d)
        begin = time.time()
        subprocess.call (command,shell=True)
        end = time.time()
        runtime_dls[i][j] = end - begin

runtime_dlsm = [[0 for _ in sample_num] for _ in depth]
for i, d in enumerate(depth):
    for j, num in enumerate(sample_num):
        outfile_name = outfile_prefix+str(num)+'.txt'
        outfile_names.append(outfile_name)

        command = "python -m clust.main -I {} -D {} -T 1".format(outfile_name, d)
        begin = time.time()
        subprocess.call (command,shell=True)
        end = time.time()
        runtime_dlsm[i][j] = end - begin

print(runtime_dls)
print(runtime_dlsm)


    
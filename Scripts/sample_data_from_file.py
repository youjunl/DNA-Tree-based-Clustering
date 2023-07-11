infile = 'testdata/toClustReads.txt'
outfile_prefix = 'testdata/toClustReads_'

sample_num = [10000, 100000, 1000000, 10000000]
outfiles = []
for num in sample_num:
    outfile_name = outfile_prefix+str(num)+'.txt'
    outfiles.append(open(outfile_name, 'w'))

for file in outfiles:
    file.truncate(0)

count = 0
with open(infile, 'r') as f:
    line = f.readline()
    while line != '':
        for i, outfile in enumerate(outfiles):
            if count < sample_num[i]:
                outfile.write(line)
        count += 1
        line = line = f.readline()
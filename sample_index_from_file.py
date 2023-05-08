infile = 'testdata/toClustReads.txt'
outfile_prefix = 'testdata/toClustReads_'

outfiles = []

for file in outfiles:
    file.truncate(0)

count = 0
with open(infile, 'r') as f:
    line = f.readline()
    while line != '':
        tag, dna = line.strip().split(' ')
        dna = dna[0:16]
        for i, outfile in enumerate(outfiles):
            outfile.write(line)
        count += 1
        line = line = f.readline()
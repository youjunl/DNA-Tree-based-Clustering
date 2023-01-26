import random
import numpy as np

# Define the channel
pi, pd, ps = 0.02, 0.02, 0.02

# Define the number of clusters
nCluster = 4

# Define the size of each clusters
minSize = 6
maxSize = 16

# Target
def lfsr(state, mask):
  result = state
  nbits = mask.bit_length()-1
  while True:
    result = (result << 1)
    xor = result >> nbits
    if xor != 0:
        result ^= mask
    yield result

def screen_homopolymers(data, pattern):
  for w in pattern:
    if w in data:
      return 0
  return 1

num2dna = ['A', 'T', 'G', 'C']
dna2num = {'A':'00', 'T':'01', 'G':'10', 'C':'11'}
pattern = ['AAA', 'TTT', 'GGG', 'CCC']
def channel(txSeqs, pi, pd, ps):
    numSeqs = len(txSeqs)
    rxSeqs = []
    for ind, tx in enumerate(txSeqs):
        i, rx = 0, ''
        while i < len(tx):
            if randsel(pi): # Insertion
                rx += num2dna[np.random.randint(0, 4)]
                continue      
            if ~randsel(pd): # Deletion
                if randsel(ps): # Substitution
                    tmp = num2dna.copy()
                    tmp.remove(tx[i])
                    rx += tmp[np.random.randint(0, 3)]
                else:
                    rx += tx[i]
            i += 1
        rxSeqs.append(rx)
    return rxSeqs

def randsel(p):
    sel = np.random.rand() < p
    return sel

# PSNR Initialization
state = 0b001
mask = 0b100000000000000000000000011000101
psnr = lfsr(state, mask)

clInds = []
with open('toClust.txt', 'w') as f:
  f.truncate(0)
  for i in range(nCluster):
    nRepeat = random.randint(minSize, maxSize)
    txSeed = '{:032b}'.format(next(psnr))
    tx = ''
    for j in range(16):
      tx += num2dna[int(txSeed[j*2:j*2+2], 2)]
    while not screen_homopolymers(tx, pattern):
      txSeed = '{:032b}'.format(next(psnr))
      tx = ''
      for j in range(16):
        tx += num2dna[int(txSeed[j*2:j*2+2], 2)]

    rxs = channel([tx for _ in range(nRepeat)], pi, pd, ps)
    clInds.extend([i+1 for _ in range(nRepeat)]) 
    for rx in rxs:
      rx = rx[0:16]
      f.write(rx + '\n')

with open('clustInd.txt', 'w') as f:
  f.truncate(0)
  for ind in clInds:
    f.write(str(ind) + '\n')
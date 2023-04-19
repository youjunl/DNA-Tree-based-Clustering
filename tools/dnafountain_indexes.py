from tqdm import tqdm

pattern = ['AAA', 'TTT', 'GGG', 'CCC']

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


def screen_gc(data):
    gc = (data.count('G') + data.count('C') + 0.0)/ len(data)
    if (gc < 0.4) or (gc > 0.6):
        return 0
    return 1

if __name__ == '__main__':
    state = 0b101
    mask = 0b100000000000000000000000011000101
    psnr = lfsr(state, mask)
    intab = "0123"
    outtab = "ACGT"
    numSeed = 72000
    fout = open('dnafountain_indexes.txt', 'w')
    fout.truncate(0)
    for _ in tqdm(range(numSeed)):
        tx = ''
        txSeed = '{:032b}'.format(next(psnr))
        for j in range(16):
            tx += outtab[int(txSeed[j*2:j*2+2], 2)]
        while not (screen_homopolymers(tx, pattern) and screen_gc(tx)):
            tx = ''
            txSeed = '{:032b}'.format(next(psnr))
            for j in range(16):
                tx += outtab[int(txSeed[j * 2:j * 2 + 2], 2)]
        fout.write(tx+'\n')
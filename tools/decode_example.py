import numpy as np
from reedsolo import RSCodec

intab = "0123"
outtab = "ACGT"
dna_to_num = {'A':0,'C':1,'G':2,'T':3}
max_hamming = 2
def dna_to_int_array(dna_str):
    #convert a string like ACTCA to an array of ints like [10, 2, 4]
    s = ''
    for ch in dna_str:
        s += '{0:02b}'.format(int(dna_to_num[ch]))
    
    data = [int(s[t:t+8],2) for t in range(0,len(s), 8)]
    return data

f = open('toClustSmall.txt', 'r')
lines = f.readlines()
f.close

codec = RSCodec(2) # 2Bytes
cnt_corrected = 0
cnt_detected = 0
for content in lines:
    tag, string = content.strip().split(' ')
    if 'N' in string:
        continue
    dna = dna_to_int_array(string)
    try:
        #First is the decoded/repaired message
        #Second is the decoded message and error correction code (because both get repaired in reality)
        #Third is the position of the errors and erasures.
        data_corrected, _, _ = codec.decode(dna)

    except:
        continue #could not correct the code
    
    #we will encode the data again to evaluate the correctness of the decoding
    data_again = list(codec.encode(data_corrected)) #list is to convert byte array to int
    if np.count_nonzero(dna != data_again) > max_hamming: #measuring hamming distance between raw input and expected raw input
        #too many errors to correct in decoding
        cnt_detected += 1                 
        continue
    cnt_corrected += 1
print('%d out of %d are corrected'%(cnt_corrected, len(lines)))
print('%d out of %d are detected'%(cnt_detected, len(lines)))

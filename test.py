# -*- coding: utf-8 -*-
import random
import cupy
import chainer
from reedsolo import RSCodec
import numpy as np
import random
import matplotlib.pyplot as plt

"""# Functions"""

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

state = 0b101
mask = 0b100000000000000000000000011000101
psnr = lfsr(state, mask)

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

"""#LCS"""

def mutual_lcs(X, length):
    X = one_hot_encoding(X, {'A':0, 'T':1, 'G':2, 'C':3}, length)
    X = cupy.array(X)
    n = len(X)
    I = np.ravel(np.broadcast_to(np.arange(n), (n, n)).T)
    J = np.ravel(np.broadcast_to(np.arange(n), (n, n)))
    d = lcs(X[I], X[J])
    d = cupy.asnumpy(d).reshape(n, n)
    return d  

def mutual_edit(X, length):
    X = one_hot_encoding(X, {'A':0, 'T':1, 'G':2, 'C':3}, length)
    X = cupy.array(X)
    n = len(X)
    I = np.ravel(np.broadcast_to(np.arange(n), (n, n)).T)
    J = np.ravel(np.broadcast_to(np.arange(n), (n, n)))
    d = edit_distance(X[I], X[J])
    d = cupy.asnumpy(d).reshape(n, n)
    return d  

def all_lcs(y, X, length):
    X = one_hot_encoding(X, {'A':0, 'T':1, 'G':2, 'C':3}, length)
    X = cupy.array(X)
    y = one_hot_encoding([y], {'A':0, 'T':1, 'G':2, 'C':3}, length)
    y = cupy.array(y)
    n = len(X)
    d = lcs(y, X)
    d = cupy.asnumpy(d).reshape(n,)
    return d

def all_edit(y, X, length):
    X = one_hot_encoding(X, {'A':0, 'T':1, 'G':2, 'C':3}, length)
    X = cupy.array(X)
    y = one_hot_encoding([y], {'A':0, 'T':1, 'G':2, 'C':3}, length)
    n = len(X)
    y = cupy.array([y[0] for _ in range(n)])
    d = edit_distance(y, X)
    d = cupy.asnumpy(d).reshape(n,)
    return d

def one_hot_encoding(X, dict_alphabet, max_seq_length, smooth=1.):
    out = np.zeros((len(X), len(dict_alphabet), max_seq_length), dtype=np.float32)
    if smooth < 1:
        out[:] = (1 - smooth) / (len(dict_alphabet) - 1)
    for i, seq in enumerate(X):
        l = len(seq)
        for j, c in enumerate(seq):
            out[i, dict_alphabet[c], j] = smooth
        out[i, :, l:] = 0
    return out

def lcs(x1, x2, l1=None, l2=None):
    if len(x1.shape) == 3:
        kernel = cupy.ElementwiseKernel(
            'raw T x1, raw T x2, raw Z l1, raw Z l2, Z max_l, Z n_symbol',
            'raw T d',
            """
            int offset = i * (max_l + 1) * (max_l + 1);
            for(int j = 1; j < l1[i] + 1; j++){
                for(int k = 1; k < l2[i] + 1; k++){
                    int index = offset + j * (max_l + 1) + k;
                    T delta = 0;
                    int offset1 = i * max_l * n_symbol + j - 1;
                    int offset2 = i * max_l * n_symbol + k - 1;
                    for(int r = 0; r < n_symbol; r++)
                        delta += max(x1[offset1 + r * max_l] - x2[offset2 + r * max_l], 0.0);
                    if(delta > 0)
                        d[index] = max(d[index - (max_l + 1)], d[index - 1]);
                    else
                        d[index] = d[index - max_l - 2] + 1;
                }  
            } 
            """,
            name='lcs_kernel_W'
        )
        if l1 is None:
            l1 = cupy.sum((cupy.sum(x1, axis=1) > 0).astype(np.int32), axis=1)
            l2 = cupy.sum((cupy.sum(x2, axis=1) > 0).astype(np.int32), axis=1)
        d = cupy.zeros((x1.shape[0], x1.shape[2] + 1, x1.shape[2] + 1), dtype=np.float32)
        kernel(x1, x2, l1, l2, x1.shape[2], x1.shape[1], d, size=len(x1))
             
    else:
        kernel = cupy.ElementwiseKernel(
            'raw T x1, raw T x2, raw Z l1, raw Z l2, Z max_l',
            'raw T d',
            """
            int offset = i * (max_l + 1) * (max_l + 1);
            for(int j = 1; j < l1[i] + 1; j++){
                for(int k = 1; k < l2[i] + 1; k++){
                    int index = offset + j * (max_l + 1) + k;
                    T delta = 0;
                    if(x1[i * max_l + j - 1] != x2[i * max_l + k - 1])
                        d[index] = max(d[index - (max_l + 1)], d[index - 1]);
                    else
                        d[index] = d[index - max_l - 2] + 1;
                }  
            }
            """,
            name='lcs_kernel_W'
        )
        if l1 is None:
            l1 = cupy.sum((x1 > 0).astype(np.int32), axis=1)
            l2 = cupy.sum((x2 > 0).astype(np.int32), axis=1)
        if x1.shape[1] < 255:
            dtype = np.uint8
        else:
            if x1.shape[1] < 65535:
                dtype = np.uint16
            else:
                dtype = np.uint32
        d = cupy.zeros((x1.shape[0], x1.shape[1] + 1, x1.shape[1] + 1), dtype=dtype)
        kernel(x1.astype(d.dtype), x2.astype(d.dtype), l1.astype(d.dtype), l2.astype(d.dtype), x1.shape[1], d,
               size=len(x1))
    
    d = d[list(range(len(l1))), l1, l2]
    return d

def edit_distance(x1, x2, l1=None, l2=None, normolized=False):
    if len(x1.shape) == 3:
        kernel = cupy.ElementwiseKernel(
            'raw T x1, raw T x2, raw Z l1, raw Z l2, Z max_l, Z n_symbol',
            'raw T d',
            """
            int offset = i * (max_l + 1) * (max_l + 1);
            for(int j = 0; j < l1[i] + 1; j++)
                d[offset + j * (max_l + 1)] = j;
            for(int k = 0; k < l2[i] + 1; k++)
                d[offset + k] = k;
            for(int j = 1; j < l1[i] + 1; j++){
                for(int k = 1; k < l2[i] + 1; k++){
                    int index = offset + j * (max_l + 1) + k;
                    T delta = 0;
                    int offset1 = i * max_l * n_symbol + j - 1;
                    int offset2 = i * max_l * n_symbol + k - 1;
                    for(int r = 0; r < n_symbol; r++)
                        delta += max(x1[offset1 + r * max_l] - x2[offset2 + r * max_l], 0.0);
                    d[index] = min(d[index - (max_l + 1)] + 1, min(d[index - 1] + 1, d[index - max_l - 2] + delta));
                }  
            } 
            """,
            name='edit_distance_kernel_W'
        )
        if l1 is None:
            l1 = cupy.sum((cupy.sum(x1, axis=1) > 0).astype(np.int32), axis=1)
            l2 = cupy.sum((cupy.sum(x2, axis=1) > 0).astype(np.int32), axis=1)
        d = cupy.zeros((x1.shape[0], x1.shape[2] + 1, x1.shape[2] + 1), dtype=np.float32)
        kernel(x1, x2, l1, l2, x1.shape[2], x1.shape[1], d, size=len(x1))

    else:
        kernel = cupy.ElementwiseKernel(
            'raw T x1, raw T x2, raw Z l1, raw Z l2, Z max_l',
            'raw T d',
            """
            int offset = i * (max_l + 1) * (max_l + 1);
            for(int j = 0; j < l1[i] + 1; j++)
                d[offset + j * (max_l + 1)] = j;
            for(int k = 0; k < l2[i] + 1; k++)
                d[offset + k] = k;
            for(int j = 1; j < l1[i] + 1; j++){
                for(int k = 1; k < l2[i] + 1; k++){
                    int index = offset + j * (max_l + 1) + k;
                    T delta = 0;
                    if(x1[i * max_l + j - 1] != x2[i * max_l + k - 1])
                        delta = 1;
                    d[index] = min(d[index - (max_l + 1)] + 1, min(d[index - 1] + 1, d[index - max_l - 2] + delta));
                }  
            } 
            """,
            name='edit_distance_kernel'
        )
        if l1 is None:
            l1 = cupy.sum((x1 > 0).astype(np.int32), axis=1)
            l2 = cupy.sum((x2 > 0).astype(np.int32), axis=1)
        if x1.shape[1] < 255:
            dtype = np.uint8
        else:
            if x1.shape[1] < 65535:
                dtype = np.uint16
            else:
                dtype = np.uint32
        d = cupy.zeros((x1.shape[0], x1.shape[1] + 1, x1.shape[1] + 1), dtype=dtype)
        kernel(x1.astype(d.dtype), x2.astype(d.dtype), l1.astype(d.dtype), l2.astype(d.dtype), x1.shape[1], d,
               size=len(x1))

    d = d[list(range(len(l1))), l1, l2]
    if normolized:
        d = d - cupy.abs(l1 - l2)
    return d

"""# Create a dataset"""

pi, pd, ps = 0.01, 0.01, 0.01
nCluster = 200
clusterSize = [25, 40]
extraLen = 16
rscode = 2
max_hamming = 2
codec = RSCodec(rscode)
num2dna = ['A', 'T', 'G', 'C']
dna2num = {'A':'00', 'T':'01', 'G':'10', 'C':'11'}
dna_to_num = {'A':0,'C':1,'G':2,'T':3}
num_to_dna = {0:'A',1:'C',2:'G',3:'T'}
pattern = ['AAA', 'TTT', 'GGG', 'CCC']
def channel(txSeqs, pi, pd, ps):
    numSeqs = len(txSeqs)
    rxSeqs = []
    simerror = [0, 0, 0]
    for ind, tx in enumerate(txSeqs):
        i, rx = 0, ''
        while i < len(tx):
            if randsel(pi): # Insertion
                rx += num2dna[random.randint(0, 3)]
                simerror[0] += 1
                continue      
            if not randsel(pd): # Deletion
                if randsel(ps): # Substitution
                    tmp = num2dna.copy()
                    tmp.remove(tx[i])
                    rx += tmp[random.randint(0, 2)]
                    simerror[2] += 1
                else:
                    rx += tx[i]
            else:
              simerror[1] += 1
            i += 1
        rxSeqs.append(rx)
    return rxSeqs, simerror

def randsel(p):
    sel = np.random.rand() < p
    return sel

def dna_to_int_array(dna_str):
    #convert a string like ACTCA to an array of ints like [10, 2, 4]
    s = ''
    for ch in dna_str:
        s += '{0:02b}'.format(int(dna_to_num[ch]))

    data = [int(s[t:t+8],2) for t in range(0,len(s), 8)]
    return data

def int_array_to_dna(data):
    s = ''
    bin_data = ''
    for num in data:
        bin_data += '{0:08b}'.format(num)
    for i in range(0, len(bin_data), 2):
        s += num_to_dna[int(bin_data[i:i+2], 2)]
    return s

psnr = lfsr(state, mask) # PSNR Initialization
clInds = []
simerror = np.array([0, 0, 0])
with open('seeds.txt', 'w') as f:
    f.truncate(0)
    for i in range(nCluster):
        nRepeat = random.randint(clusterSize[0], clusterSize[1])
        txSeed = '{:032b}'.format(next(psnr))
        tx = ''
        # for j in range(16):
        #     tx += num2dna[int(txSeed[j*2:j*2+2], 2)]
        # for j in range(extraLen):
        #     tx += num2dna[random.randint(0, 3)]
        # while not (screen_homopolymers(tx, pattern) and screen_gc(tx)):
        #     txSeed = '{:032b}'.format(next(psnr))
        #     tx = ''
        #     for j in range(16):
        #         tx += num2dna[int(txSeed[j*2:j*2+2], 2)]
        #     for j in range(extraLen):
        #         tx += num2dna[random.randint(0, 3)]
        for j in range(16+extraLen):
            tx += num2dna[random.randint(0, 3)]
        if rscode > 0:
            txInt = dna_to_int_array(tx)
            txCode = codec.encode(txInt)
            tx = int_array_to_dna(txCode)
        rxs, e = channel([tx for _ in range(nRepeat)], pi, pd, ps)
        simerror = simerror + np.array(e)
        clInds.extend([i+1 for _ in range(nRepeat)]) 
        for rx in rxs:
            rx = rx[0:16+extraLen]
            # rxSeedStr = ''
            # for c in rx:
            #   rxSeedStr += dna2num[c]
            f.write(rx + '\n')

with open('clustInd.txt', 'w') as f:
    f.truncate(0)
    for ind in clInds:
        f.write(str(ind) + '\n')

"""# Data preprocessing"""

# Data preprocessing
data = np.genfromtxt('seeds.txt' ,delimiter='\n',dtype=str)
clInds = np.genfromtxt('clustInd.txt' ,delimiter='\n',dtype='uint64')
nRead = len(data)
print('Import %d PRNs'%nRead)
# Random Shuffling
inds = np.arange(nRead)
Shuffling = False
if Shuffling:
  np.random.shuffle(inds)
  data = data[inds]
  clInds = clInds[inds]
dnaData = [s for s in data]
# dnaData = []
# for bin in binData:
#   seq = ''
#   for i in range(16):
#     seq += num2dna[int(bin[i*2:i*2+2], 2)]
#   dnaData.append(seq)

with open('refInd.txt', 'w') as f:
  f.truncate(0)
  for ind in clInds:
    f.write(str(ind) + '\n')

"""# Data Structure (Improved with local greedy strategy)"""

"""Tree Structure Module
This module provides retrieval algorithms based on tree structures, 
including fuzzy search algorithms that allow for horizontal drift.
"""

class Trie:
    """Tree Structure Class
    This class is used to initialize a tree structure without the 
    need to additionally specify the depth of the tree.
    Attributes:
        dna_dict: dict,The elements contained in the sequence.
        node_nums: int,The number of elements contained in the sequence.
        children: int,Number of tree branches.
        isEnd: int,The value to determine if the tree is terminated.
    
    """

    def __init__(self):
        self.dna_dict = {"A":0,"T":1,"G":2,"C":3}
        self.node_nums = len(self.dna_dict)
        self.children = [None] * self.node_nums
        #self.freqs = [0] * self.node_nums
        self.isEnd = False
        self.maxOptimDepth = 3
    #Node retrieval function without drift.
    def searchPrefix(self, prefix: str) -> "Trie":
        dict=self.dna_dict
        node = self
        for ch in prefix:
            ch = dict[ch]
            if not node.children[ch]:
                return None
            node = node.children[ch]
        return node.isEnd

    def insert(self, word: str,label:str) -> None:
        """Add a branch to the tree.
        
        Args:
            word: str,The sequence added to the tree.
            label: str,Sequence of labels.
        """
        node = self
        dict=self.dna_dict
        for ch in word:
            ch = dict[ch]
            if not node.children[ch]:
                node.children[ch] = Trie()
            node = node.children[ch]
        node.isEnd = label

    
    def delete(self,word:str):
        """Deletes a branch from the tree.
        
        Args:
            word: str, Sequences deleted from the tree.
        """
        node = self
        dict=self.dna_dict
        for ch in word:
            ch = dict[ch]
            if not node.children[ch]:
                node.children[ch] = Trie()
            node = node.children[ch]
        node.isEnd = False

    def fuzz_align(self, dna):
        """Horizontal drift function
        Args:
            word: Sequence of fuzzy retrieval.
        
        return:
            Returns a list with the positions that need to be drifted 
            laterally and the nodes that can be drifted laterally.
        
        """
        word, e = dna
        currentErr = sum(e)
        node = self 
        dict = self.dna_dict
        num=0
        tmp_list=[]
        sub_list=[]
        ins_list=[]
        del_list=[]
        len_=len(word)

        for pos, ch in enumerate(word) :

            ch = dict[ch]
            if not node.children[ch]:
                
                for i in range(self.node_nums):
                    if not node.children[i]:
                        pass
                    else:
                        tmp_list.append(i)
                
                maxOptimDepth = self.maxOptimDepth
                traverseNum = min(len_-num-1, maxOptimDepth)
                for k in tmp_list:
                    #Deletion
                    depth = 0
                    tmp = node.children[k]
                    while depth < traverseNum:
                        if not tmp.children[dict[word[num+depth]]]:
                            break
                        tmp = tmp.children[dict[word[num+depth]]]
                        depth += 1

                    if depth == traverseNum:
                        del_list.append(k)

                    #Insertion
                    depth = 0
                    tmp = node
                    while depth < traverseNum-1:
                        if not tmp.children[dict[word[num+depth+1]]]:
                            break
                        tmp = tmp.children[dict[word[num+depth+1]]]
                        depth += 1
                    
                    if depth == traverseNum-1:
                        ins_list.append(k)                     

                    #Substitution
                    depth = 0
                    tmp = node.children[k]
                    while depth < traverseNum:                                                
                        if not tmp.children[dict[word[num+depth+1]]]:
                            break
                        tmp = tmp.children[dict[word[num+depth+1]]]
                        depth += 1

                    if depth == traverseNum:
                        sub_list.append(k)

                return [num, sub_list, ins_list, del_list]

            else:
                node = node.children[ch]
            num = num + 1

        if type(node.isEnd) == int:
          return node.isEnd

        elif type(node.isEnd) == bool:
          tmp_list = []
          for i in range(self.node_nums):
              if not node.children[i]:
                  pass
              else:
                  tmp_list.append(i)
          return [num, [], [], tmp_list]

    def fuzz_fin(self,word,max_value):
        """Fuzzy search with horizontal drift.
        Args:
            word: str,Sequence of search
            max_value: int,The maximum number of horizontal drifts.
        return:
            Returns a list, the first element of which is the index of the final matched 
            core sequence, and the second element is the number of horizontal drifts.
        """

        queue=[[word,[0, 0, 0]]]
        
        fin_list=["",[1000, 1000, 1000]]
        num2dnaDict = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
        adjacent = []
        while True :
            if queue == [] or fin_list[1] == 0 :
                break
            dna = queue.pop(0)
            errorSum = sum(dna[1])
            if errorSum > max_value :
              # Current number of error is larger than the threshold
                continue

            result = self.fuzz_align(dna)
            if type(result) == int :
                adjacent.append([result, dna[1]])
                if errorSum < sum(fin_list[1]) :
                    fin_list=[result, dna[1]]

            elif result[0] == len(dna[0])-1:
                for i in range(len(result[1])):
                    chNum = result[1][i]
                    k = dna[0][:result[0]]+num2dnaDict[chNum]
                    dna[1][1] += 1
                    queue.append([k,dna[1]])
            else:
                # Substitution Fix
                for i in range(len(result[1])):
                    chNum = result[1][i]
                    k = dna[0][:result[0]]+num2dnaDict[chNum]+dna[0][result[0]-len(dna[0])+1:]
                    errorList = [dna[1][0], dna[1][1], dna[1][2]+1]
                    queue.append([k,errorList])

                # Insertion Fix
                for i in range(len(result[2])):
                    k = dna[0][:result[0]]+dna[0][result[0]-len(dna[0])+1:]+'A'
                    errorList = [dna[1][0]+1, dna[1][1], dna[1][2]]
                    queue.append([k,errorList])

                # Deletion Fix
                for i in range(len(result[3])):
                    chNum = result[3][i]
                    k = dna[0][:result[0]]+num2dnaDict[chNum]+dna[0][result[0]-len(dna[0]):]
                    k = k[:16]
                    errorList = [dna[1][0], dna[1][1]+1, dna[1][2]]
                    queue.append([k,errorList])
        fin_list.append(adjacent)
        return fin_list

"""# Clustering"""

from pyparsing.helpers import List
import time
from collections import defaultdict
def clust(tree, dnaData, indexBegin = 0):
  test_num = 0
  tree_threshold = 6
  h_drift = 6
  v_drift = 4
  indexList = []
  dna_number = indexBegin
  timer = time.time()
  error = np.array([0, 0, 0])
  freqs = [0 for _ in range(len(dnaData))]
  adjacent_list = []
  core_set = dict()
  for i, seq in enumerate(dnaData):
    test_num += 1
    align = tree.fuzz_fin(seq, tree_threshold)
    if sum(align[1]) < h_drift:
      indexList.append((clInds[indexBegin+i], align[0]))
      error = error + np.array(align[1])
    else:
      dna_number += 1
      tree.insert(seq, dna_number)
      core_set[dna_number] = seq
      indexList.append((clInds[indexBegin+i], dna_number)) # Fix Clover's problem
  # Output Result
  t = (time.time() - timer)
  sortedList = sorted(indexList, key = lambda k: k[0])
  prob = error / len(dnaData) / len(dnaData[0])
  return indexList

def clustMP(tree, dnaData, indexList, indexBegin = 0):
  test_num = 0
  tree_threshold = 6
  h_drift = 6
  v_drift = 4
  dna_number = indexBegin
  timer = time.time()
  error = np.array([0, 0, 0])
  freqs = [0 for _ in range(len(dnaData))]
  adjacent_list = []
  core_set = dict()
  for i, seq in enumerate(dnaData):
    test_num += 1
    align = tree.fuzz_fin(seq, tree_threshold)
    if sum(align[1]) < h_drift:
      indexList[i+indexBegin] = (clInds[indexBegin+i], align[0])
      error = error + np.array(align[1])
    else:
      dna_number += 1
      tree.insert(seq, dna_number)
      core_set[dna_number] = seq
      indexList[i+indexBegin] = (clInds[indexBegin+i], dna_number) # Fix Clover's problem

def computeAccuracy(indexList):
    f = open('refInd.txt', 'r').readlines()
    refind = [int(t.strip()) for t in f]
    algind = indexList
    # caculate freqs
    clustNum = defaultdict(int)
    results = indexList
    for ind in refind:
        clustNum[ind] += 1

    results = sorted(results, key=lambda k: k[1])
    nClust = len(clustNum.keys())

    maxClustNum = results[-1][1] if len(results) > 0 else 0
    clusters = [[] for _ in range(maxClustNum)]
    for pair in results:
        clusters[pair[1]-1].append(pair[0])

    clusters = [c for c in clusters if c != []]
    WrongClusterNum = 0
    score = [0] * max(clustNum.keys())
    for cluster in clusters:
        # Check if all the tags in a clusters are the same.
        tags = set(cluster)
        if len(tags) > 1:
            WrongClusterNum += 1
            # This cluster is invalid
            continue
        # Maximize the size of the clusters.
        tag = int(cluster[0])
        score[tag-1] = max(score[tag-1], len(cluster))

    # Compute accuracy
    gamma = [i*0.05 for i in range(0, 21)]
    acc = [0] * len(gamma)
  
    for i, g in enumerate(gamma):
        cnt = 0
        for tag in clustNum.keys():
            if score[tag-1] / clustNum[tag] >= g:
                cnt += 1
        acc[i] = cnt / nClust

    print('Incorrect cluster: %d'%WrongClusterNum)
    return gamma, acc

"""# SSC"""

def dna_to_seed(dna_str):
    # revert seed from the DNA sequence
    s = ''
    for ch in dna_str:
        s += '{0:02b}'.format(int(dna_to_num[ch]))    
    return int(s, 2)

def seed_to_dna(seed):
    bin_str = '{:032b}'.format(seed)
    dna_str = ''
    for t in range(0, len(bin_str), 2):
        dna_str += num_to_dna[int(bin_str[t:t+2],2)]
    return dna_str

def SSC(indexList, dnaData):
    begin = time.time()
    maxClusterIndex = max(indexList, key=lambda k: k[1])
    tags = [index[1] for index in indexList]
    clusters = [[] for _ in range(maxClusterIndex[1])]
    for i, read in enumerate(dnaData):
        clustInd, tag = indexList[i]
        clusters[tag-1].append(read)
    cnt_corrected = 0
    cnt_detected = 0
    seeds = []
    seedMap = dict()
    clusterMap = dict()
    for i, cluster in enumerate(clusters):
        # Skip empty clusters
        if len(cluster) == 0:
            continue
        
        # Record the frequencies of seeds
        freq = dict()
        for read in cluster:
            seed = dna_to_seed(read[:16])
            if seed not in freq.keys():
                freq[seed] = 1
            else:
                freq[seed] += 1
        
        # Majority selection for the center seed
        maxSeed = -1
        maxFreq = 0
        for k in freq.keys():
            if freq[k] > maxFreq:
                maxFreq = freq[k]
                maxSeed = k
        if maxSeed == -1:
            continue
        
        # Create a mapping
        clusterInd = i + 1
        if maxSeed not in seeds:
            seeds.append(maxSeed)
            seedMap[maxSeed] = clusterInd # For merging the other clusters
            clusterMap[clusterInd] = clusterInd
        else:
            clusterMap[clusterInd] = seedMap[maxSeed]
            
    for tag in set(tags):
        if tag not in clusterMap.keys():
            clusterMap[tag] = tag # Map to itself
        
    newIndexList = []
    for clustInd, tag in indexList:
        newIndexList.append((clustInd, clusterMap[tag]))
    tags = [index[1] for index in newIndexList]
    return newIndexList

def SSCRS(indexList, dnaData):
    begin = time.time()
    maxClusterIndex = max(indexList, key=lambda k: k[1])
    tags = [index[1] for index in indexList]
    clusters = [[] for _ in range(maxClusterIndex[1])]
    for i, read in enumerate(dnaData):
        clustInd, tag = indexList[i]
        clusters[tag-1].append(read)
    cnt_corrected = 0
    cnt_detected = 0
    seeds = []
    seedMap = dict()
    clusterMap = dict()
    for i, cluster in enumerate(clusters):
        # Skip empty clusters
        if len(cluster) == 0:
            continue
        
        # Record the frequencies of seeds
        freq = dict()
        for read in cluster:
            if 'N' in read:
                continue
            dna = dna_to_int_array(read)
            detected = False
            try:
                #First is the decoded/repaired message
                #Second is the decoded message and error correction code (because both get repaired in reality)
                #Third is the position of the errors and erasures.
                data_corrected, _, _ = codec.decode(dna)
                detected = True
                
            except:
                detected = False #could not correct the code

            if detected:
                #we will encode the data again to evaluate the correctness of the decoding
                data_again = list(codec.encode(data_corrected)) #list is to convert byte array to int
                if np.count_nonzero(dna != data_again) > max_hamming: #measuring hamming distance between raw input and expected raw input
                    #too many errors to correct in decoding
                    cnt_detected += 1                 
                    seed = dna_to_seed(read[:16])
                
                else:
                    cnt_corrected += 1
                    dna_str_corrected = int_array_to_dna(data_again)
                    seed = dna_to_seed(dna_str_corrected[:16])
                

            else:
                seed = dna_to_seed(read[:16])

            if seed not in freq.keys():
                freq[seed] = 1
            else:
                freq[seed] += 1
        
        # Majority selection for the center seed
        maxSeed = -1
        maxFreq = 0
        for k in freq.keys():
            if freq[k] > maxFreq:
                maxFreq = freq[k]
                maxSeed = k
        if maxSeed == -1:
            continue
        
        # Create a mapping
        clusterInd = i + 1
        if maxSeed not in seeds:
            seeds.append(maxSeed)
            seedMap[maxSeed] = clusterInd # For merging the other clusters
            clusterMap[clusterInd] = clusterInd
        else:
            clusterMap[clusterInd] = seedMap[maxSeed]
            
    for tag in set(tags):
        if tag not in clusterMap.keys():
            clusterMap[tag] = tag # Map to itself

    newIndexList = []
    for clustInd, tag in indexList:
        newIndexList.append((clustInd, clusterMap[tag]))
    tags = [index[1] for index in newIndexList]
    return newIndexList

def SSCR(indexList, dnaData, tree_threshold, tree_depth, filter=False):
    maxClusterIndex = max(indexList, key=lambda k: k[1])
    tags = [index[1] for index in indexList]
    clusters = [[] for _ in range(maxClusterIndex[1])]
    for i, read in enumerate(dnaData):
        clustInd, tag = indexList[i]
        clusters[tag-1].append(read)
    cnt_corrected = 0
    cnt_detected = 0
    seeds = []
    seedMap = dict()
    clusterMap = dict()
    core = []
    tree = Trie()
    clustNum = 0
    for i, cluster in enumerate(clusters):
        # Skip empty clusters
        if len(cluster) == 0:
            continue
        clustNum += 1
        # Record the frequencies of seeds
        freq = dict()
        for read in cluster:
            prefix = read[:tree_depth]
            seed = prefix
            if seed not in freq.keys():
                freq[seed] = 1
            else:
                freq[seed] += 1
        
        # Majority selection for the center seed
        maxSeed = -1
        maxFreq = 0
        for k in freq.keys():
            if freq[k] > maxFreq:
                maxFreq = freq[k]
                maxSeed = k
        if maxSeed == -1:
            continue
        
        # Create a mapping
        clusterInd = i + 1

        # Clustering with Tree Structure
        seedStr = maxSeed
        align = tree.fuzz_fin(seedStr, tree_threshold)
        if sum(align[1]) < tree_threshold:
            # print('%s Merge To %s Err: %d'%(seedStr[-tree_depth:], readMap[align[0]], sum(align[1])))
            clusterMap[clusterInd] = align[0] # Merging
            continue

        if core and filter:
            # Prefiltering with LCS
            distance = all_edit(seedStr, core, tree_depth)
            ind = np.argmin(distance)
            if distance[ind] < tree_threshold:
                clusterMap[clusterInd] = clusterMap[seedMap[core[ind]]]
                continue

        tree.insert(seedStr, clusterInd)
        seedMap[seedStr] = clusterInd # For merging the other clusters
        clusterMap[clusterInd] = clusterInd
        core.append(seedStr)
            
    for tag in set(tags):
        if tag not in clusterMap.keys():
            clusterMap[tag] = tag # Map to itself

    newIndexList = []
    for clustInd, tag in indexList:
        newIndexList.append((clustInd, clusterMap[tag]))
    tags = [index[1] for index in newIndexList]
    return newIndexList

# def SSCML(indexList, dnaData):
#     begin = time.time()
#     maxClusterIndex = max(indexList, key=lambda k: k[1])
#     tags = [index[1] for index in indexList]
#     clusters = [[] for _ in range(maxClusterIndex[1])]
#     for i, read in enumerate(dnaData):
#         clustInd, tag = indexList[i]
#         clusters[tag-1].append(read)
#     cnt_corrected = 0
#     cnt_detected = 0
#     seeds = []
#     seedMap = dict()
#     clusterMap = dict()
#     for i, cluster in enumerate(clusters):
#         # Skip empty clusters
#         if len(cluster) == 0:
#             continue
        
#         # Record the frequencies of seeds
#         freq = dict()
#         for read in cluster:
#             seed = dna_to_seed(read[:16])
#             if seed not in freq.keys():
#                 freq[seed] = 1
#             else:
#                 freq[seed] += 1
        
#         # Majority selection for the center seed
#         maxSeed = -1
#         maxFreq = 0
#         for k in freq.keys():
#             if freq[k] > maxFreq:
#                 maxFreq = freq[k]
#                 maxSeed = k
#         if maxSeed == -1:
#             # Select random
#             maxSeed = cluster[random.randint(0, len(cluster)-1)][:16]
        
#         # Create a mapping
#         clusterInd = i + 1
#         if maxSeed not in seeds:
#             seeds.append(maxSeed)
#             seedMap[maxSeed] = clusterInd # For merging the other clusters
#             clusterMap[clusterInd] = clusterInd
#         else:
#             clusterMap[clusterInd] = seedMap[maxSeed]
            
#     for tag in set(tags):
#         if tag not in clusterMap.keys():
#             clusterMap[tag] = tag # Map to itself

#     # Clustering with Kmeans
#     alphabet = np.array(['A', 'C', 'G', 'T'])
#     alg = SoftSeqKmeans(n_centroid = nCluster, centroid_length = 16, alphabet=alphabet)
#     data = []
#     dataClusterInd = []
#     for seed in seeds:
#         dnaSeed = seed_to_dna(seed)
#         data.append(dnaSeed)
#         dataClusterInd.append(seedMap[seed])
    
#     lcurve = alg.fit(np.array(data), n_iter=20)
#     labels = alg.transform(data)

#     print('Number of clusters before KMeans: %d'%len(data))
#     print('Number of clusters after KMeans: %d'%len(set(labels)))

#     # Create mapping
#     MLMap = dict()
#     for i, label in enumerate(labels):
#         if label not in MLMap.keys():
#             MLMap[label] = dataClusterInd[i] # Set a root for merging

#         else:
#             # Update merging
#             clusterMap[dataClusterInd[i]] = MLMap[label]

#     newIndexList = []
#     for clustInd, tag in indexList:
#         newIndexList.append((clustInd, clusterMap[tag]))

#     tags = [index[1] for index in newIndexList]
#     return newIndexList

def MSCR(indexList, dnaData, tree_threshold, start_tree_depth, iter = 3, reduceSize = 2, filter = False):
    indexListMSCR = indexList
    tree_depth = start_tree_depth
    for i in range(iter):
        indexListMSCR = SSCR(indexListMSCR, dnaData, tree_threshold, tree_depth, filter=filter)
        tree_depth -= reduceSize
    return indexListMSCR

"""# Clustering

"""

threshold = 3
print('--------------------------------')
begin = time.time()
tree = Trie()
indexList = clust(tree, dnaData)
print('Time: %f'%(time.time()-begin))
gamma, acc = computeAccuracy(indexList)
print('--------------------------------')
indexListSSC = SSC(indexList, dnaData)
gamma, accSSC = computeAccuracy(indexListSSC)
print('--------------------------------')
indexListSSCRS = SSCRS(indexList, dnaData)
gamma, accSSCRS = computeAccuracy(indexListSSCRS)
# indexListSSCML = SSCML(indexList, dnaData)
# gamma, accSSCML = computeAccuracy(indexListSSCML)
print('--------------------------------')
indexListSSCR = SSCR(indexList, dnaData, threshold, 14)
gamma, accSSCR = computeAccuracy(indexListSSCR)
print('--------------------------------')
indexListSSCRF = SSCR(indexList, dnaData, threshold, 14, filter=True)
gamma, accSSCRF = computeAccuracy(indexListSSCRF)
print('--------------------------------')

# MSCR
begin = time.time()
indexListMSCR = MSCR(indexList, dnaData, threshold, 14, iter=2, reduceSize=2, filter=False)
print('Time: %f'%(time.time()-begin))
_, accMSCR = computeAccuracy(indexListMSCR)
print('--------------------------------')

# MSCRF
begin = time.time()
indexListMSCRF = MSCR(indexList, dnaData, threshold, 14, iter=2, reduceSize=2, filter=True)
print('Time: %f'%(time.time()-begin))
_, accMSCRF = computeAccuracy(indexListMSCRF)
print('--------------------------------')

# MSCR Multi Process
import math
import threading

begin = time.time()
#split the arr into N chunks
def chunks(arr, m):
    n = int(math.ceil(len(arr) / float(m)))
    return [arr[i:i + n] for i in range(0, len(arr), n)], [i for i in range(0, len(arr), n)]

nProcess = 8
splitData, indexBegin = chunks(dnaData, nProcess)
indexListMSCRMP = [[] for _ in range(len(dnaData))]
Process = []
for i in range(nProcess):
    tree = Trie()
    Process.append(threading.Thread(target=clustMP, args=(tree, splitData[i], indexListMSCRMP, indexBegin[i])))
    Process[i].start()
for i in range(nProcess):
    Process[i].join()
indexListMSCRMP = MSCR(indexListMSCRMP, dnaData, threshold, 14, iter=2, reduceSize=2, filter=False)
print('Time: %f'%(time.time()-begin))
_, accMSCRMP = computeAccuracy(indexListMSCRMP)

# sortedList = sorted(indexListMSCR, key = lambda k: k[0])
with open('result.txt', 'w') as f:
    f.truncate(0)
    for t in indexListSSCRF:
        f.write(str(t[0]) + ' ' + str(t[1]) + '\n')

plt.figure(dpi=250)
plt.plot(gamma, acc, '-', label='Single-stage')
plt.plot(gamma, accSSC, '--^', label='Two-stage')
plt.plot(gamma, accSSCRS, '--o', label='Two-stage RS')
# plt.plot(gamma, accSSCML, '--*')
plt.plot(gamma, accSSCR, '--*', label='Two-stage Repeat')
plt.plot(gamma, accSSCRF, '--s', label='Two-stage Repeat with Filter')
plt.plot(gamma, accMSCR, '--^', label='Multi-stage Repeat')
plt.plot(gamma, accMSCRF, '--o', label='Multi-stage Repeat with Filter')
plt.plot(gamma, accMSCRMP, '--*', label='Multi-stage Repeat MP')
plt.legend()

plt.show()

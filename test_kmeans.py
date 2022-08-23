from kmeans_pytorch import KMeans
import numpy as np
from channel import *
import torch
import matplotlib.pyplot as plt
num_x = 3000
pi = 0.2
pd = 0.2
ps = 0.2
alphabet = 'ATGC'
sequences = ['ATTGCA', 'GCTACC',
             'GCGAAT']
# Repeat sequences n times
tx_strands = sequences * num_x
# Get output from IDS channel
rx_strands = ids_channel(tx_strands, pi, pd, ps)

def onehot(string, alphabet, max_len):
    n = len(string)
    m = len(alphabet)
    result = np.zeros((max_len, m))
    for i in range(n):
        for j in range(m):
            if string[i] == alphabet[j]:
                result[i, j] = 1
    return result

# K-Means params
mini_batch = 50
n_iter = 100
n_centroid = len(sequences)
centroid_length = len(sequences[0])
L = np.array([len(seq) for seq in rx_strands])
max_len = np.max(L)
tau = -2
# Initialize k centroid points randomly
random_indexes = np.where(L == centroid_length)[0]
random_indexes = random_indexes[np.random.choice(len(random_indexes), n_centroid)]

init = []
for i in random_indexes:
    init.append(rx_strands[i])
centroid_str = np.unique(init)
# Encoding with one-hot coding
centroid = []
for c in centroid_str:
    centroid.append(onehot(c, alphabet, max_len))

X = []
for seq in rx_strands:
    X.append(onehot(seq, alphabet, max_len))
X = torch.Tensor(np.array(X))
print(X.shape)
model = KMeans(n_clusters = 3, max_iter=100, minibatch=50, verbose = 2)
model.fit_predict(X)
max_sim_v, max_sim_i = model.predict(X)


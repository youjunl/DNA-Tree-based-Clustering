# from collections import defaultdict
# from ossaudiodev import control_labels
# from pickletools import optimize
# from turtle import forward
import numpy as np
from distance import *
from channel import *
import torch
num_x = 400
pi = 0.2
pd = 0.2
ps = 0.2
alphabet = 'ATGC'
sequences = ['AAAATTTT', 'CCCCCCCC',
             'GGGGGGGG', 'TTTTAAAA']
# Repeat sequences n times
tx_strands = sequences * num_x
# Get output from IDS channel
rx_strands = ids_channel(tx_strands, pi, pd, ps)
X = rx_strands

def onehot(string, alphabet, max_len):
    n = len(string)
    m = len(alphabet)
    result = np.zeros((max_len, m))
    for i in range(n):
        for j in range(m):
            if string[i] == alphabet[j]:
                result[i, j] = 1
    return result

def softEditDistance(X1: np.array, X2: np.array, tau):
    #print(type(X2))
    s_len = X1.shape[0]
    t_len = X2.shape[0]
    alpha = np.zeros((s_len + 1, t_len + 1))
    beta = np.zeros((s_len + 1, t_len + 1))
    # Initialization
    for i in range(0, s_len + 1):
        alpha[i, 0] = i*np.exp(tau*i)
        beta[i, 0] = i*np.exp(tau*i)
    for j in range(0, t_len + 1):
        alpha[0, j] = j*np.exp(tau*j)
        beta[0, j] = j*np.exp(tau*j)
    for j in range(t_len):
        for i in range(s_len):
            soft_sigma = 0.5*sum(abs(X1[i, :]-X2[j, :]))
            alpha[i+1, j+1] = np.exp(tau)*(alpha[i, j+1] + alpha[i+1, j] + beta[i, j+1] + beta[i+1, j]) + np.exp(
                tau*soft_sigma)*(alpha[i, j] + beta[i, j]*soft_sigma) - np.exp(2*tau)*(alpha[i, j] + 2*beta[i, j])
            beta[i+1, j+1] = np.exp(tau)*(beta[i, j+1] + beta[i+1, j]) + \
                beta[i, j]*(np.exp(tau*soft_sigma)-np.exp(2*tau))
    return alpha[-1, -1]/beta[-1, -1]

# K-Means
mini_batch = 50
n_iter = 10
n_centroid = len(sequences)
centroid_length = len(sequences[0])

L = np.array([len(seq) for seq in rx_strands])
max_len = np.max(L)
alphaLen = len(alphabet)
X = [onehot(sample, alphabet, max_len) for sample in rx_strands]

tau = -4
# Initialize k centroid points randomly
random_indexes = np.where(L == centroid_length)[0]
random_indexes = random_indexes[np.random.choice(len(random_indexes), n_centroid)]

init = []
for i in random_indexes:
    init.append(X[i])
#init = np.unique(init)

# Encoding with one-hot coding
centroid = []
for c in init:
    centroid.append(c)
print('Center _______________________')
for c in centroid:
    seq = ''
    for ind in np.argmax(c, axis=1):
        seq += alphabet[ind]
    print(seq)    

for i in range(n_iter):
    # Ramdomly pick M sample from X
    indexes = np.random.choice(len(X), mini_batch)
    x = []
    for index_i in indexes:
        x.append(X[index_i])
    # Compute the soft edit distance of a mini-batch x
    centroid_indexes = np.zeros((mini_batch,))
    clusters = [[c] for c in centroid]
    min_sed = np.zeros((mini_batch,))
    for sample_index, sample in enumerate(x):
        sed = np.zeros(len(centroid))
        for centroid_index, c in enumerate(centroid):
           
            sed[centroid_index] = softEditDistance(sample, c, tau)
        # Assign the sequence to the nearest centroid
        est_centroid_index = int(np.argmin(sed))
        # Save the soft edit distance of the point
        min_sed[sample_index] = sed.min()
        # Save the index of nearest centroid
        centroid_indexes[sample_index] = est_centroid_index
        # Save the samples in the cluster
        #print(clusters)
        clusters[est_centroid_index].append(sample)

    # For each cluster, update the centroid
    for saved_centroid_index, cluster in enumerate(clusters):
        mean_centroid = np.mean(cluster, axis=0)
        #learning_rate = 1/len(cluster)
        indexes = mean_centroid.argmax(axis=1)
        mask = np.zeros((max_len, alphaLen))
        for mask_ind in range(mask.shape[0]):
            mask[mask_ind, indexes[mask_ind]] = 1
        centroid[saved_centroid_index] = mask

    # Calculate Loss for this iteration
    Loss = 0
    for sample_index, sample in enumerate(X):
        sed = np.zeros(len(centroid))
        for centroid_index, c in enumerate(centroid):
            sed[centroid_index] = softEditDistance(sample, c, tau)
        Loss += sed.min()
    Loss /= len(X)
    print('iter %d, Loss %2f'%(i, Loss))
    # Print centrod
    print('Center _______________________')
    for c in centroid:
        seq = ''
        for ind in np.argmax(c, axis=1):
            seq += alphabet[ind]
        print(seq)    

print(centroid)

from traceback import print_tb
import matplotlib.pyplot as plt
import numpy as np
from channel import ids_channel
from sklearn.manifold import TSNE
from itertools import cycle, islice

def onehot(string, alphabet, max_len):
    n = len(string)
    m = len(alphabet)
    result = np.zeros((max_len, m))
    for i in range(n):
        for j in range(m):
            if string[i] == alphabet[j]:
                result[i, j] = 1
    return result

# Test soft edit distance
def ed(x, center, tau=-4):
    """
        Compute soft edit distance of one-hot encoded sequences
        Parameters:
        x: torch.Tensor, shape: [n, maxLen, alphabetLen]
        center: torch.Tensor, shape: [n, maxLen, alphabetLen]
        tau: float
    """
    maxLen = x.shape[1]
    nSample = x.shape[0]
    nCenter = center.shape[0]
    alpha = np.zeros((nSample, nCenter, maxLen + 1, maxLen + 1))
    beta = np.zeros((nSample, nCenter, maxLen + 1, maxLen + 1))
    # Initialization
    for s_i in range(nSample):
        for c_i in range(nCenter):
            for i in range(0, maxLen + 1):
                alpha[s_i, c_i, i, 0] = i*np.exp(tau*i)
                beta[s_i, c_i, i, 0] = i*np.exp(tau*i)
            for j in range(0, maxLen + 1):
                alpha[s_i, c_i, 0, j] = j*np.exp(tau*j)
                beta[s_i, c_i, 0, j] = j*np.exp(tau*j)
            for j in range(maxLen):
                for i in range(maxLen):
                    soft_sigma = 0.5*sum(abs(x[s_i, i, :]-center[c_i, j, :]))
                    alpha[s_i, c_i, i+1, j+1] = np.exp(tau)*(alpha[s_i, c_i, i, j+1] + alpha[s_i, c_i, i+1, j] + beta[s_i, c_i, i, j+1] + beta[s_i, c_i, i+1, j]) + np.exp(
                        tau*soft_sigma)*(alpha[s_i, c_i, i, j] + beta[s_i, c_i, i, j]*soft_sigma) - np.exp(2*tau)*(alpha[s_i, c_i, i, j] + 2*beta[s_i, c_i, i, j])
                    beta[s_i, c_i, i+1, j+1] = np.exp(tau)*(beta[s_i, c_i, i, j+1] + beta[s_i, c_i, i+1, j]) + \
                        beta[s_i, c_i, i, j] * \
                        (np.exp(tau*soft_sigma)-np.exp(2*tau))
    return alpha[:, :, -1, -1]/beta[:, :, -1, -1]

# Generate noisy DNA sequence
num_x = 300
pi = 0.1
pd = 0.1
ps = 0.1
alphabet = np.array(['T', 'A', 'G', 'C'])
sequences = ['ATTGCATA', 'GCTACCCA',
             'GCGAATCG', 'AGATTAAC']
# Repeat sequences n times
tx_strands = sequences * num_x
# Get output from IDS channel
X = ids_channel(tx_strands, pi, pd, ps)

centroids = sequences
labels = [i for i in range(len(centroids))] * num_x

alphabet_dict = {alphabet[i]: i + 1 for i in range(len(alphabet))}

max_length = np.max([len(seq) for seq in X] + [len(seq) for seq in centroids])
encoded_x = np.zeros((len(X), max_length, len(alphabet)), dtype=np.uint8)
encoded_c = np.zeros((len(centroids), max_length, len(alphabet)), dtype=np.uint8)
# Encode sequence with one-hot encoding
for i, s in enumerate(X):
    encoded_x[i, :, :] = onehot(s, alphabet, max_length)
for i, s in enumerate(centroids):
    encoded_c[i, :, :] = onehot(s, alphabet, max_length)

encoded_x = np.vstack((encoded_x, encoded_c))
encoded_x = np.array(encoded_x)
dist = ed(encoded_x, encoded_c)
print(dist)

tsne = TSNE(n_iter=1000, perplexity=50)
points = tsne.fit_transform(dist)
labels = np.concatenate((labels, np.full(len(centroids), len(centroids), np.int32)), axis=0)
colors = np.array(list(islice(cycle(['#377eb8', '#ff7f00', '#4daf4a',
                                        '#f781bf', '#a65628', '#984ea3',
                                        '#999999', '#e41a1c', '#dede00']),
                                int(len(centroids) + 1))))
colors = colors[labels]
print(points.shape)
indexes = np.arange(0, len(X))
plt.figure()
plt.scatter(points[indexes, 0], points[indexes, 1], c=colors[indexes], s=5, alpha=0.8)
for i in indexes:
    plt.text(points[i, 0], points[i, 1], X[i], color=colors[i], horizontalalignment='center',
            verticalalignment='bottom', fontsize=7, alpha=0.8)
plt.scatter(points[len(X):, 0], points[len(X):, 1], c='black', s=20)

for i in range(len(centroids)):
    plt.text(points[len(X) + i, 0], points[len(X) + i, 1], centroids[i], color='black',
            horizontalalignment='center', verticalalignment='bottom', fontsize=11)
plt.show()



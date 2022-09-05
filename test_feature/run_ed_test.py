
from http.cookiejar import CookieJar
from traceback import print_tb
import matplotlib.pyplot as plt
import numpy as np
from channel import ids_channel
from sklearn.manifold import TSNE
from itertools import cycle, islice

# Test edit distance
def levenshteinDistanceDP(s, s_len, t,  t_len):
    edit_m = np.zeros((s_len + 1, t_len + 1))
    edit_m[:, 0] = np.arange(s_len + 1)
    edit_m[0, :] = np.arange(t_len + 1)
    for j in range(t_len):
        for i in range(s_len):
            if s[i] == t[j]:
                edit_m[i+1, j+1] = edit_m[i, j]
            else:
                edit_m[i+1, j+1] = min([edit_m[i, j+1]+1,
                                        edit_m[i+1, j]+1,
                                        edit_m[i, j]+1, ])
    # print(edit_m)
    return edit_m[-1, -1]

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

X = X + centroids

dist = np.zeros((len(X), len(centroids)))
for i, x in enumerate(X):
    for j, c in enumerate(centroids):
        dist[i, j] = levenshteinDistanceDP(x, len(x), c, len(c))

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
indexes = np.arange(0, len(tx_strands))
plt.figure()
plt.scatter(points[indexes, 0], points[indexes, 1], c=colors[indexes], s=5, alpha=0.8)
for i in indexes:
    plt.text(points[i, 0], points[i, 1], X[i], color=colors[i], horizontalalignment='center',
            verticalalignment='bottom', fontsize=7, alpha=0.8)
plt.scatter(points[len(X):, 0], points[len(X):, 1], c='black', s=20)

for i in range(len(centroids)):
    plt.text(points[len(tx_strands) + i, 0], points[len(tx_strands) + i, 1], centroids[i], color='black',
            horizontalalignment='center', verticalalignment='bottom', fontsize=11)
plt.show()



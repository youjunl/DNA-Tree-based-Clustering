# Adapted from https://github.com/DeMoriarty/fast_pytorch_kmeans/blob/master/fast_pytorch_kmeans/kmeans.py

import math
import torch
from time import time
import numpy as np


class KMeans():
    '''
    Kmeans clustering algorithm implemented with PyTorch
    Parameters:
      n_clusters: int, 
        Number of clusters
      max_iter: int, default: 100
        Maximum number of iterations
      tol: float, default: 0.0001
        Tolerance

      verbose: int, default: 0
        Verbosity
      mode: {'euclidean', 'cosine'}, default: 'euclidean'
        Type of distance measure
      minibatch: {None, int}, default: None
        Batch size of MinibatchKmeans algorithm
        if None perform full KMeans algorithm

    Attributes:
      centroids: torch.Tensor, shape: [n_clusters, n_features]
        cluster centroids
    '''

    def __init__(self, n_clusters, max_iter=100, tol=0.0001, verbose=0, minibatch=None):
        self.n_clusters = n_clusters
        self.max_iter = max_iter
        self.tol = tol
        self.verbose = verbose
        self.minibatch = minibatch
        self._loop = False
        self._show = False

        try:
            import PYNVML
            self._pynvml_exist = True
        except ModuleNotFoundError:
            self._pynvml_exist = False

        self.centroids = None

    def sed(self, x, center, tau=-2):
        """
          Compute soft edit distance of one-hot encoded sequences
          Parameters:
          x: torch.Tensor, shape: [n, maxLen, alphabetLen]
          center: torch.Tensor, shape: [n, maxLen, alphabetLen]
          tau: float
        """
        maxLen = center.size(dim=1)
        nSample = x.size(dim=0)
        nCenter = center.size(dim=0)
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
        return torch.Tensor(alpha[:, :, -1, -1]/beta[:, :, -1, -1])

    def nobias_sed(self, x, center, tau=-2):
        mul_sed = self.sed(x, center)
        x_sed = self.sed(x, center)
        y_sed = self.sed(x, center)
        return mul_sed - 0.5*(x_sed+y_sed)

    def remaining_memory(self):
        """
          Get remaining memory in gpu
        """
        torch.cuda.synchronize()
        torch.cuda.empty_cache()
        if self._pynvml_exist:
            pynvml.nvmlInit()
            gpu_handle = pynvml.nvmlDeviceGetHandleByIndex(0)
            info = pynvml.nvmlDeviceGetMemoryInfo(gpu_handle)
            remaining = info.free
        else:
            remaining = torch.cuda.memory_allocated()
        return remaining

    def max_sim(self, a, b):
        """
          Compute maximum similarity (or minimum distance) of each vector
          in a with all of the vectors in b
          Parameters:
          a: torch.Tensor, shape: [nSample, maxLen, alphabetLen]
          b: torch.Tensor, shape: [nCenter, maxLen, alphabetLen]
        """
        device = a.device.type
        batch_size = a.shape[0]
        if device == 'cpu':
            sim = self.nobias_sed(a, b)
            max_sim_v, max_sim_i = sim.max(dim=-1)
            return max_sim_v, max_sim_i
        else:
            if a.dtype == torch.float:
                expected = a.shape[0] * a.shape[1] * b.shape[0] * 4
            elif a.dtype == torch.half:
                expected = a.shape[0] * a.shape[1] * b.shape[0] * 2
            ratio = math.ceil(expected / self.remaining_memory())
            subbatch_size = math.ceil(batch_size / ratio)
            msv, msi = [], []
            for i in range(ratio):
                if i*subbatch_size >= batch_size:
                    continue
                sub_x = a[i*subbatch_size: (i+1)*subbatch_size]
                sub_sim = self.nobias_sed(sub_x, b)
                sub_max_sim_v, sub_max_sim_i = sub_sim.min(dim=-1)
                del sub_sim
                msv.append(sub_max_sim_v)
                msi.append(sub_max_sim_i)
            if ratio == 1:
                max_sim_v, max_sim_i = msv[0], msi[0]
            else:
                max_sim_v = torch.cat(msv, dim=0)
                max_sim_i = torch.cat(msi, dim=0)
            return max_sim_v, max_sim_i

    def fit_predict(self, X, centroids=None):
        """
          Combination of fit() and predict() methods.
          This is faster than calling fit() and predict() seperately.
          Parameters:
          X: torch.Tensor, shape: [n_samples, n_features]
          centroids: {torch.Tensor, None}, default: None
            if given, centroids will be initialized with given tensor
            if None, centroids will be randomly chosen from X
          Return:
          labels: torch.Tensor, shape: [n_samples]
        """
        batch_size, maxLen, alphabetLen = X.shape
        device = X.device.type
        start_time = time()
        if centroids is None:
            self.centroids = X[np.random.choice(
                batch_size, size=[self.n_clusters], replace=False)]
        else:
            self.centroids = centroids
        num_points_in_clusters = torch.ones(self.n_clusters, device=device)
        closest = None
        for iter_num in range(self.max_iter):
            iter_time = time()
            if self.minibatch is not None:
                x = X[np.random.choice(
                    batch_size, size=[self.minibatch], replace=False)]
            else:
                x = X
            closest = self.max_sim(a=x, b=self.centroids)[1]
            matched_clusters, counts = closest.unique(return_counts=True)
            c_grad = torch.zeros_like(self.centroids)
            if self._loop:
                for j, count in zip(matched_clusters, counts):
                    c_grad[j] = x[closest == j].sum(dim=0) / count
            else:
                if self.minibatch is None:
                    expanded_closest = closest[None].expand(
                        self.n_clusters, -1)
                    mask = (expanded_closest == torch.arange(
                        self.n_clusters, device=device)[:, None]).float()
                    c_grad = mask @ x / mask.sum(-1)[..., :, None]
                    c_grad[c_grad != c_grad] = 0  # remove NaNs
                else:
                    expanded_closest = closest[None].expand(
                        len(matched_clusters), -1)
                    # mask: one-hot encoded clustering result
                    mask = (expanded_closest ==
                            matched_clusters[:, None]).float()
                    # Modified: update the gradient using clustered samples
                    for i in range(len(matched_clusters)):
                        tmp = 0
                        for j in range(len(mask[i,:])):
                            if mask[i, j]:
                                tmp += x[j, :, :]
                        c_grad[matched_clusters[i]] = tmp

            error = (c_grad - self.centroids).pow(2).sum()
            print(c_grad)
            print(self.centroids)
            if self.minibatch is not None:
                lr = 1/num_points_in_clusters[:, None] * 0.9 + 0.1
                # lr = 1/num_points_in_clusters[:,None]**0.1
            else:
                lr = 1
            num_points_in_clusters[matched_clusters] += counts
            # Update centroid 
            for i in range(len(self.centroids)):
                self.centroids[i] = self.centroids[i] * (1-lr[i]) + c_grad[i] * lr[i]
            if self.verbose >= 2:
                print('iter:', iter_num, 'error:', error.item(),
                      'time spent:', round(time()-iter_time, 4))
            if error <= self.tol:
                break

        # END SCATTER
        if self.verbose >= 1:
            print(
                f'used {iter_num+1} iterations ({round(time()-start_time, 4)}s) to cluster {batch_size} items into {self.n_clusters} clusters')
        return closest

    def predict(self, X):
        """
          Predict the closest cluster each sample in X belongs to
          Parameters:
          X: torch.Tensor, shape: [n_samples, n_features]
          Return:
          labels: torch.Tensor, shape: [n_samples]
        """
        return self.max_sim(a=X, b=self.centroids)[1]

    def fit(self, X, centroids=None):
        """
          Perform kmeans clustering
          Parameters:
          X: torch.Tensor, shape: [n_samples, n_features]
        """
        self.fit_predict(X, centroids)

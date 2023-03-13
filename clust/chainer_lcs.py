import numpy as np
import cupy

alphabet = {'A':0, 'T':1, 'G':2, 'C':3}

def mutual_lcs(X, length):
    X = one_hot_encoding(X, alphabet, length)
    X = cupy.array(X)
    n = len(X)
    I = np.ravel(np.broadcast_to(np.arange(n), (n, n)).T)
    J = np.ravel(np.broadcast_to(np.arange(n), (n, n)))
    d = lcs(X[I], X[J])
    d = cupy.asnumpy(d).reshape(n, n)
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
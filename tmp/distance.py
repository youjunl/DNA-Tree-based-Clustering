import numpy as np


def levenshteinDistance(s, s_len, t,  t_len):
    # Paper:Clustering Billions of Reads for DNA Data Storage
    if s_len == 0 or t_len == 0:
        return max(s_len, t_len)
    if s[s_len - 1] == t[t_len - 1]:
        cost = 0
    else:
        cost = 1
    return min([levenshteinDistance(s, s_len - 1, t, t_len) + 1,
                levenshteinDistance(s, s_len, t, t_len - 1) + 1,
                levenshteinDistance(s, s_len - 1, t, t_len - 1) + cost])


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


def soft_levenshteinDistance(s, s_len, t,  t_len, tau):
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
            if s[i] == t[j]:
                soft_sigma = 0
            else:
                soft_sigma = 1
            alpha[i+1, j+1] = np.exp(tau)*(alpha[i, j+1] + alpha[i+1, j] + beta[i, j+1] + beta[i+1, j]) + np.exp(
                tau*soft_sigma)*(alpha[i, j] + beta[i, j]*soft_sigma) - np.exp(2*tau)*(alpha[i, j] + 2*beta[i, j])
            beta[i+1, j+1] = np.exp(tau)*(beta[i, j+1] + beta[i+1, j]) + \
                beta[i, j]*(np.exp(tau*soft_sigma)-np.exp(2*tau))
    return alpha[-1, -1]/beta[-1, -1]


def nobiased_soft_levenshteinDistance(s, s_len, t,  t_len, tau):
    mul_sed = soft_levenshteinDistance(s, s_len, t,  t_len, tau)
    x_sed = soft_levenshteinDistance(s, s_len, s, s_len, tau)
    y_sed = soft_levenshteinDistance(t,  t_len, t,  t_len, tau)
    return mul_sed - 0.5*(x_sed+y_sed)

def onehot(string, alphabet, max_length):
    n = len(string)
    m = len(alphabet)
    result = np.zeros((max_length, m))
    for i in range(n):
        for j in range(m):
            if string[i] == alphabet[j]:
                result[i, j] = 1
            else:
                result[i, j] = 0
    return result


def softEditDistance(X1: np.array, X2: np.array, tau):
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

if __name__ == '__main__':
    print('--Performing testing for Levenshtein distance function--')
    seq1 = '011023'
    seq2 = '010121'

    print(levenshteinDistance(seq1, len(seq1), seq2, len(seq2)))
    print(levenshteinDistanceDP(seq1, len(seq1), seq2, len(seq2)))
    print(nobiased_soft_levenshteinDistance(seq1, len(seq1), seq2, len(seq2), -2))
    print(soft_levenshteinDistance(seq1, len(seq1), seq2, len(seq2), -2))
    print(softEditDistance(onehot(seq1, '0123', 8), onehot(seq2, '0123', 8), -2))

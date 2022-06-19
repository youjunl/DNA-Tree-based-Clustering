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


if __name__ == '__main__':
    print('--Performing testing for Levenshtein distance function--')
    seq1 = 'GTTGCA'
    seq2 = 'GATCCA'

    print(levenshteinDistance(seq1, len(seq1), seq2, len(seq2)))
    print(levenshteinDistanceDP(seq1, len(seq1), seq2, len(seq2)))

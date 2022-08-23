import numpy as np


def accuracy(est_c, c, gamma):
    # Paper:Clustering Billions of Reads for DNA Data Storage
    est_c_len = len(est_c)
    c_len = len(c)
    tmp = 0
    for i in range(est_c_len):
        print(set(est_c[i]).issubset(c[i]))
        print(len(set(est_c[i]).intersection(c[i])))
        if set(est_c[i]).issubset(c[i]) and len(set(est_c[i]).intersection(c[i])) >= gamma*len(c[i]):
            tmp = tmp + 1
    acc = tmp/c_len
    return acc


if __name__ == '__main__':
    print('--Performing testing for accuracy function--')
    est_c = [[0, 1, 2, 3, 4],
            [5, 6],
            [7, 8, 9, 10]]
    c = [[0, 1, 2, 3, 4],
         [5, 6],
         [7, 8, 9],
         [10]]
    print(accuracy(est_c, c, gamma = 0.9))
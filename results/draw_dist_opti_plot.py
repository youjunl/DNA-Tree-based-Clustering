import matplotlib
import matplotlib.pyplot as plt
import scienceplots
import numpy as np
plt.style.use('ieee')

files = ['acc_dist_1.csv', 'acc_dist_2.csv', 'acc_dist_3.csv']
def compute(infile):
    f = open(infile, 'r')
    lines = f.readlines()
    f.close()
    gamma = list(map(float, lines[0].strip().split(',')[1:-1]))
    th_0_acc = list(map(float, lines[1].strip().split(',')[1:-1]))
    th_1_acc = list(map(float, lines[2].strip().split(',')[1:-1]))
    th_2_acc = list(map(float, lines[3].strip().split(',')[1:-1]))
    th_3_acc = list(map(float, lines[4].strip().split(',')[1:-1]))
    return gamma, th_0_acc, th_1_acc, th_2_acc, th_3_acc

if __name__ == '__main__':
    fig, ax = plt.subplots(figsize=[4, 4])
    gamma_pos = -5
    dist = [1, 2, 3]
    th_0_acc_list = [0, 0, 0]
    th_1_acc_list = [0, 0, 0]
    th_2_acc_list = [0, 0, 0]
    th_3_acc_list = [0, 0, 0]
    for i in range(1, 4):
        gamma, th_0_acc, th_1_acc, th_2_acc, th_3_acc = compute(files[i-1])
        th_0_acc_list[i-1] = th_0_acc[gamma_pos]
        th_1_acc_list[i-1] = th_1_acc[gamma_pos]
        th_2_acc_list[i-1] = th_2_acc[gamma_pos]
        th_3_acc_list[i-1] = th_3_acc[gamma_pos]
    ax.plot(dist, th_0_acc_list, '-*', label='DLS')
    ax.plot(dist, th_1_acc_list, '--*', label='DLSM $\\theta_s$=1')
    ax.plot(dist, th_2_acc_list, '--*', label='DLSM $\\theta_s$=2')
    ax.plot(dist, th_3_acc_list, '--*', label='DLSM $\\theta_s$=3')
    # fig.text(0.25, 0.1, '(a)', horizontalalignment='center')
    # fig.text(0.5, 0.1, '(b)', horizontalalignment='center')
    # fig.text(0.75, 0.1, '(c)', horizontalalignment='center')
    ax.set_xlabel('Minimum Levenshtein Distance of Index')
    ax.set_ylabel('Accuracy (γ=%.2f)'%gamma[gamma_pos])
    ax.set_xticks(dist)
    ax.legend()
    ax.grid(which='both')
    # plt.gcf().set_dpi(600)
    plt.savefig('dist_optimization.pdf')

    fig, ax = plt.subplots(figsize=[4, 4])
    gamma_pos = -5
    dist = [1, 2, 3]
    th_0_acc_list = [0, 0, 0]
    th_1_acc_list = [0, 0, 0]
    th_2_acc_list = [0, 0, 0]
    th_3_acc_list = [0, 0, 0]
    for i in range(1, 4):
        gamma, th_0_acc, th_1_acc, th_2_acc, th_3_acc = compute(files[i-1])
        th_0_acc_list[i-1] = th_0_acc[gamma_pos]
        th_1_acc_list[i-1] = th_1_acc[gamma_pos]
        th_2_acc_list[i-1] = th_2_acc[gamma_pos]
        th_3_acc_list[i-1] = th_3_acc[gamma_pos]

    arr = np.array([th_0_acc_list, th_1_acc_list, th_2_acc_list, th_3_acc_list])
    th = [0, 1, 2, 3]
    print(arr)
    ax.plot(th, arr[:, 0], '-*', label='Min. LD = 1')
    ax.plot(th, arr[:, 1], '--*', label='Min. LD = 2')
    ax.plot(th, arr[:, 2], '--*', label='Min. LD = 3')
    # fig.text(0.25, 0.1, '(a)', horizontalalignment='center')
    # fig.text(0.5, 0.1, '(b)', horizontalalignment='center')
    # fig.text(0.75, 0.1, '(c)', horizontalalignment='center')
    ax.set_xlabel('Levenshtein Distance Threshold')
    ax.set_ylabel('Accuracy (γ=%.2f)'%gamma[gamma_pos])
    ax.set_xticks(th)
    ax.legend()
    ax.grid(which='both')
    # plt.gcf().set_dpi(600)
    plt.savefig('dist_optimization_1.pdf')
    plt.show()
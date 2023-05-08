import matplotlib.pyplot as plt
import scienceplots
from scipy.special import comb
import numpy as np
plt.style.use('ieee')

if __name__ == '__main__':
    file = 'number_of_clusters_10k.csv'
    check_points = []
    dls_data = []
    dlsm_data = []
    lines = open(file, 'r').readlines()
    check_points = list(map(int, lines[0].strip().split(',')[1:-1]))
    dls_data = list(map(int, lines[1].strip().split(',')[1:-1]))
    dlsm_data_1 = list(map(int, lines[2].strip().split(',')[1:-1]))
    dlsm_data_2 = list(map(int, lines[3].strip().split(',')[1:-1]))
    dlsm_data_3 = list(map(int, lines[4].strip().split(',')[1:-1]))

    fig, ax = plt.subplots(figsize=[5,4])
    ax.plot(check_points, dls_data, label = 'DLS')
    ax.plot(check_points, dlsm_data_1, label = 'DLSM $\\theta_s$=1')
    ax.plot(check_points, dlsm_data_2, label = 'DLSM $\\theta_s$=2')
    ax.plot(check_points, dlsm_data_3, label = 'DLSM $\\theta_s$=3')
    ax.get_xaxis().get_major_formatter().set_scientific(True)
    ax.get_yaxis().get_major_formatter().set_scientific(True)
    ax.grid()
    ax.set_xlabel('Number of reads')
    ax.set_ylabel('Number of clusters')
    ax.legend()
    ax.legend()
    plt.gcf().set_dpi(300)
    plt.savefig('num_cl_plot.pdf')
    plt.show()
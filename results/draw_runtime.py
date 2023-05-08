import matplotlib.pyplot as plt
import scienceplots
plt.style.use('ieee')
# DLS configurations:
# config_dict={
#     "end_tree_len" : 20,
#     "tree_threshold" : 4,
#     "sub_tree_threshold" : 2,
#     "h_index_nums" : 0,
#     "train_start": 0,
#     "train_end": 0
# }

# DLSM configurations:
# config_dict={
#     "end_tree_len" : 20,
#     "tree_threshold" : 4,
#     "sub_tree_threshold" : 2,
#     "h_index_nums" : 0,
#     "train_start": 0,
#     "train_end": 10000000
# }
if __name__ == '__main__':
    lines = open('runtime.csv', 'r').readlines()
    size = list(map(int, lines[0].strip().split(',')[1:]))
    l1 = list(map(float, lines[1].strip().split(',')[1:]))
    l2 = list(map(float, lines[2].strip().split(',')[1:]))
    l3 = list(map(float, lines[3].strip().split(',')[1:]))
    l4 = list(map(float, lines[4].strip().split(',')[1:]))
    l5 = list(map(float, lines[5].strip().split(',')[1:]))
    l6 = list(map(float, lines[6].strip().split(',')[1:]))
    fig, ax = plt.subplots(figsize=(5,5))
    ax.plot(size, l1, '-ro', label='DLS L=16')
    ax.plot(size, l2, '-go', label='DLS L=18')
    ax.plot(size, l3, '-bo', label='DLS L=20')
    ax.plot(size, l4, '--r*', label='DLSM L=16')
    ax.plot(size, l5, '--g*', label='DLSM L=18')
    ax.plot(size, l6, '--b*', label='DLSM L=20')
    ax.grid(which='both')
    ax.set_xlabel('Number of reads (log scale)')
    ax.set_ylabel('Runtime (seconds, log scale)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.legend()
    plt.gcf().set_dpi(300)
    plt.savefig('runtime.pdf')
    plt.show()
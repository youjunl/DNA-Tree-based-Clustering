import matplotlib.pyplot as plt
import scienceplots
plt.style.use('ieee')
# DLS configurations:
# config_dict={
#     "end_tree_len" : 21,
#     "tree_threshold" : 4,
#     "sub_tree_threshold" : 2,
#     "h_index_nums" : 0,
#     "train_start": 0,
#     "train_end": 0
# }

# DLSM configurations:
# config_dict={
#     "end_tree_len" : 21,
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
    fig, ax = plt.subplots(figsize=(5,5))
    ax.plot(size, l1, '-o', label='DLS')
    ax.plot(size, l2, '-*', label='DLSM')
    ax.grid(which='both')
    ax.set_xlabel('Number of reads (log scale)')
    ax.set_ylabel('Runtime (seconds, log scale)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.legend()
    plt.gcf().set_dpi(300)
    plt.savefig('runtime.pdf')
    plt.show()
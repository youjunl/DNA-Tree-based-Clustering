import matplotlib.pyplot as plt
import scienceplots
plt.style.use('ieee')

if __name__ == '__main__':
    dls_files = ['mem_usage_dls_13.csv','mem_usage_dls_17.csv','mem_usage_dls_21.csv']
    dlsm_files = ['mem_usage_dlsm_13.csv','mem_usage_dlsm_17.csv','mem_usage_dlsm_21.csv']
    check_points = []
    dls_data = []
    dlsm_data = []
    for i in range(3):
        lines = open(dls_files[i], 'r').readlines()
        check_points = list(map(float, lines[0].strip().split(',')[:-1]))
        data = list(map(float, lines[1].strip().split(',')[:-1]))
        dls_data.append(data)

        lines = open(dlsm_files[i], 'r').readlines()
        check_points = list(map(float, lines[0].strip().split(',')[:-1]))
        data = list(map(float, lines[1].strip().split(',')[:-1]))
        dlsm_data.append(data)

    fig, ax = plt.subplots(1, 2, figsize=[8,4])
    ax[0].plot(check_points, dls_data[0], label = 'DLS L=13')
    ax[0].plot(check_points, dls_data[1], label = 'DLS L=17')
    ax[0].plot(check_points, dls_data[2], label = 'DLS L=21')

    ax[1].plot(check_points, dlsm_data[0], label = 'DLSM L=13')
    ax[1].plot(check_points, dlsm_data[1], label = 'DLSM L=17')
    ax[1].plot(check_points, dlsm_data[2], label = 'DLSM L=21')

    ax[0].grid()
    ax[1].grid()
    ax[0].set_xlabel('Number of reads')
    ax[1].set_xlabel('Number of reads')
    ax[0].set_ylabel('Memory usage (gigabyte)')
    ax[0].legend()
    ax[1].legend()
    plt.gcf().set_dpi(300)
    plt.savefig('mem_plot.pdf')
    plt.show()
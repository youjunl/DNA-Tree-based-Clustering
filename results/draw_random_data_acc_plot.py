import matplotlib
import matplotlib.pyplot as plt
import scienceplots
plt.style.use('ieee')

if __name__ == '__main__':
    fig, ax = plt.subplots(figsize=[5,4])
    lines = open('acc_dist_1.csv', 'r').readlines()
    gamma = list(map(float, lines[0].strip().split(',')[1:-1]))
    dls = list(map(float, lines[1].strip().split(',')[1:-1]))
    dlsm_1 = list(map(float, lines[2].strip().split(',')[1:-1]))
    dlsm_2 = list(map(float, lines[3].strip().split(',')[1:-1]))
    dlsm_3 = list(map(float, lines[4].strip().split(',')[1:-1]))

    ax.plot(gamma, dls, label='DLS')
    ax.plot(gamma, dlsm_1, label='DLSM $\\theta_s$=1')
    ax.plot(gamma, dlsm_2, label='DLSM $\\theta_s$=2')
    ax.plot(gamma, dlsm_3, label='DLSM $\\theta_s$=3')
    ax.grid(which='both')
    ax.set_xlim(0, 1)
    ax.legend()
    ax.set_xlabel('γ')
    ax.set_ylabel('Accuracy')
    plt.gcf().set_dpi(300)
    plt.savefig('random_data_accuracy.pdf')
    # plt.show()
import matplotlib
import matplotlib.pyplot as plt
import scienceplots
plt.style.use('ieee')
matplotlib.use('pgf')

files = ['ERR1816980_acc.csv', 'P10_acc.csv', 'id20_acc.csv']
def compute(infile):
    f = open(infile, 'r')
    line = f.readline()
    gamma = list(map(float, line.strip().split(',')[1:-1]))
    line = f.readline()
    dls_acc = list(map(float, line.strip().split(',')[1:-1]))
    line = f.readline()
    dlsm_acc = list(map(float, line.strip().split(',')[1:-1]))
    line = f.readline()
    clover_acc = list(map(float, line.strip().split(',')[1:-1]))
    line = f.readline()
    starcode_acc = list(map(float, line.strip().split(',')[1:-1]))
    line = f.readline()
    mmseqs_acc = list(map(float, line.strip().split(',')[1:-1]))
    f.close()
    return gamma, dls_acc, dlsm_acc, clover_acc, starcode_acc, mmseqs_acc

if __name__ == '__main__':
    fig, ax = plt.subplots(1,3, sharey=True, figsize=[8, 3])
    labels = ['(a)', '(b)', '(c)']
    for i in range(3):
        gamma, dls_acc, dlsm_acc, clover_acc, starcode_acc, mmseqs_acc = compute(files[i])
        ax[i].plot(gamma, dls_acc, label='DLS')
        ax[i].plot(gamma, dlsm_acc, label='DLSM')
        ax[i].plot(gamma, clover_acc, label='Clover')
        ax[i].plot(gamma, starcode_acc, label='Starcode')
        ax[i].plot(gamma, mmseqs_acc, '--y', label='MMSeqs2')
        ax[i].grid(which='both')
        ax[i].set_xlim(0, 1)
        ax[i].set_aspect(1)
        ax[i].legend()
        ax[i].set_xlabel(labels[i])

    # fig.text(0.25, 0.1, '(a)', horizontalalignment='center')
    # fig.text(0.5, 0.1, '(b)', horizontalalignment='center')
    # fig.text(0.75, 0.1, '(c)', horizontalalignment='center')
    fig.supxlabel('γ')
    fig.supylabel('Accuracy')
    fig.tight_layout()
    # plt.gcf().set_dpi(600)
    plt.savefig('accuracy.pdf')
    # plt.show()
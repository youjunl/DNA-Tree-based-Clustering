import matplotlib.pyplot as plt
import scienceplots
import numpy as np
plt.style.use('ieee')

def plot_hist(ax, x, y, label):
    ax.bar(x, y, width=0.002, alpha=0.7, label=label, log=True)

if __name__ == '__main__':
    lines = open('number_of_errors.csv', 'r').readlines()
    name = ['ERR1816980', 'P10-5-BDDP210000009', 'ID20']
    x =  list(map(float, lines[0].strip().split(',')[:-1]))
    l1 = list(map(float, lines[1].strip().split(',')[:-1]))
    l2 = list(map(float, lines[2].strip().split(',')[:-1]))
    l3 = list(map(float, lines[3].strip().split(',')[:-1]))
    fig, ax = plt.subplots(figsize=[4,4])

    ax.grid()
    plot_hist(ax, x, l1, 'ERR1816980')
    plot_hist(ax, x, l2, 'P10-5-BDDP210000009')
    plot_hist(ax, x, l3, 'ID20')
    ax.set_xlabel('Read error rate')
    ax.set_ylabel('Number of references')
    ax.set_xlim(0, 1)
    plt.legend()
    plt.gcf().set_dpi(300)
    plt.savefig('number_of_errors.pdf')
    plt.show()
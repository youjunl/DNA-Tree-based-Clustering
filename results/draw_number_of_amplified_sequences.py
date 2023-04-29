import matplotlib.pyplot as plt
import scienceplots
plt.style.use('ieee')

def plot_hist(ax, x, label):
    bins = list(range(1, len(x) + 1))
    ax.bar(bins, x, alpha=0.7, label=label, log=True)

if __name__ == '__main__':
    lines = open('number_of_amplified_sequences.csv', 'r').readlines()
    name = ['ERR1816980', 'P10-5-BDDP210000009', 'ID20']
    l1 = list(map(int, lines[1].strip().split(',')[:-1]))
    l2 = list(map(int, lines[2].strip().split(',')[:-1]))
    l3 = list(map(int, lines[3].strip().split(',')[:-1]))
    fig, ax = plt.subplots(figsize=[4,4])
    ax.grid()
    plot_hist(ax, l1, 'ERR1816980')
    plot_hist(ax, l2, 'P10-5-BDDP210000009')
    plot_hist(ax, l3, 'ID20')
    ax.set_xlabel('Number of amplified sequences')
    ax.set_ylabel('Number of references')
    plt.legend()
    plt.gcf().set_dpi(300)
    plt.savefig('number_of_amplified_sequences.pdf')
    plt.show()
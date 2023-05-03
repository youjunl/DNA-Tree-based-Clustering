import matplotlib.pyplot as plt
import scienceplots
plt.style.use('ieee')

def plot_hist_amp(ax, x, label):
    bins = list(range(1, len(x) + 1))
    ax.bar(bins, x, alpha=0.7, label=label, log=True)


lines = open('number_of_amplified_sequences.csv', 'r').readlines()
name = ['ERR1816980', 'P10-5-BDDP210000009', 'ID20']
l1 = list(map(int, lines[1].strip().split(',')[:-1]))
l2 = list(map(int, lines[2].strip().split(',')[:-1]))
l3 = list(map(int, lines[3].strip().split(',')[:-1]))
fig, ax = plt.subplots(1,2,figsize=[8,4])
ax[0].grid()
plot_hist_amp(ax[0], l1, 'ERR1816980')
plot_hist_amp(ax[0], l2, 'P10-5-BDDP210000009')
plot_hist_amp(ax[0], l3, 'ID20')
ax[0].set_xlabel('Number of amplified sequences')
ax[0].set_ylabel('Number of references')
ax[0].legend()
# plt.gcf().set_dpi(300)
# plt.savefig('number_of_amplified_sequences.pdf')
# plt.show()


def plot_hist_rate(ax, x, y, label):
    ax.bar(x, y, width=0.002, alpha=0.7, label=label, log=True)


lines = open('number_of_errors.csv', 'r').readlines()
name = ['ERR1816980', 'P10-5-BDDP210000009', 'ID20']
x =  list(map(float, lines[0].strip().split(',')[:-1]))
l1 = list(map(float, lines[1].strip().split(',')[:-1]))
l2 = list(map(float, lines[2].strip().split(',')[:-1]))
l3 = list(map(float, lines[3].strip().split(',')[:-1]))

ax[1].grid()
plot_hist_rate(ax[1], x, l1, 'ERR1816980')
plot_hist_rate(ax[1], x, l2, 'P10-5-BDDP210000009')
plot_hist_rate(ax[1], x, l3, 'ID20')
ax[1].set_xlabel('Read error rate')
ax[1].set_ylabel('Number of references')
ax[1].set_xlim(0, 1)
ax[1].legend()
# plt.gcf().set_dpi(600)
plt.savefig('data_characterization.pdf')
# plt.show()
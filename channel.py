import numpy as np


def ids_channel(tx_seqs, pi, pd, ps):
    numSeqs = len(tx_seqs)
    rx_seqs = [''] * numSeqs
    alphabet = ['A', 'T', 'G', 'C']
    # Markov process
    for ind, tx in enumerate(tx_seqs):
        N = len(tx)
        rx = ''
        i = 0  # state number
        while i < N:
            # Insertion
            if randsel(pi):
                rx = rx + alphabet[np.random.randint(0, 4)]
            else:
                # Deletion
                if ~randsel(pd):
                    # Substitution
                    if randsel(ps):
                        tmp = alphabet.copy()
                        tmp.remove(tx[i])
                        rx = rx + tmp[np.random.randint(0, 3)]
                    else:
                        rx = rx + tx[i]
                i = i + 1
        rx_seqs[ind] = rx
    return rx_seqs


def randsel(p):
    sel = np.random.rand() < p
    return sel


if __name__ == '__main__':
    # Test
    tx_seqs = ['AATGC', 'TTTTTTTTTTTTTTTT']
    rx_seqs = ids_channel(tx_seqs, 0.1, 0.1, 0.1)
    for tx, rx in zip(tx_seqs, rx_seqs):
        print('TX: ' + tx)
        print('RX: ' + rx)

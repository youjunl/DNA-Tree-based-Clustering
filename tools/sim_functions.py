import random
import numpy as np

dna_to_num = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
num_to_dna = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
num_to_dna_array = ['A', 'T', 'G', 'C']

# Implementation of IDS channel
def channel(txSeqs, pi, pd, ps):
    rxSeqs = []
    simErr = 0
    for ind, tx in enumerate(txSeqs):
        errorFlag = ''
        i, rx = 0, ''
        while i < len(tx):
            if randsel(pi):  # Insertion
                rx += num_to_dna_array[random.randint(0, 3)]
                errorFlag += 'e'
                continue
            if not randsel(pd):  # Deletion
                if randsel(ps):  # Substitution
                    tmp = num_to_dna_array.copy()
                    tmp.remove(tx[i])
                    rx += tmp[random.randint(0, 2)]
                    errorFlag += 'e'
                else:
                    rx += tx[i]
                    errorFlag += '_'
            else:
                errorFlag += 'e'
            i += 1
        rxSeqs.append(rx)
        if 'ee' in errorFlag[:n]:
            simErr += 1

    return rxSeqs, simErr

# Random selection
def randsel(p):
    sel = np.random.rand() < p
    return sel

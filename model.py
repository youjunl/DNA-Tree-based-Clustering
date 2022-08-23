from turtle import distance
import numpy as np
from accuracy import accuracy
from distance import levenshteinDistanceDP
from channel import ids_channel

class DNA_StrandSimulation():
    def __init__(
        self,
        seqs: list,
        num_x: int = 1000,
        num_class: int = 10,
        len_x: int = 10,
        pi: float = 0.1,
        pd: float = 0.1,
        ps: float = 0.1):
        # Generate random number copies of sequences
        # self.tx_strands = self.bits[np.random.randint(
        #     0, num_class, size=(1, num_x))][0]
        self.tx_strands = []
        self.rx_strands = []

    def get_tx(self):
        return self.tx_strands

    def get_rx(self):
        return self.rx_strands

    def get_distance(self):
        # Caculate distance from each sequence to each centers
        self.distance = np.zeros((num_x, num_class))
        for i in range(num_x):
            print('%d/%d' % (i+1, num_x), end='\r')
            tx = self.tx_strands[i]
            rx = ids_channel(tx, pi, pd, ps)
            self.rx_strands.append(rx)
            for j in range(num_class):
                distance = levenshteinDistanceDP(
                    self.bits[j], len_x, rx, len(rx))
                self.distance[i, j] = distance

        return self.distance


if __name__ == '__main__':
    c = DNA_StrandSimulation(
        num_x=1000,
        num_class=10,
        len_x=10,
        pi=0.1,
        pd=0.1,
        ps=0.1
    )

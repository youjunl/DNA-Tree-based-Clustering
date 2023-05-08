import matplotlib.pyplot as plt
import random
from scipy.stats import binom
import numpy as np
from scipy.special import comb
import scienceplots
plt.style.use('ieee')
p = 1-(1-0.01-0.01)*(1-0.01)
def S(L, D):
    if L == 1:
        return 0
    elif L == 2:
        return p**2 
    elif L <= D:
        tmp = 0
        for i in range(2, D+1):
            tmp += comb(D, i)*p**i*(1-p)**(D-i)
        return tmp
    else:
        return S(L-1, D)+p**2*(1-S(L-2, D))

def P(L, theta):
    # Calculate the cumulative probability P(X ≤ θ)
    cumulative_prob = binom.cdf(theta, L, p)
    # Calculate the probability P(X > θ)
    return 1 - cumulative_prob

# if __name__ == '__main__':
#     xL = [i for i in range(12, 22)]
#     yL_D2 = [S(L, 2) for L in xL]
#     yL_D3 = [S(L, 3) for L in xL]
#     yL_D4 = [S(L, 4) for L in xL]
    
#     fig, ax = plt.subplots(figsize=[4,4])
#     ax.plot(xL, yL_D2, label='D=2')
#     ax.plot(xL, yL_D3, label='D=3')
#     ax.plot(xL, yL_D4, label='D=4')
#     ax.grid()
#     plt.legend()
#     plt.gcf().set_dpi(300)
#     plt.savefig('lower_bound.pdf')
#     plt.show()
def randsel(p):
    sel = np.random.rand() < p
    return sel

simulations = 100000
fail_count = 0
success_count = 0
depth_limit = 3
for i in range(simulations):
    ind = 0
    count = 0
    for j in range(21):
        if ind == j-1:
            success_count += 1
        if count > 4:
            fail_count += 1
            break
        ind = j
        depth = 0
        if randsel(p):
            traverse_num = min(depth_limit, 21 - j - 1)
            while depth < depth_limit:
                if randsel(p):
                    break
                depth+=1
            if depth==depth_limit:
                count+=1
            else:
                fail_count+=1
                break
print(fail_count)
fail_rate = fail_count/simulations
print(fail_rate)
print(S(21, 3)+P(21, 4))
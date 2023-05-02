import matplotlib.pyplot as plt
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
    
if __name__ == '__main__':
    xL = [i for i in range(13, 22)]
    yL_D2 = [S(L, 2) for L in xL]
    yL_D3 = [S(L, 3) for L in xL]
    yL_D4 = [S(L, 4) for L in xL]
    
    fig, ax = plt.subplots(figsize=[4,4])
    ax.plot(xL, yL_D2, label='D=2')
    ax.plot(xL, yL_D3, label='D=3')
    ax.plot(xL, yL_D4, label='D=4')
    ax.grid()
    plt.legend()
    plt.gcf().set_dpi(300)
    plt.savefig('lower_bound.pdf')
    plt.show()
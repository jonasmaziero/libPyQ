import numpy as np
import math
import matplotlib.pyplot as plt
import coherence
import rpvg
from rpvg import rpv_zhsl
from rug import ru_gram_schmidt
from distances import normHS

def test():
    np.random.seed()
    ns = 10**3  # number of samples for the average
    nqb = 5  # maximum number of qubits regarded
    Cavg1 = np.zeros(nqb)
    Cavg2 = np.zeros(nqb)
    d = np.zeros(nqb, dtype=int)
    for j in range(0, nqb):
        d[j] = 2**(j+1)
        rs1 = np.zeros((d[j], d[j]), dtype=complex)
        rs2 = np.zeros((d[j], d[j]), dtype=complex)
        Cavg1[j] = 0.0; Cavg2[j] = 0.0
        for k in range(0, ns):
            rs1 = rdm(d[j],'gin')
            rs2 = rdm(d[j],'pos_')
            Cavg1[j] += coherence.coh(rs1,'l1')
            Cavg2[j] += coherence.coh(rs2,'l1')
        Cavg1[j] /= ns;  Cavg2[j] /= ns
        print(Cavg1[j], Cavg2[j])
    plt.plot(d, Cavg1, label='gin')
    plt.plot(d, Cavg2, label='pos_')
    plt.xlabel('d')
    plt.ylabel('C')
    plt.legend()
    plt.show()

def rdm(d,method):
    if method == 'std':
        return rdm_std(d)
    elif method == 'gin':
        return rdm_ginibre(d)
    elif method == 'pos':
        return rdm_pos(d)
    elif method == 'pos_':
        return rdm_pos_(d)

def rdm_pos_(d): # did not work
    rpv = np.zeros(d)
    rpv = rpvg.rpv_zhsl(d)
    rrho = np.zeros((d,d), dtype = complex)
    for j in range(0,d):
        rrho[j][j] = rpv[j]
    for j in range(0,d-1):
        for k in range(j+1,d):
            x, y = rand_circle(math.sqrt(abs(rrho[j][j]*rrho[k][k])))
            rrho[j][k] = x + 1j*y
            rrho[k][j] = x - 1j*y
    return rrho

def rdm_pos(d):
    rpv = np.zeros(d)
    rpv = rpvg.rpv_zhsl(d)
    rrho = np.zeros((d,d), dtype = complex)
    for j in range(0,d):
        rrho[j][j] = rpv[j]
    r = 0.0
    for j in range(0,d-1):
        for k in range(j+1,d):
            r += rpv[j]*rpv[k]
    r = math.sqrt(r)
    for j in range(0,d-1):
        for k in range(j+1,d):
            x, y = rand_circle(r)
            rrho[j][k] = x + 1j*y
            rrho[k][j] = x - 1j*y
            r -= (math.pow(x,2)+math.pow(y,2))
    return rrho

'''
def test_rcircle():
    n = 10000
    r = 1
    #th = np.zeros(n)
    #rh = np.zeros(n)
    #th = 2*math.pi*np.random.random(n)
    #rh = r*np.sqrt(np.random.random(n))
    #x = rh*np.cos(th)
    #y = rh*np.sin(th)
    x = np.zeros(n)
    y = np.zeros(n)
    for j in range(0,n):
        x[j], y[j] = rand_circle(r)
    plt.scatter(x, y, label='scatter', color='blue', s=5, marker='*')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()
    return
'''

def rand_circle(r):
    th = 2*math.pi*np.random.random()
    rh = r*math.sqrt(np.random.random())
    x = rh*math.cos(th);  y = rh*math.sin(th)
    return x, y

def rdm_std(d):
    rdm = np.zeros((d, d), dtype=complex)
    rpv = rpv_zhsl(d)
    ru = ru_gram_schmidt(d)
    for j in range(0, d):
        for k in range(j, d):
            for l in range(0, d):
                rdm[j][k] = rdm[j][k] + rpv[l]*ru[j][l].real*ru[k][l].real
                rdm[j][k] = rdm[j][k] + rpv[l]*ru[j][l].imag*ru[k][l].imag
                rdm[j][k] = rdm[j][k] + rpv[l]*(1j)*ru[j][l].imag*ru[k][l].real
                rdm[j][k] = rdm[j][k] - rpv[l]*(1j)*ru[j][l].real*ru[k][l].imag
                if j != k:
                    rdm[k][j] = rdm[j][k].real - (1j)*rdm[j][k].imag
    return rdm

def rdm_ginibre(d):
    rdm = np.zeros((d, d), dtype=complex)
    G = ginibre(d)
    N2 = (normHS(d, G))**2.0
    for j in range(0, d):
        for k in range(j, d):
            for l in range(0, d):
                rdm[j][k] = rdm[j][k] + (G[j][l].real)*(G[k][l].real) \
                    + (G[j][l].imag)*(G[k][l].imag) \
                    - (1j)*((G[j][l].real)*(G[k][l].imag)
                            - (G[j][l].imag)*(G[k][l].real))
            rdm[j][k] = rdm[j][k]/N2
            if j != k:
                rdm[k][j] = rdm[j][k].real - (1j)*rdm[j][k].imag
    return rdm

def ginibre(d):
    G = np.zeros((d, d), dtype=complex)
    mu, sigma = 0.0, 1.0
    for j in range(0, d):
        grn = np.random.normal(mu, sigma, 2*d)
        for k in range(0, d):
            G[j][k] = grn[k] + (1j)*grn[k+d]
    return G

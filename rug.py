import numpy as np


'''
def test():
    from numpy import sort
    from math import pi
    from cmath import phase
    import matplotlib.pyplot as plt
    d = 20
    ns = 10**4
    ni = 100
    egv = np.zeros(d, dtype=complex)
    egp = np.zeros(d)
    sp = np.zeros(d-1)
    x = np.zeros(ni)
    y1 = np.zeros(ni)
    y2 = np.zeros(ni)
    ct = np.zeros(ni, dtype=int)
    cts = np.zeros(ni, dtype=int)
    delta = (2.0*pi)/ni  # step for the eigenphases
    for j in range(0, ns):
        ru = ru_gram_schmidt(d)
        egv = np.linalg.eigvals(ru)
        for m in range(0, d):  # Puts the eigenphases in [0,2*pi] and sort them
            egp[m] = phase(egv[m]) + pi
        egp = sort(egp)
        for m in range(0, d-1):
            sp[m] = egp[m+1]-egp[m]
        spavg = sum(sp)/(d-1)
        sp = sp/spavg  # Normalize the spacings by the average spacing
        for k in range(0, d):
            for l in range(0, ni):
                if egp[k] >= (l-1.0)*delta and egp[k] < l*delta:
                    ct[l] = ct[l] + 1
                if k < d-1:
                    if sp[k] >= (l-1.0)*delta and sp[k] < l*delta:
                        cts[l] = cts[l] + 1
    for l in range(0, ni):
        x[l] = l*delta
        y1[l] = ct[l]/(ns*d)
        y2[l] = cts[l]/(ns*(d-1))
    plt.plot(x, y1, label='ev')
    plt.plot(x, y2, label='sp')
    plt.xlabel('')
    plt.ylabel('y')
    axes = plt.gca()
    axes.set_xlim([0.06, 2.5])
    axes.set_ylim([0, 0.065])
    plt.legend()
    plt.show()
'''


def ru_gram_schmidt(d):
    from rdmg import ginibre
    G = ginibre(d)
    return gram_schmidt_modified(d, G)


def gram_schmidt_modified(d, G):
    from distances import norm, inner
    B = np.zeros((d, d), dtype=complex)
    for j in range(0, d):
        B[:][j] = G[:][j]/norm(d, G[:][j])
        if j < (d-1):
            for k in range(j+1, d):
                G[:][k] = G[:][k] - inner(d, B[:][j], G[:][k])*B[:][j]
    return B

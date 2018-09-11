import numpy as np


def purity(d, rho):
    purity = 0.0
    j = -1
    while (j < d-1):
        j = j + 1
        k = -1
        while (k < d-1):
            k = k + 1
            purity += (rho[j][k].real)**2 + (rho[j][k].imag)**2
    return purity


def shannon(d, pv):
    from math import log
    SE = 0.0
    j = -1
    while (j < d-1):
        j = j + 1
        if pv[j] > 1.e-15 and pv[j] < (1.0-1.e-15):
            SE -= pv[j]*log(pv[j], 2)
    return SE


def neumann(d, rho):
    import scipy.linalg.lapack as lapak
    b = lapak.zheevd(rho)
    VnE = shannon(d, b[0])
    return VnE

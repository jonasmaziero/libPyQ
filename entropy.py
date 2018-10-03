from math import log
import scipy.linalg.lapack as lapak
import pTrace


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
    SE = 0.0
    j = -1
    while (j < d-1):
        j = j + 1
        if pv[j] > 1.e-15 and pv[j] < (1.0-1.e-15):
            SE -= pv[j]*log(pv[j], 2)
    return SE


def neumann(d, rho):
    b = lapak.zheevd(rho)
    VnE = shannon(d, b[0])
    return VnE


def mutual_info(rho, dl, dr):
    rhor = pTrace.pTraceL(dl, dr, rho)
    rhol = pTrace.pTraceR(dl, dr, rho)
    return neumann(dr, rhor) + neumann(dl, rhol) - neumann(dl*dr, rho)

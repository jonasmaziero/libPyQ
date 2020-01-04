import pTrace
import math
import numpy as np

def coh(rho, name):
    if name == 'l1':
        return coh_l1(rho)
    elif name == 're':
        return coh_re(rho)
    elif name == 'hs':
        return coh_hs(rho)

def coh_hs(rho):
    d = rho.shape[0]
    hsc = 0.0
    for j in range(0,d-1):
        for k in range(j+1,d):
            hsc += (rho[j][k].real)**2.0 + (rho[j][k].imag)**2.0
    return 2*hsc

def coh_l1(rho):  # normalized to [0,1]
    d = rho.shape[0]
    coh = 0.0
    for j in range(0, d-1):
        for k in range(j+1, d):
            coh += math.sqrt((rho[j][k].real)**2.0 + (rho[j][k].imag)**2.0)
    return 2.0*coh/(d-1)

def coh_re(rho):
    d = rho.shape[0]
    pv = np.zeros(d)
    for j in range(0,d):
        pv[j] = rho[j][j].real
    from entropy import shannon, von_neumann
    coh = shannon(pv) - von_neumann(rho)
    return coh/math.log(d,2)


def coh_nl(da, db, rho):
    rhoa = pTrace.pTraceL(da, db, rho)
    rhob = pTrace.pTraceR(da, db, rho)
    return coh_l1(da*db, rho)-coh_l1(da, rhoa)-coh_l1(db, rhob)

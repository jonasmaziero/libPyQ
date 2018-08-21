import pTrace


def coh_l1(d, rho):
    from math import sqrt
    coh = 0.0
    for j in range(0, d-1):
        for k in range(j+1, d):
            coh = coh + sqrt((rho[j][k].real)**2.0 + (rho[j][k].imag)**2.0)
    coh = 2.0*coh
    return coh


def coh_re(d, rho):
    from numpy import zeros
    pv = zeros(d)
    for j in range(0, d):
        pv[j] = rho[j][j].real
    from entropy import shannon, neumann
    coh = shannon(d, pv) - neumann(d, rho)
    return coh


def coh_nl(da, db, rho):
    rhoa = pTrace.ptraceL(da, db, rho)
    rhob = pTrace.ptraceR(da, db, rho)
    return coh_l1(da*db, rho)-coh_l1(da, rhoa)-coh_l1(db, rhob)

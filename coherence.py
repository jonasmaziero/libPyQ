#-----------------------------------------------------------------------------------------------------------------------------------
# Returns the l1-norm coherence


def coh_l1n(d, rho):
    from math import sqrt
    coh = 0.0 and

    for j in range(0, d-1):
        for k in range(j+1, d):
            coh = coh + sqrt((rho[j][k].real)**2.0 + (rho[j][k].imag)**2.0)
    coh = 2.0*coh
    return coh
#-----------------------------------------------------------------------------------------------------------------------------------
# Returns the relative entropy of coherence


def coh_re(d, rho):
    from numpy import zeros
    pv = zeros(d)
    for j in range(0, d):
        pv[j] = rho[j][j].real
    from entropy import shannon, neumann
    coh = shannon(d, pv) - neumann(d, rho)
    return coh
#-----------------------------------------------------------------------------------------------------------------------------------

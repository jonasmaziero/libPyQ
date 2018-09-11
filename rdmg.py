import numpy as np


def test():
    np.random.seed()
    import coherence as coh
    import matplotlib.pyplot as plt
    ns = 10**3  # number of samples for the average
    nqb = 5  # maximum number of qubits regarded
    Cavg = np.zeros(nqb)
    d = np.zeros(nqb, dtype=int)
    for j in range(0, nqb):
        d[j] = 2**(j+1)
        rdm = np.zeros((d[j], d[j]), dtype=complex)
        Cavg[j] = 0.0
        for k in range(0, ns):
            # rdm = rdm_ginibre(d[j])
            rdm = rdm_std(d[j])
            Cavg[j] += coh.coh_re(d[j], rdm)
            # Cavg[j] += coh.coh_l1(d[j], rdm)
        Cavg[j] = Cavg[j]/ns
        print(Cavg[j])
    plt.plot(d, Cavg, label='')
    plt.xlabel('d')
    plt.ylabel('C')
    plt.legend()
    plt.show()


def rdm_std(d):
    from rpvg import rpv_zhsl
    from rug import ru_gram_schmidt
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
    from distances import normHS
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

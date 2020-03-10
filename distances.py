import numpy as np


def fidelity_pp(psi, phi):
    psiD = np.conj(psi)
    return (abs(np.inner(psiD, phi)))**2


def fidelity_pm(psi, rho):
    psiD = np.conj(psi)
    phi = np.matmul(rho, psi)
    return abs(np.inner(psiD, phi))


def fidelity_mm(d, rho, zeta):
    from mat_func import mat_sqrt
    from pTrace import trace
    A = np.zeros((d, d), dtype=complex)
    B = np.zeros((d, d), dtype=complex)
    A = mat_sqrt(d, rho)
    B = np.matmul(np.matmul(A, zeta), A)
    A = mat_sqrt(d, B)
    return trace(d, A)


def norm(d, psi):
    N = 0
    for j in range(0, d):
        N += (psi[j].real)**2 + (psi[j].imag)**2
    return np.sqrt(N)


def normr(d, psi):
    N = 0
    for j in range(0, d):
        N += psi[j]**2
    return np.sqrt(N)


def inner(d, psi, phi):
    csi = np.zeros(d, dtype=complex)
    csi = np.conj(psi)
    ip = 0
    for j in range(0, d):
        ip += csi[j]*phi[j]
    return ip


def normTr(d, A):
    ev = np.linalg.eigvalsh(A)
    Ntr = 0
    for j in range(0, d):
        Ntr += abs(ev[j])
    return Ntr


def normHS(d, A):
    N = 0
    for j in range(0, d):
        for k in range(0, d):
            N += (A[j][k].real)**2.0 + (A[j][k].imag)**2.0
    return np.sqrt(N)

def normHS2(d,A):
    N2 = 0
    for j in range(0, d):
        for k in range(0, d):
            N2 += (A[j][k].real)**2 + (A[j][k].imag)**2
    return N2


'''
def test_distances():
    from numpy import zeros
    from math import sqrt

    # fidelity pure - pure
    z = zeros(2, dtype=complex)
    w = zeros(2, dtype=complex)
    z[0] = 1.0
    z[1] = 0.0
    w[0] = 1.0/sqrt(2.0)
    w[1] = 1.0/sqrt(2.0)
    print(fidelity_pp(z, w))

# fidelity pure - mixed
    x = zeros(2, dtype=complex)
    y = zeros(4, dtype=complex)
    x = [1.0/sqrt(2.0), 1.0/sqrt(2.0)*(1j)]
    y = [[1, 0], [0, 0]]
    print(fidelity_pm(x, y))

    # fidelity mixed - mixed
    x = zeros((2, 2), dtype=complex)
    y = zeros((2, 2), dtype=complex)
    x = [[1.0, 0.0], [0.0, 0.0]]
    y = [[1.0/2.0, 1.0/2.0], [1.0/2.0, 1.0/2.0]]
    print(fidelity_mm(2, x, y))


def fidelity_mm(d, rho, zeta):
    from numpy import linalg as LA
    val, evec = LA.eigh(rho)
    from numpy import kron, zeros
    from math import sqrt
    import mat_func as mf
    A = zeros((d, d), dtype=complex)
    phi = zeros(d, dtype=complex)
    psi = zeros(d, dtype=complex)
    for j in range(0, d):
        for l in range(0, d):
            phi[l] = evec[l][j]
        for k in range(0, d):
            for m in range(0, d):
                psi[m] = evec[m][k]
                A = A + (sqrt(val[j]*val[k])*mf.sandwich(d, phi, zeta, psi))
                *mf.outer(d, phi, psi)
    eigA = LA.eigvalsh(A)
    F = 0.0
    for n in range(0, d):
        F = F + sqrt(eigA[n])
    return F
'''

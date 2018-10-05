import numpy as np
from numpy import linalg as LA
from math import sqrt


def mat_sqrt(d, A):
    w, v = LA.eigh(A)
    Asr = np.zeros((d, d), dtype=complex)
    psi = np.zeros((d, 1), dtype=complex)
    for j in range(0, d):
        psi = v[:, j]
        Asr += sqrt(w[j])*proj(d, psi)
    return Asr


def adjunct(nr, nc, A):
    Aa = np.zeros((nc, nr), dtype=complex)
    for j in range(0, nr):
        for k in range(0, nc):
            Aa[k, j] = np.conj(A[j, k])
    return Aa


def transpose(nr, nc, A):
    At = np.zeros((nc, nr))
    for j in range(0, nr):
        for k in range(0, nc):
            At[k, j] = A[j, k]
    return At


def outer(d, psi, phi):
    op = np.zeros((d, d), dtype=complex)
    for j in range(0, d):
        for k in range(0, d):
            op[j, k] = psi[j]*np.conj(phi[k])
    return op


def outerr(d, psi, phi):
    op = np.zeros((d, d))
    for j in range(0, d):
        for k in range(0, d):
            op[j, k] = psi[j]*phi[k]
    return op


def proj(d, psi):
    return outer(d, psi, psi)


def sandwich(d, phi, A, psi):
    sd = 0
    for j in range(0, d):
        for k in range(0, d):
            sd += np.conj(phi[j])*A[j, k]*psi[k]
    return sd


'''
import su
A = np.dot(su.Pauli(1), su.Pauli(0))
w, v = LA.eigh(A)
d = 2
print(2, mat_sqrt(d, A))
'''

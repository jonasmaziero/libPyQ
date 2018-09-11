import numpy as np
from numpy import linalg as LA


def Dagger(A):
    return np.conjugate(np.transpose(A))


def transpose(m, n, A):
    At = np.zeros((n, m))
    for j in range(0, m):
        for k in range(0, n):
            At[k, j] = A[j, k]
    return At


def outer(psi, phi):
    return psi*np.conjugate(np.transpose(phi))


def outerr(psi, phi):
    return psi*np.transpose(phi)


def proj(psi):
    return np.outer(psi, np.conj(psi))


def sandwich(phi, A, psi):
    return np.dot(np.dot(Dagger(phi), A), phi)


def mat_sqrt(d, A):
    w, v = LA.eigh(A)
    Asr = np.zeros((d, d), dtype=complex)
#    Proj = np.zeros((d,d),dtype=complex)
    for j in range(0, d):
        Proj = proj(v[:, j])
        Asr += (np.sqrt(w[j])/np.trace(Proj))*Proj
    return Asr


'''
import su
A = np.dot(su.Pauli(1), su.Pauli(0))
w, v = LA.eigh(A)
d = 2
print(2, mat_sqrt(d, A))
'''

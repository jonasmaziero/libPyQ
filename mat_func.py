import numpy as np
from numpy import linalg as LA
import math


def mat_sqrt(d, A):
    w, v = LA.eigh(A)
    Asr = np.zeros((d, d), dtype=complex)
    psi = np.zeros((d, 1), dtype=complex)
    for j in range(0, d):
        psi = v[:, j]
        Asr += math.sqrt(abs(w[j]))*proj(d, psi)
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

def ip_c(d,v,w):
    ipc = 0.0 + (1j)*0.0
    for j in range(0,d):
        ipc += np.conj(v[j])*w[j]
    return ipc

def vnorm_c(d,v):
    return math.sqrt(ip_c(d,v,v).real)

def vnorm2_c(d, v):
    vn = 0.0
    for j in range(0,d):
        vn += (v[j].real)**2 + (v[j].imag)**2
    return vn

def versor_c(d,v):
    return v/vnorm_c(d,v)

'''
import su
A = np.dot(su.Pauli(1), su.Pauli(0))
w, v = LA.eigh(A)
d = 2
print(2, mat_sqrt(d, A))
'''

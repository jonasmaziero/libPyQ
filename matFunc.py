import numpy as np


def Dagger(A):
    return np.conjugate(np.transpose(A))


def proj(psi):
    return np.dot(psi, Dagger(psi))


def sandwich(phi, A, psi):
    return np.dot(np.dot(Dagger(phi), A), phi)

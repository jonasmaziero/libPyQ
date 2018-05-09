from math import sqrt
import numpy as np
from su import Pauli
from matFunc import proj  # , Dagger


def Bell(j, k):
    if j == 0 and k == 0:
        return np.array([[1/sqrt(2)], [0], [0], [1/sqrt(2)]])  # phi+
    elif j == 0 and k == 1:
        return np.array([[0], [1/sqrt(2)], [1/sqrt(2)], [0]])  # psi+
    elif j == 1 and k == 0:
        return np.array([[1/sqrt(2)], [0], [0], [-1/sqrt(2)]])  # phi-
    elif j == 1 and k == 1:
        return np.array([[0], [1/sqrt(2)], [-1/sqrt(2)], [0]])  # psi-


def Werner(w):
    return ((1-w)/4)*np.eye(4) + w*proj(Bell(1, 1))


def rhoBD(c1, c2, c3):
    return (1/4)*(np.eye(4) + c1*np.kron(Pauli(1), Pauli(1))
                  + c2*np.kron(Pauli(2), Pauli(2)) + c3*np.kron(Pauli(3), Pauli(3)))


def cb(d, j):
    cbs = np.zeros(d)
    cbs[j] = 1
    return cbs


def psi1qb(theta, phi):
    return np.cos(theta/2.0)*cb(2, 0) + (np.cos(phi)
                                         + np.sin(phi)*1j)*np.sin(theta/2.0)*cb(2, 1)


def rho1qb(r1, r2, r3):
    return 0.5*(Pauli(0) + r1*Pauli(1) + r2*Pauli(2) + r3*Pauli(3))


def rho2qb(CM):
    return (1/4)*(CM[0][0]*np.kron(Pauli(0), Pauli(0)) + CM[1][0]*np.kron(Pauli(1), Pauli(0))
                  + CM[2][0]*np.kron(Pauli(2), Pauli(0)) + CM[3][0]*np.kron(Pauli(3), Pauli(0))
                  + CM[0][1]*np.kron(Pauli(0), Pauli(1)) + CM[1][1]*np.kron(Pauli(1), Pauli(1))
                  + CM[2][1]*np.kron(Pauli(2), Pauli(1)) + CM[3][1]*np.kron(Pauli(3), Pauli(1))
                  + CM[0][2]*np.kron(Pauli(0), Pauli(2)) + CM[1][2]*np.kron(Pauli(1), Pauli(2))
                  + CM[2][2]*np.kron(Pauli(2), Pauli(2)) + CM[3][2]*np.kron(Pauli(3), Pauli(2))
                  + CM[0][3]*np.kron(Pauli(0), Pauli(3)) + CM[1][3]*np.kron(Pauli(1), Pauli(3))
                  + CM[2][3]*np.kron(Pauli(2), Pauli(3)) + CM[3][3]*np.kron(Pauli(3), Pauli(3)))


#print(Dagger(Bell(1, 1)))
#print(rho1qb(0, 0, 0))

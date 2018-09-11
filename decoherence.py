import numpy as np
from su import Pauli


def Kpd(j, p):
    if j == 0:
        return np.array([[0, 0], [0, np.sqrt(p)]])
    elif j == 1:
        return np.array([[1, 0], [0, np.sqrt(1-p)]])


def Kgad(j, g, a):
    if j == 0:
        return np.sqrt(g)*np.array([[1, 0], [0, np.sqrt(1-a)]])
    elif j == 1:
        return np.sqrt(g)*np.array([[0, np.sqrt(a)], [0, 0]])
    elif j == 2:
        return np.sqrt(1-g)*np.array([[np.sqrt(1-a), 0], [0, 1]])
    elif j == 3:
        return np.sqrt(1-g)*np.array([[0, 0], [np.sqrt(a), 0]])


def Kpdgad(j, p, g, a):
    if j == 0:
        return np.matmul(Kpd(0, p), Kgad(0, g, a))
    elif j == 1:
        return np.matmul(Kpd(0, p), Kgad(1, g, a))
    elif j == 2:
        return np.matmul(Kpd(0, p), Kgad(2, g, a))
    elif j == 3:
        return np.matmul(Kpd(0, p), Kgad(3, g, a))
    elif j == 4:
        return np.matmul(Kpd(1, p), Kgad(0, g, a))
    elif j == 5:
        return np.matmul(Kpd(1, p), Kgad(1, g, a))
    elif j == 6:
        return np.matmul(Kpd(1, p), Kgad(2, g, a))
    elif j == 7:
        return np.matmul(Kpd(1, p), Kgad(3, g, a))


def werner_pdad(w, p, a):
    return np.array([[(1-w+a+a*w)/4, 0, 0, 0],
                     [0, (1+w+a-a*w)/4, -(w*np.sqrt(1-p)*np.sqrt(1-a))/2, 0],
                     [0, -(w*np.sqrt(1-p)*np.sqrt(1-a))/2, (1+w-a-a*w)/4, 0],
                     [0, 0, 0, (1-w-a+a*w)/4]])

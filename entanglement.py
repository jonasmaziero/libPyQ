import numpy as np
from numpy import linalg as LA
from su import Pauli
from math import sqrt


def test_entanglement():
    #from pTranspose import Ta
    Ec = np.zeros(100)
    #EF = zeros(100)
    #En = zeros(100)
    #Eln = zeros(100)
    x = np.zeros(100)
    from states import Werner
    dw = 1.01/100
    w = -dw
    for j in range(0, 100):
        w = w + dw
        if w > 1.0:
            break
        rho = Werner(w)
        Ec[j] = concurrence(rho)
        #EF[j] = EoF(rho)
        #rhoTp = Ta(2, 2, rho)
        #En[j] = negativity(4, rhoTp)
        #Eln[j] = log_negativity(4, rhoTp)
        x[j] = w
    import matplotlib.pyplot as plt
    plt.plot(x, Ec, label='Ec')
    #plt.plot(x, EF, label='EoF')
    #plt.plot(x, En, label='En')
    #plt.plot(x, Eln, label='Eln')
    plt.xlabel('x')
    plt.ylabel('')
    plt.legend(loc=4)
    plt.show()


def concurrence(rho):
    ev = np.zeros(4, dtype='float')
    R = np.zeros((4, 4), dtype=complex)
    R = np.dot(rho, np.kron(Pauli(2), Pauli(2)))
    R = np.dot(R, np.conj(rho))
    R = np.dot(R, np.kron(Pauli(2), Pauli(2)))
    ev = LA.eigvalsh(R)
    evm = max(abs(ev[0]), abs(ev[1]), abs(ev[2]), abs(ev[3]))
    C = 2.0*sqrt(abs(evm)) - sqrt(abs(ev[0]))
    C = C - sqrt(abs(ev[1])) - sqrt(abs(ev[2])) - sqrt(abs(ev[3]))
    if C < 0.0:
        C = 0.0
    return C


'''
def concurrence(rho):
    s2s2 = np.zeros((4, 4), dtype=complex)
    rhoc = np.zeros((4, 4), dtype=complex)
    rhoa = np.zeros((4, 4), dtype=complex)
    rhot = np.zeros((4, 4), dtype=complex)
    R = np.zeros((4, 4), dtype=complex)
    ev = np.zeros(4)
    s2s2 = [[0.0, 0.0, 0.0, -1.0], [0.0, 0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 0.0]]
    rhoc = np.conjugate(rho)
    rhoa = np.dot(s2s2, rhoc)
    rhot = np.dot(rhoa, s2s2)
    R = np.dot(rho, rhot)
    ev = LA.eigvalsh(R)
    evm = max(ev[0], ev[1], ev[2], ev[3])
    C = 2.0*sqrt(abs(evm)) - sqrt(abs(ev[0])) - \
        sqrt(abs(ev[1])) - sqrt(abs(ev[2])) - sqrt(abs(ev[3]))
    if C < 0.0:
        C = 0.0
    return C
'''

'''
def EoF(rho):
    from numpy import zeros
    pv = zeros(2)
    from math import sqrt
    Ec = concurrence(rho)
    pv[0] = (1.0 + sqrt(1.0 - Ec**2.0))/2.0
    pv[1] = 1.0 - pv[0]
    from entropy import shannon
    EF = shannon(2, pv)
    return EF


def negativity(d, rhoTp):
    from distances import normTr
    En = 0.5*(normTr(d, rhoTp) - 1.0)
    return En


def log_negativity(d, rhoTp):
    from math import log
    En = negativity(d, rhoTp)
    Eln = log(2.0*En+1.0, 2)
    return Eln
'''

# test_entanglement()

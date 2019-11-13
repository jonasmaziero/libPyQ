import numpy as np
import matplotlib.pyplot as plt
import rdmg
import gell_mann as gm
from math import log
from mat_func import mat_sqrt
from pTrace import trace
import platform
import entropy


def coh_pop():
    d = 4
    rho = np.zeros((d, d), dtype=complex)
    rhosr = np.zeros((d, d), dtype=complex)
    bv = np.zeros(d**2-1)
    ns = 5*10**3
    coh = np.zeros(ns)
    # sl = np.zeros(ns)
    svn = np.zeros(ns)
    # maxx = 0
    for j in range(0, ns):
        rho = rdmg.rdm_std(d)
        # bv = gm.bloch_vector(d, rho)
        rhosr = mat_sqrt(d, rho)
        bv = gm.bloch_vector(d, rhosr)
        # coh[j] = coh_hs(d, bv)
        coh[j] = coh_he(d, bv)
        # sl[j] = Sl(d, bv)
        trsr = trace(d, rhosr)
        # sl[j] = SlHe(d, bv, trsr)
        svn[j] = SvnHe(d, bv, trsr)
    x = np.zeros(11)
    y = np.zeros(11)
    dx = (1-1/d)/10
    # dx = maxx/10
    # dx = ((d-1)/2)/10
    # dx = log(d)/10
    xx = -dx
    for j in range(0, 11):
        xx += dx
        x[j] = xx
        y[j] = x[j]
    plt.plot(x, y, color='black', label=r'$d=4$')
    # plt.scatter(sl, coh, color='blue', s=2, marker='.')
    # plt.scatter(sl, coh, color='blue', s=2, marker='.')
    plt.scatter(svn, coh, color='blue', s=2, marker='.')
    # plt.xlabel(r'$S_{l}$')
    plt.xlabel(r'$\Omega$')
    plt.ylabel(r'$C_{he}$')
    plt.legend()
    if platform.system() == 'Linux':
        path1 = '/home/jonasmaziero/Dropbox/Research/qnesses/coherence/'
        path = path1 + 'coeh_pop_tradeoff/calc/hevnd4.eps'
        plt.savefig(path, format='eps', dpi=100)
    else:
        path1 = '/Users/jonasmaziero/Dropbox/Research/qnesses/coherence/'
        path = path1 + 'coeh_pop_tradeoff/calc/hsd8.eps'
        plt.savefig(path, format='eps', dpi=100)
    plt.show()


def coh_hs(d, bv):
    coh = 0
    for j in range(d-1, d**2-1):
        coh += bv[j]**2
    return coh/2


def Sl(d, bv):
    s = 0
    for j in range(0, d-1):
        s -= bv[j]**2
    s /= 2
    s += (d-1)/d
    return s


def Svn(d, bv):
    nbv = np.zeros(d**2-1)
    for j in range(0, d-1):
        nbv[j] = bv[j]
    for j in range(d-1, d**2-1):
        nbv[j] = 0
    dm = gm.rho(d, nbv)
    pd = np.zeros(d)
    for j in range(0, d):
        pd[j] = dm[j, j].real
    return entropy.shannon(pd)


def coh_he(d, bvsr):
    coh = 0
    for j in range(d-1, d**2-1):
        coh += bvsr[j]**2
    return coh/2


def SlHe(d, bvsr, trsr):
    s = 0
    for j in range(0, d-1):
        s -= bvsr[j]**2
    s /= 2
    s += ((trsr)**2)*(1-1/d)
    return s


def SvnHe(d, bvsr, trsr):
    nbv = np.zeros(d**2-1)
    for j in range(0, d-1):
        nbv[j] = bvsr[j]
    for j in range(d-1, d**2-1):
        nbv[j] = 0
    dm = gm.rho(d, nbv)
    pd = np.zeros(d)
    for j in range(0, d):
        pd[j] = dm[j, j].real
    svn = entropy.shannon(d, pd) + trsr*(trsr-1)
    return svn

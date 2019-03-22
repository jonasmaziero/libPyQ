import numpy as np
from math import sin, cos
from entropy import mutual_info
from constants import pi
from numpy import linalg as LA
from math import sqrt
import gell_mann as gm
from mat_func import mat_sqrt, transpose, outerr
from distances import normr
from pTrace import trace, pTraceL, pTraceR


def hellinger(da, db, rho):  # arXiv:1510.06995
    daa = da**2-1
    M = mat_sqrt(da*db, rho)
    A = pTraceR(da, db, M)
    bva = gm.bloch_vector(da, A)/sqrt(2*db)
    B = pTraceL(da, db, M)
    bvb = gm.bloch_vector(db, B)/2
    cm = gm.corr_mat(da, db, M)/2
    ev = np.zeros(3)
    ev = LA.eigvalsh(outerr(daa, bva, bva)+np.matmul(cm, transpose(3, 3, cm)))
    mev = max(ev[0], ev[1], ev[2])
    bvbn = normr(2, bvb)
    nor = 1-1/sqrt(da)
    return max(0, (1-sqrt((trace(da, A)/sqrt(2*db))**2+bvbn**2+mev))/nor)


def oz_2qb(rho):
    return mutual_info(rho, 2, 2) - ccorr_hv_2qb(rho)


def ccorr_hv_2qb(rho):
    st = 0.01
    cc = 0
    th = -st
    while th < pi():
        th += st
        ph = -st
        while ph < 2*pi():
            ph += st
            mi = mutual_info(rhoMb(rho, th, ph), 2, 2)
            if mi > cc:
                cc = mi
    return cc


def rhoMb(rho, th, ph):
    rm = np.matmul(np.matmul(P1(th, ph), rho), P1(th, ph))
    rm += np.matmul(np.matmul(P2(th, ph), rho), P2(th, ph))
    return rm


def P1(th, ph):
    ep = cos(ph) + 1j*sin(ph)
    en = cos(ph) - 1j*sin(ph)
    return np.array([[cos(th/2)**2, en*sin(th/2)*cos(th/2), 0, 0],
                     [ep*sin(th/2)*cos(th/2), sin(th/2)**2, 0, 0],
                     [0, 0, cos(th/2)**2, en*sin(th/2)*cos(th/2)],
                     [0, 0, ep*sin(th/2)*cos(th/2), sin(th/2)**2]])


def P2(th, ph):
    ep = cos(ph) + 1j*sin(ph)
    en = cos(ph) - 1j*sin(ph)
    return np.array([[sin(th/2)**2, -en*sin(th/2)*cos(th/2), 0, 0],
                     [-ep*sin(th/2)*cos(th/2), cos(th/2)**2, 0, 0],
                     [0, 0, sin(th/2)**2, -en*sin(th/2)*cos(th/2)],
                     [0, 0, -ep*sin(th/2)*cos(th/2), cos(th/2)**2]])


def discord_oz_bds(rho):
    # Returns the OLLIVIER-ZUREK discord for 2-qubit Bell-diagonal states
    D = mi_bds(rho) - ccorr_hv_bds(rho)
    return D


def ccorr_hv_bds(rho):
    # Returns the Henderson-Vedral classical correlation for 2-qubit
    # Bell-diagonal states
    from math import log
    c3 = 4.0*rho[0][0] - 1.0
    c1 = 2.0*(rho[0][3] + rho[1][2])
    c2 = 2.0*(rho[1][2] - rho[0][3])
    c = max(abs(c1), abs(c2), abs(c3))
    cc = float(0.5*((1.0-c)*log(1.0-c, 2) + (1.0+c)*log(1.0+c, 2)))
    return cc


def mi_bds(rho):
    # Returns the mutual information for 2-qubit Bell-diagonal states
    from math import log
    c3 = 4.0*rho[0][0] - 1.0
    c1 = 2.0*(rho[0][3] + rho[1][2])
    c2 = 2.0*(rho[1][2] - rho[0][3])
    l00 = (1.0 + c1 - c2 + c3)/4.0
    l01 = (1.0 + c1 + c2 - c3)/4.0
    l10 = (1.0 - c1 + c2 + c3)/4.0
    l11 = (1.0 - c1 - c2 - c3)/4.0
    mi = float(l00*log(4.0*l00, 2) + l01*log(4.0*l01, 2) +
               l10*log(4.0*l10, 2) + l11*log(4.0*l11, 2))
    return mi

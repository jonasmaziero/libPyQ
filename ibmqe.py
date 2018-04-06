import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg.lapack as lapak
import numpy as np


def bdsCorr():
    from numpy import zeros
    from math import acos, sqrt
    c1 = -0.9
    c2 = -0.8
    c3 = -0.8
    from states import rho_bds
    rho = zeros((4, 4), dtype=complex)
    rho = rho_bds(c1, c2, c3)
    eig = lapak.zheevd(rho)
    print('eigens:', eig[0][0], eig[0][1], eig[0][2], eig[0][3])
    import discord
    dp = 0.05
    d = 20*(1-0)
    di = zeros(d)
    cc = zeros(d)
    mi = zeros(d)
    diE = zeros(d)
    ccE = zeros(d)
    miE = zeros(d)
    rhopE = zeros((4, 4), dtype=complex)
    pv = zeros(d)
    rhop = zeros((4, 4), dtype=complex)
    import kraus
    import tomography as tomo
    p = -dp
    for j in range(0, d):
        p = p + dp
        # rhop = kraus.rho_pd(p, rho)  # ; print(rhop) # teste para ver a cara do rhop
        c1p = c1*(1.0-p)**2
        c2p = c2*(1.0-p)**2
        c3p = c3
        print("p = ", p, ", c1p=", c1p, ", c2p=", c2p, ", c3p=", c3p)
        p00 = (1.0 + c1p - c2p + c3p)/4.0
        p01 = (1.0 + c1p + c2p - c3p)/4.0
        p10 = (1.0 - c1p + c2p + c3p)/4.0
        p11 = (1.0 - c1p - c2p - c3p)/4.0
        print("p00 = ", p00, ", p01=", p01, ", p10=", p10, ", p11=", p11)
        theta = 2.0*acos(sqrt(p00 + p01))
        alpha = 2.0*acos(sqrt(p00 + p10))
        print("theta = ", theta, ", alpha", alpha)
        print("")
        rhop = rho_bds(c1p, c2p, c3p)
        pv[j] = p
        di[j] = discord.discord_oz_bds(rhop)
        cc[j] = discord.ccorr_hv_bds(rhop)
        mi[j] = discord.mi_bds(rhop)
        '''sj = str(j)
        path1 = '/home/jonasmaziero/Dropbox/Research/IBM_QC/'
        path2 = 'tomography/BDS/bell_diagonal/'
        path = path1 + path2 + sj + '/'
        rhopE = tomo.tomo_2qb(path)
        diE[j] = discord.discord_oz_bds(rhopE)
        ccE[j] = discord.ccorr_hv_bds(rhopE)
        miE[j] = discord.mi_bds(rhopE)
        eigE = lapak.zheevd(rhopE)
        print('eigensE:', eigE[0][0], eigE[0][1], eigE[0][2], eigE[0][3])
        print('di =', di[j], '  diE = ', diE[j])'''

    plt.plot(pv, di, label='di')
    plt.plot(pv, cc, label='cc')
    plt.plot(pv, mi, label='im')
    plt.xlabel('p')
    plt.legend()
    plt.show()


def Werner():
    from numpy import zeros
    import tomography as tomo
    import pTranspose as Tp
    rhoE = zeros((4, 4), dtype=complex)
    import entanglement as ent
    N = 14
    wv = np.array([0, 0.1, 0.2, 0.3, 0.32,
                   0.34, 0.36, 0.38, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    Ev = zeros(N)
    for j in range(0, N):
        sj = str(j)
        path1 = '/home/jonasmaziero/Dropbox/Research/IBM_QC/'
        path2 = 'tomography/Werner/'
        path = path1 + path2 + sj + '/'
        rhoE = tomo.tomo_2qb(path)
        #Ev[j] = ent.concurrence(rhoE)
        Ev[j] = 2.*ent.negativity(4, Tp.pTransposeL(2, 2, rhoE))
        '''if j == 13:
            print(wv[j], Ev[j])
            tomo.plot_rho2qb(rhoE)'''
    plt.plot(wv, Ev, label='E')
    plt.xlabel('w')
    plt.legend()
    plt.show()

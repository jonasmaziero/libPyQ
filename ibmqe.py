import matplotlib.pyplot as plt
import numpy as np
#import scipy.linalg.lapack as lapak
from mpl_toolkits.mplot3d import Axes3D
import pTranspose as pT
import discord


def werner():
    import tomography as tomo
    import coherence as coh
    import entanglement as ent
    from distances import fidelity_mm
    from states import Werner
    Ne = 11
#    we = np.array([0, 0.1, 0.2, 0.3, 0.32,
#                   0.34, 0.36, 0.38, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    we = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    Ee = np.zeros(Ne)
    Cnle = np.zeros(Ne)
    Nle = np.zeros(Ne)
    Se = np.zeros(Ne)
    De = np.zeros(Ne)
    F = np.zeros(Ne)
    for j in range(0, Ne):
        sj = str(j)
        path1 = '/home/jonasmaziero/Dropbox/Research/ibm/bds/'
        path2 = 'werner_qx2/dados_plot/'
        path = path1 + path2 + sj + '/'
        rhoe = tomo.tomo_2qb(path)
#        Ee[j] = ent.concurrence(rhoe)
        Ee[j] = 2*ent.negativity(4, pT.pTransposeL(2, 2, rhoe))
        Cnle[j] = coh.coh_nl(2, 2, rhoe)
        Nle[j] = ent.chsh(rhoe)
        Se[j] = ent.steering(rhoe)
#        De[j] = discord.hellinger(2, 2, rhoe)
        De[j] = discord.oz_2qb(rhoe)
        F[j] = fidelity_mm(4, Werner(j*0.1), rhoe)
    Nt = 100
    Et = np.zeros(Nt)
    Nlt = np.zeros(Nt)
    St = np.zeros(Nt)
    Cnlt = np.zeros(Nt)
    wt = np.zeros(Nt)
    Dt = np.zeros(Nt)
    dw = 1.01/Nt
    w = -dw
    for j in range(0, Nt):
        w = w + dw
        if w > 1.0:
            break
        rho = Werner(w)
#        Et[j] = ent.concurrence(rho)
        Et[j] = 2*ent.negativity(4, pT.pTransposeL(2, 2, rho))
        Cnlt[j] = coh.coh_nl(2, 2, rho)
        Nlt[j] = ent.chsh(rho)
        St[j] = ent.steering(rho)
#        Dt[j] = discord.hellinger(2, 2, rho)
        Dt[j] = discord.oz_2qb(rho)
        wt[j] = w
    plt.plot(wt, Cnlt, '.', label='$C$', color='gray')
    plt.plot(we, Cnle, '*', label=r'$C_{e}$', color='gray')
    plt.plot(wt, Dt, '-', label='D', color='magenta')
    plt.plot(we, De, 'o', label=r'$D_{e}$', color='magenta')
    plt.plot(wt, Et, '-.', label='E', color='blue')
    plt.plot(we, Ee, 's', label=r'$E_{e}$', color='blue')
    plt.plot(wt, St, ':', label='$S$', color='red')
    plt.plot(we, Se, '^', label=r'$S_{e}$', color='red')
    plt.plot(wt, Nlt, '--', label='$N$', color='cyan')
    plt.plot(we, Nle, 'h', label=r'$N_{e}$', color='cyan')
    plt.plot(we, F, 'X', label=r'$F$', color='black')
    plt.xlabel('w')
    plt.legend(loc=6)
    import platform
    if platform.system() == 'Linux':
        plt.savefig('/home/jonasmaziero/Dropbox/Research/ibm/bds/calc/qcorr.eps',
                    format='eps', dpi=100)
    else:
        plt.savefig('/Users/jonasmaziero/Dropbox/Research/ibm/bds/calc/qcorr.eps',
                    format='eps', dpi=100)
    plt.show()


def werner_decoh():
    import coherence as coh
    import entanglement as ent
    from distances import fidelity_mm
    from states import Werner
    from decoherence import werner_pdad
    Nt = 100
    wt = np.zeros(Nt)
    Et = np.zeros(Nt)
    Nlt = np.zeros(Nt)
    St = np.zeros(Nt)
    Cnlt = np.zeros(Nt)
    Dt = np.zeros(Nt)
    Etd = np.zeros(Nt)
    Nltd = np.zeros(Nt)
    Std = np.zeros(Nt)
    Cnltd = np.zeros(Nt)
    Dtd = np.zeros(Nt)
    p = 0.3
    a = p
    dw = 1.01/Nt
    w = -dw
    for j in range(0, Nt):
        w = w + dw
        if w > 1.0:
            break
        rho = Werner(w)
#        Et[j] = ent.concurrence(rho)
        Et[j] = 2*ent.negativity(4, pT.pTransposeL(2, 2, rho))
        Cnlt[j] = coh.coh_nl(2, 2, rho)
        Nlt[j] = ent.chsh(rho)
        St[j] = ent.steering(rho)
#        Dt[j] = discord.hellinger(2, 2, rho)
        Dt[j] = discord.oz_2qb(rho)
        wt[j] = w
        rhod = werner_pdad(w, p, a)
#        Etd[j] = ent.concurrence(rhod)
        Etd[j] = 2*ent.negativity(4, pT.pTransposeL(2, 2, rhod))
        Cnltd[j] = coh.coh_nl(2, 2, rhod)
        Nltd[j] = ent.chsh(rhod)
        Std[j] = ent.steering(rhod)
#        Dtd[j] = discord.hellinger(2, 2, rhod)
        Dtd[j] = discord.oz_2qb(rhod)
    plt.plot(wt, Cnlt, '.', label='$C$', color='gray')
    plt.plot(wt, Cnltd, '*', markersize=4, label=r'$C_{d}$', color='gray')
    plt.plot(wt, Dt, '-', label='D', color='magenta')
    plt.plot(wt, Dtd, 'o', markersize=4, label=r'$D_{d}$', color='magenta')
    plt.plot(wt, Et, '-.', label='E', color='blue')
    plt.plot(wt, Etd, 's', markersize=4, label=r'$E_{d}$', color='blue')
    plt.plot(wt, St, ':', label='$S$', color='red')
    plt.plot(wt, Std, '^', markersize=4, label=r'$S_{d}$', color='red')
    plt.plot(wt, Nlt, '--', label='$N$', color='cyan')
    plt.plot(wt, Nltd, 'h', markersize=4, label=r'$N_{d}$', color='cyan')
    plt.xlabel('w')
    plt.legend(loc=6)
    plt.show()


'''
    import platform
    if platform.system() == 'Linux':
        plt.savefig('/home/jonasmaziero/Dropbox/Research/ibm/bds/calc/qcorr_decoh025.eps',
                    format='eps', dpi=100)
    else:
        plt.savefig('/Users/jonasmaziero/Dropbox/Research/ibm/bds/calc/qcorr_decoh099.eps',
                    format='eps', dpi=100)
    plt.show()
'''
# --------------
'''
    def bdsCorr():
        from math import acos, sqrt
        c1 = -0.9
        c2 = -0.8
        c3 = -0.8
        from states import rhoBD
        rho = np.zeros((4, 4), dtype=complex)
        rho = rhoBD(c1, c2, c3)
        eig = lapak.zheevd(rho)
        print('eigens:', eig[0][0], eig[0][1], eig[0][2], eig[0][3])
        import discord
        dp = 0.05
        d = 20*(1-0)
        di = np.zeros(d)
        cc = np.zeros(d)
        mi = np.zeros(d)
    #    diE = np.zeros(d)
    #    ccE = np.zeros(d)
    #    miE = np.zeros(d)
    #    rhopE = np.zeros((4, 4), dtype=complex)
        pv = np.zeros(d)
        rhop = np.zeros((4, 4), dtype=complex)
    #    import kraus
    #    import tomography as tomo
        p = -dp
        for j in range(0, d):
            p = p + dp
    #        rhop = kraus.rho_pd(p, rho) # ; print(rhop) # teste para ver rhop
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
            rhop = rhoBD(c1p, c2p, c3p)
            pv[j] = p
            di[j] = discord.discord_oz_bds(rhop)
            cc[j] = discord.ccorr_hv_bds(rhop)
            mi[j] = discord.mi_bds(rhop)
            sj = str(j)
            path1 = '/home/jonasmaziero/Dropbox/Research/IBM_QC/'
            path2 = 'tomography/BDS/bell_diagonal/'
            path = path1 + path2 + sj + '/'
            rhopE = tomo.tomo_2qb(path)
            diE[j] = discord.discord_oz_bds(rhopE)
            ccE[j] = discord.ccorr_hv_bds(rhopE)
            miE[j] = discord.mi_bds(rhopE)
            eigE = lapak.zheevd(rhopE)
            print('eigensE:', eigE[0][0], eigE[0][1], eigE[0][2], eigE[0][3])
            print('di =', di[j], '  diE = ', diE[j])
        plt.plot(pv, di, label='di')
        plt.plot(pv, cc, label='cc')
        plt.plot(pv, mi, label='im')
        plt.xlabel('p')
        plt.legend()
        plt.show()
'''

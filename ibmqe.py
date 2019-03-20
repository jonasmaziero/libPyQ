import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pTranspose as pT
import discord
import coherence as coh
import entanglement as ent
from distances import fidelity_mm
from states import Werner
from decoherence import werner_pdad
import tomography as tomo
from math import sqrt


def werner():
    Nw = 11 # no. of experiments of each configuration
    Nr = 5 # no. of rounds of the experiment
    we = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    Ee = np.zeros(Nw)
    Eerr = np.zeros(Nw)  # for the standard deviation
    Cnle = np.zeros(Nw)
    Cnlerr = np.zeros(Nw)
    Nle = np.zeros(Nw)
    Nlerr = np.zeros(Nw)
    Se = np.zeros(Nw)
    Serr = np.zeros(Nw)
    De = np.zeros(Nw)
    Derr = np.zeros(Nw)
    Fe = np.zeros(Nw)
    Ferr = np.zeros(Nw)
    for j in range(0, Nw):
        sj = str(j)
        Em = 0.0;  E2m = 0.0;  Fm = 0.0;  F2m = 0.0;  Cnlm = 0.0;  Cnl2m = 0.0
        Nlm = 0.0;  Nl2m = 0.0;  Sm = 0.0;  S2m = 0.0;  Dm = 0.0;  D2m = 0.0
        for k in range(0,Nr):
            sk = str(k)
            path1 = '/home/jonas/Dropbox/Research/ibm/bds'
            #path1 = '/Users/jonasmaziero/Dropbox/Research/ibm/bds'
            path2 = '/calc_mauro/dados_plot/dados_plot'
            path = path1 + path2 + sk + '/' + sj + '/'
            if k == 0:
                rhoe = tomo.tomo_2qb(path)
            else:
                rhoe = tomo.tomo_2qb_(path)
            E = 2*ent.negativity(4, pT.pTransposeL(2, 2, rhoe))
            Em += E;  E2m += pow(E,2)
            F = fidelity_mm(4, Werner(j*0.1), rhoe);  Fm += F;  F2m += pow(F,2)
            Cnl = coh.coh_nl(2, 2, rhoe);  Cnlm += Cnl;  Cnl2m += pow(Cnl,2)
            Nl = ent.chsh(rhoe);  Nlm += Nl;  Nl2m += pow(Nl,2)
            S = ent.steering(rhoe);  Sm += S;  S2m += pow(S,2)
            #D = discord.oz_2qb(rhoe);  Dm += D;  D2m += pow(D,2)
        Em = Em/Nr;  E2m = E2m/Nr;  Eerr[j] = sqrt(E2m - pow(Em,2));  Ee[j] = Em
        Fm = Fm/Nr;  F2m = F2m/Nr;  Ferr[j] = sqrt(F2m - pow(Fm,2));  Fe[j] = Fm
        Cnlm = Cnlm/Nr;  Cnl2m = Cnl2m/Nr;  Cnlerr[j] = sqrt(Cnl2m - pow(Cnlm,2));  Cnle[j] = Cnlm
        Nlm = Nlm/Nr;  Nl2m = Nl2m/Nr;  Nlerr[j] = sqrt(Nl2m - pow(Nlm,2));  Nle[j] = Nlm
        Sm = Sm/Nr;  S2m = S2m/Nr;  Serr[j] = sqrt(S2m - pow(Sm,2));  Se[j] = Sm
        Dm = Dm/Nr;  D2m = D2m/Nr;  Derr[j] = sqrt(D2m - pow(Dm,2));  De[j] = Dm
        #F = fidelity_mm(4, Werner(j*0.1), rhoe)
        #Ee = 2*ent.negativity(4, pT.pTransposeL(2, 2, rhoe))
        #Cnle = coh.coh_nl(2, 2, rhoe)
        #Nle[j] = ent.chsh(rhoe)
        #Se[j] = ent.steering(rhoe)
        # De[j] = discord.hellinger(2, 2, rhoe)
        #De[j] = discord.oz_2qb(rhoe)
        #F[j] = fidelity_mm(4, Werner(j*0.1), rhoe)
    Nt = 110
    Et = np.zeros(Nt)
    Nlt = np.zeros(Nt)
    St = np.zeros(Nt)
    Cnlt = np.zeros(Nt)
    wt = np.zeros(Nt)
    Dt = np.zeros(Nt)
    Etd = np.zeros(Nt)
    Nltd = np.zeros(Nt)
    Std = np.zeros(Nt)
    Cnltd = np.zeros(Nt)
    Dtd = np.zeros(Nt)
    p = 0.2
    a = p
    dw = 1.01/Nt
    w = -dw
    for j in range(0, Nt):
        w = w + dw
        if w > 1.01:
            break
        rho = Werner(w)
        Et[j] = 2*ent.negativity(4, pT.pTransposeL(2, 2, rho))
        Cnlt[j] = coh.coh_nl(2, 2, rho)
        Nlt[j] = ent.chsh(rho)
        St[j] = ent.steering(rho)
        #Dt[j] = discord.oz_2qb(rho)
        # Dt[j] = discord.hellinger(2, 2, rho)
        #rhod = werner_pdad(w, p, a)
        # Etd[j] = ent.concurrence(rhod)
        #Etd[j] = 2*ent.negativity(4, pT.pTransposeL(2, 2, rhod))
        #Cnltd[j] = coh.coh_nl(2, 2, rhod)
        #Nltd[j] = ent.chsh(rhod)
        #Std[j] = ent.steering(rhod)
        # Dtd[j] = discord.hellinger(2, 2, rhod)
        #Dtd[j] = discord.oz_2qb(rhod)
        wt[j] = w
    plt.errorbar(we, Fe, Ferr, marker='x', label=r'$F$', color='black', markersize=5)
    #plt.plot(we, F, 'x', label=r'$F$', color='black')
    plt.plot(wt, Cnlt, '.', label='$C$', color='gray')
    plt.errorbar(we, Cnle, Cnlerr, marker='*', label=r'$C_{e}$', color='gray', markersize=5)
    #plt.plot(wt, Cnltd, 'H', label='$C_{d}$', color='gray', markersize=3)
    #plt.plot(we, Cnle, '*', label=r'$C_{e}$', color='gray', markersize=8)
    plt.plot(wt, Dt, '-', label='D', color='magenta')
    plt.errorbar(we, De, Derr, marker='o', label=r'$D_{e}$', color='magenta', markersize=5)
    #plt.plot(wt, Dtd, 'X', label='$D_{d}$', color='magenta', markersize=3)
    #plt.plot(we, De, 'o', label=r'$D_{e}$', color='magenta', markersize=8)
    plt.plot(wt, Et, '-.', label='E', color='blue')
    plt.errorbar(we, Ee, Eerr, marker='s', label=r'$E_{e}$', color='blue', markersize=5)
    #plt.plot(wt, Etd, '4', label='$E_{d}$', color='blue', markersize=3)
    #plt.plot(we, Ee, 's', label=r'$E_{e}$', color='blue', markersize=8)
    #plt.errorbar(we, Ee, errE, xerr=None)
    plt.plot(wt, St, ':', label='$S$', color='red')
    plt.errorbar(we, Se, Serr, marker='^', label=r'$S_{e}$', color='red', markersize=5)
    #plt.plot(wt, Std, 'd', label='$S_{d}$', color='red', markersize=3)
    #plt.plot(we, Se, '^', label=r'$S_{e}$', color='red', markersize=8)
    plt.plot(wt, Nlt, '--', label='$N$', color='cyan')
    plt.errorbar(we, Nle, Nlerr, marker='h', label=r'$N_{e}$', color='cyan', markersize=5)
    #plt.plot(wt, Nltd, '+', label='$N_{d}$', color='cyan', markersize=3)
    #plt.plot(we, Nle, 'h', label=r'$N_{e}$', color='cyan', markersize=8)
    plt.xlabel('w')
    plt.legend(loc=6)
    plt.xlim(-0.2,1.02)
    plt.ylim(-0.02,1.02)
    import platform
    if platform.system() == 'Linux':
        plt.savefig('/home/jonas/Dropbox/Research/ibm/bds/calc/qcorr.eps',
                    format='eps', dpi=100)
    else:
        plt.savefig('/Users/jonas/Dropbox/Research/ibm/bds/calc/qcorr.eps',
                    format='eps', dpi=100)
    plt.show()


def werner_decoh():
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
        # Et[j] = ent.concurrence(rho)
        Et[j] = 2*ent.negativity(4, pT.pTransposeL(2, 2, rho))
        Cnlt[j] = coh.coh_nl(2, 2, rho)
        Nlt[j] = ent.chsh(rho)
        St[j] = ent.steering(rho)
        # Dt[j] = discord.hellinger(2, 2, rho)
        Dt[j] = discord.oz_2qb(rho)
        wt[j] = w
        rhod = werner_pdad(w, p, a)
        # Etd[j] = ent.concurrence(rhod)
        Etd[j] = 2*ent.negativity(4, pT.pTransposeL(2, 2, rhod))
        Cnltd[j] = coh.coh_nl(2, 2, rhod)
        Nltd[j] = ent.chsh(rhod)
        Std[j] = ent.steering(rhod)
        # Dtd[j] = discord.hellinger(2, 2, rhod)
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

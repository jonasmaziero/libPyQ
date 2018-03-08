import numpy as np
#-------------------------------------------------------------------------------
def pTransposeTest():
    dl = 2; dm = 2; dr = 2; d = dl*dm*dr
    rho = np.zeros((d,d), dtype = complex)
    rho[1][1] = 1.0/3.0;  rho[1][2] = 1.0/3.0;  rho[1][4] = 1.0/3.0;
    rho[2][1] = 1.0/3.0;  rho[2][2] = 1.0/3.0;  rho[2][4] = 1.0/3.0;
    rho[4][1] = 1.0/3.0;  rho[4][2] = 1.0/3.0;  rho[4][4] = 1.0/3.0;
    rhoTL = np.zeros((d,d), dtype = complex)
    rhoTL = pTransposeL(dl, dm*dr, rho)
    print((np.matrix(rhoTL)).real)
    print('')
#-------------------------------------------------------------------------------
def pTransposeL(dl, dr, rhoLR):
    # Returns the left partial transpose over the 'left' subsystem of rhoLR
    d = dl*dr
    rhoTL = np.zeros((d,d), dtype = complex)
    for jl in range(0,dl):
        for kl in range(0,dl):
            for jr in range(0,dr):
                for kr in range(0,dr):
                    rhoTL[kl*dr+jr][jl*dr+kr] = rhoLR[jl*dr+jr][kl*dr+kr]
    return rhoTL
#-------------------------------------------------------------------------------
def pTransposeR(dl, dr, rhoLR):
    # Returns the left partial transpose over the 'right' subsystem of rhoLR
    d = dl*dr
    rhoTR = np.zeros((d,d), dtype = complex)
    for jl in range(0,dl):
        for kl in range(0,dl):
            for jr in range(0,dr):
                for kr in range(0,dr):
                    rhoTR[jl*dr+kr][kl*dr+jr] = rhoLR[jl*dr+jr][kl*dr+kr]
    return rhoTR
#-------------------------------------------------------------------------------
def pTransposeM(dl, dm, dr, rhoLMR):
    # Returns the PT with relation to system B, for a 3-partite state
    d = dl*dm*dr
    rhoTM = np.zeros((d,d), dtype = complex)
    for jl in range(0,dl):
        for kl in range(0,dl):
            for jm in range(0,dm):
                for km in range(0,dm):
                    for jr in range(0,dr):
                        for kr in range(0,dc):
                            a = jl*dm*dr + km*dr + jr
                            b = kl*dm*dr + jm*dr + kr
                            c = jl*dm*dr + jm*dr + jr
                            d = kl*dm*dr + km*dr + kr
                            rhoTM[a][b] = rhoLMR[c][d]
    return rhoTM
#-------------------------------------------------------------------------------

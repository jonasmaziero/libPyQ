import numpy as np


def pTraceTest():
    dl = 2
    dm = 2
    dr = 2
    d = dl*dm*dr
    rho = np.zeros((d, d), dtype=complex)
    rho[1][1] = 1.0/3.0
    rho[1][2] = 1.0/3.0
    rho[1][4] = 1.0/3.0
    rho[2][1] = 1.0/3.0
    rho[2][2] = 1.0/3.0
    rho[2][4] = 1.0/3.0
    rho[4][1] = 1.0/3.0
    rho[4][2] = 1.0/3.0
    rho[4][4] = 1.0/3.0
    rhoR = np.zeros((dr, dr), dtype=complex)
    rhoR = pTraceL(dl*dm, dr, rho)
    print(np.matrix(rhoR))
    print('from the wrapper')
    import ptrace
    rhoR = ptrace.partial_trace_a(rho, dl*dm, dr)
    print(np.matrix(rhoR))
    print('')
    rhoL = np.zeros((dl, dl), dtype=complex)
    rhoL = pTraceR(dl, dm*dr, rho)
    print(np.matrix(rhoL))
    print('from the wrapper')
    rhoL = ptrace.partial_trace_b(rho, dl, dm*dr)
    print(np.matrix(rhoL))
    print('')
    rhoLR = np.zeros((dl*dr, dl*dr), dtype=complex)
    rhoLR = pTraceM(dl, dm, dr, rho)
    print(np.matrix(rhoLR))
    print('from the wrapper')
    rhoLR = ptrace.partial_trace_3(rho, dl, dm, dr)
    print(np.matrix(rhoLR))


def pTraceL(dl, dr, rhoLR):
    # Returns the left partial trace over the 'left' subsystem of rhoLR
    rhoB = np.zeros((dr, dr), dtype=complex)
    for j in range(0, dr):
        for k in range(j, dr):
            for l in range(0, dl):
                rhoB[j][k] += rhoLR[l*dr+j][l*dr+k]
            if j != k:
                rhoB[k][j] = np.conj(rhoB[j][k])
    return rhoB
#-------------------------------------------------------------------------------


def pTraceR(dl, dr, rhoLR):
    # Returns the right partial trace over the 'right' subsystem of rhoLR
    rhoA = np.zeros((dl, dl), dtype=complex)
    for j in range(0, dl):
        for k in range(j, dl):
            for l in range(0, dr):
                rhoA[j][k] += rhoLR[j*dr+l][k*dr+l]
        if j != k:
            rhoA[k][j] = np.conj(rhoA[j][k])
    return rhoA
#-------------------------------------------------------------------------------


def pTraceM(dl, dm, dr, rhoLMR):
    # Returns the partial trace over the middle subsystem of rhoLMR
    dlr = dl*dr
    rhoLR = np.zeros((dlr, dlr), dtype=complex)
    for j in range(0, dl):
        for l in range(0, dr):
            cj = j*dr + l
            ccj = j*dm*dr + l
            for m in range(0, dl):
                for o in range(0, dr):
                    ck = m*dr + o
                    cck = m*dm*dr + o
                    for k in range(0, dm):
                        rhoLR[cj][ck] += rhoLMR[ccj+k*dr][cck+k*dr]
    return rhoLR
#-------------------------------------------------------------------------------

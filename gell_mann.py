import numpy as np


def rho(d, bv):
    dm = gell_mann(d, 'i', 0, 0)/d
    for j in range(0, d-1):   # diagonal
        dm += (bv[j]/2)*gell_mann(d, 'd', j+1, j+1)
    for k in range(1, d):  # symmetric
        for l in range(k+1, d+1):
            j += 1
            dm += (bv[j]/2)*gell_mann(d, 's', k, l)
    for k in range(1, d):  # antisymmetric
        for l in range(k+1, d+1):
            j += 1
            dm += (bv[j]/2)*gell_mann(d, 'a', k, l)
    return dm


def gell_mann(d, g, j, k):
    gm = np.zeros((d, d), dtype=complex)
    if g == 'i':  # identity
        for m in range(0, d):
            gm[m, m] = 1
    elif g == 'd':  # diagonal
        nt = np.sqrt(2/(j*(j+1)))
        for m in range(0, j):
            gm[m, m] = nt
        gm[j, j] = -j*nt
    elif g == 's':  # symmetric
        gm[j-1, k-1] = 1
        gm[k-1, j-1] = 1
    elif g == 'a':  # anti-symmetric
        gm[j-1, k-1] = -1j
        gm[k-1, j-1] = 1j
    return gm


def bloch_vector(d, A):
    bv = np.zeros(d**2-1)
    for j in range(1, d):   # diagonal
        bv[j-1] = 0
        for k in range(1, j+1):
            bv[j-1] += A[k-1, k-1]
        bv[j-1] -= j*A[j, j]
        bv[j-1] *= np.sqrt(2/(j*(j+1)))
    for k in range(1, d):  # symmetric
        for l in range(k+1, d+1):
            j += 1
            bv[j-1] = A[l-1, k-1] + A[k-1, l-1]
    for k in range(1, d):  # anti-symmetric
        for l in range(k+1, d+1):
            j += 1
            bv[j-1] = -1j*(A[l-1, k-1] - A[k-1, l-1])
    return bv


def corr_mat_dd(da, db, M):
    cmdd = np.zeros((da-1, db-1), dtype=complex)
    for j in range(1, da):
        for k in range(1, db):
            for m in range(1, j+1):
                for n in range(1, k+1):
                    cmdd[j-1, k-1] += M[(m-1)*db+(n-1),(m-1)*db+(n-1)]
            m = j+1
            for n in range(1, k+1):
                cmdd[j-1, k-1] -= j*M[(m-1)*db+(n-1), (m-1)*db+(n-1)]
            n = k+1
            for m in range(1, j+1):
                cmdd[j-1, k-1] -= k*M[(m-1)*db+(n-1), (m-1)*db+(n-1)]
            m = j+1
            n = k+1
            cmdd[j-1, k-1] += j*k*M[(m-1)*db+(n-1), (m-1)*db+(n-1)]
            cmdd[j-1, k-1] *= 2/np.sqrt(j*(j+1)*k*(k+1))
    return cmdd.real


def corr_mat_ds(da, db, M):
    cmds = np.zeros((da-1, db*(db-1)//2), dtype=complex)
    for j in range(1, da):
        n = 0
        for k in range(1, db):
            for l in range(k+1, db+1):
                n += 1
                for m in range(1, j+1):
                    cmds[j-1, n-1] += (M[(m-1)*db+(l-1), (m-1)*db+(k-1)] +
                                       M[(m-1)*db+(k-1), (m-1)*db+(l-1)])
                m = j+1
                cmds[j-1, n-1] -= j*(M[(m-1)*db+(l-1), (m-1)*db+(k-1)] +
                                     M[(m-1)*db+(k-1), (m-1)*db+(l-1)])
                cmds[j-1, n-1] *= np.sqrt(2/(j*(j+1)))
    return cmds.real


def corr_mat_da(da, db, M):
    cmda = np.zeros((da-1, db*(db-1)//2), dtype=complex)
    for j in range(1, da):
        n = 0
        for k in range(1, db):
            for l in range(k+1, db+1):
                n += 1
                for m in range(1, j+1):
                    cmda[j-1, n-1] += (M[(m-1)*db+(l-1), (m-1)*db+(k-1)] -
                                       M[(m-1)*db+(k-1), (m-1)*db+(l-1)])
                m = j+1
                cmda[j-1, n-1] -= j*(M[(m-1)*db+(l-1), (m-1)*db+(k-1)] -
                                     M[(m-1)*db+(k-1), (m-1)*db+(l-1)])
                cmda[j-1, n-1] *= -1j*np.sqrt(2/(j*(j+1)))
    return cmda.real


def corr_mat_sd(da, db, M):
    cmsd = np.zeros((da*(da-1)//2, db-1), dtype=complex)
    n = 0
    for k in range(1, da):
        for l in range(k+1, da+1):
            n += 1
            for j in range(1, db):
                for m in range(1, j+1):
                    cmsd[n-1, j-1] += (M[(l-1)*db+(m-1), (k-1)*db+(m-1)] +
                                       M[(k-1)*db+(m-1), (l-1)*db+(m-1)])
                m = j+1
                cmsd[n-1, j-1] -= j*(M[(l-1)*db+(m-1), (k-1)*db+(m-1)] +
                                     M[(k-1)*db+(m-1), (l-1)*db+(m-1)])
                cmsd[n-1, j-1] *= np.sqrt(2/(j*(j+1)))
    return cmsd.real


def corr_mat_ad(da, db, M):
    cmad = np.zeros((da*(da-1)//2, db-1), dtype=complex)
    n = 0
    for k in range(1, da):
        for l in range(k+1, da+1):
            n += 1
            for j in range(1, db):
                for m in range(1, j+1):
                    cmad[n-1, j-1] += (M[(l-1)*db+m-1, (k-1)*db+m-1] -
                                       M[(k-1)*db+m-1, (l-1)*db+m-1])
                m = j+1
                cmad[n-1, j-1] -= j*(M[(l-1)*db+m-1, (k-1)*db+m-1]
                                     - M[(k-1)*db+m-1, (l-1)*db+m-1])
                cmad[n-1, j-1] *= -1j*np.sqrt(2/(j*(j+1)))
    return cmad.real


def corr_mat_ss(da, db, M):
    cmss = np.zeros((da*(da-1)//2, db*(db-1)//2), dtype=complex)
    p = 0
    for k in range(1, da):
        for l in range(k+1, da+1):
            p += 1
            q = 0
            for m in range(1, db):
                for n in range(m+1, db+1):
                    q += 1
                    cmss[p-1, q-1] += (M[(l-1)*db+n-1, (k-1)*db+m-1] +
                                       M[(k-1)*db+m-1, (l-1)*db+n-1])
                    cmss[p-1, q-1] += (M[(l-1)*db+m-1, (k-1)*db+n-1] +
                                       M[(k-1)*db+n-1, (l-1)*db+m-1])
    return cmss.real


def corr_mat_sa(da, db, M):
    cmsa = np.zeros((da*(da-1)//2, db*(db-1)//2), dtype=complex)
    p = 0
    for k in range(1, da):
        for l in range(k+1, da+1):
            p += 1
            q = 0
            for m in range(1, db):
                for n in range(m+1, db+1):
                    q += 1
                    cmsa[p-1, q-1] -= 1j*(M[(l-1)*db+n-1, (k-1)*db+m-1] -
                                          M[(k-1)*db+m-1, (l-1)*db+n-1])
                    cmsa[p-1, q-1] -= 1j*(M[(k-1)*db+n-1, (l-1)*db+m-1] -
                                          M[(l-1)*db+m-1, (k-1)*db+n-1])
    return cmsa.real


def corr_mat_as(da, db, M):
    cmas = np.zeros((da*(da-1)//2, db*(db-1)//2), dtype=complex)
    p = 0
    for k in range(1, da):
        for l in range(k+1, da+1):
            p += 1
            q = 0
            for m in range(1, db):
                for n in range(m+1, db+1):
                    q += 1
                    cmas[p-1, q-1] -= 1j*(M[(l-1)*db+n-1, (k-1)*db+m-1] -
                                          M[(k-1)*db+m-1, (l-1)*db+n-1])
                    cmas[p-1, q-1] -= 1j*(M[(l-1)*db+m-1, (k-1)*db+n-1] -
                                          M[(k-1)*db+n-1, (l-1)*db+m-1])
    return cmas.real


def corr_mat_aa(da, db, M):
    cmaa = np.zeros((da*(da-1)//2, db*(db-1)//2), dtype=complex)
    p = 0
    for k in range(1, da):
        for l in range(k+1, da+1):
            p += 1
            q = 0
            for m in range(1, db):
                for n in range(m+1, db+1):
                    q += 1
                    cmaa[p-1, q-1] += (M[(l-1)*db+m-1, (k-1)*db+n-1] +
                                       M[(k-1)*db+n-1, (l-1)*db+m-1])
                    cmaa[p-1, q-1] -= (M[(l-1)*db+n-1, (k-1)*db+m-1] +
                                       M[(k-1)*db+m-1, (l-1)*db+n-1])
    return cmaa.real


def corr_mat(da, db, M):
    dda = (da*(da-1))//2
    ddb = (db*(db-1))//2
    cm = np.zeros((da**2-1, db**2-1))
    k = -1
    l = -1
    cmdd = corr_mat_dd(da, db, M)
    cmds = corr_mat_ds(da, db, M)
    cmda = corr_mat_da(da, db, M)
    for m in range(0, da-1):
        k += 1
        for n in range(0, db-1):  # diagonal-diagonal
            l += 1
            cm[k, l] = cmdd[m, n]
        for n in range(0, ddb):  # diagonal-symmetric
            l += 1
            cm[k, l] = cmds[m, n]
        for n in range(0, ddb):  # diagonal-antisymmetric
            l += 1
            cm[k, l] = cmda[m, n]
    cmsd = corr_mat_sd(da, db, M)
    cmss = corr_mat_ss(da, db, M)
    cmsa = corr_mat_sa(da, db, M)
    l = -1
    for m in range(0, dda):
        k += 1
        for n in range(0, db-1):  # diagonal-diagonal
            l += 1
            cm[k, l] = cmsd[m, n]
        for n in range(0, ddb):  # diagonal-symmetric
            l += 1
            cm[k, l] = cmss[m, n]
        for n in range(0, ddb):  # diagonal-antisymmetric
            l += 1
            cm[k, l] = cmsa[m, n]
    cmad = corr_mat_ad(da, db, M)
    cmas = corr_mat_as(da, db, M)
    cmaa = corr_mat_aa(da, db, M)
    l = -1
    for m in range(0, dda):
        k += 1
        for n in range(0, db-1):  # diagonal-diagonal
            l += 1
            cm[k, l] = cmad[m, n]
        for n in range(0, ddb):  # diagonal-symmetric
            l += 1
            cm[k, l] = cmas[m, n]
        for n in range(0, ddb):  # diagonal-antisymmetric
            l += 1
            cm[k, l] = cmaa[m, n]
    return cm

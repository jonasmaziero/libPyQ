
def rho_pd(p, rho):
    # Returns the two-qubit evolved state for local PHASE DAMPING channels
    from numpy import zeros, matmul, kron
    from math import sqrt
    K0 = zeros((2, 2))
    K0[0][0] = sqrt(1.0-p)
    K0[1][1] = K0[0][0]
    K1 = zeros((2, 2))
    K1[0][0] = sqrt(p)
    K2 = zeros((2, 2))
    K2[1][1] = K1[0][0]
    rhop = zeros((4, 4))
    tp = zeros((4, 4))
    tp = kron(K0, K0)
    rhop = rhop + matmul(matmul(tp, rho), tp)
    tp = kron(K0, K1)
    rhop = rhop + matmul(matmul(tp, rho), tp)
    tp = kron(K0, K2)
    rhop = rhop + matmul(matmul(tp, rho), tp)
    tp = kron(K1, K0)
    rhop = rhop + matmul(matmul(tp, rho), tp)
    tp = kron(K1, K1)
    rhop = rhop + matmul(matmul(tp, rho), tp)
    tp = kron(K1, K2)
    rhop = rhop + matmul(matmul(tp, rho), tp)
    tp = kron(K2, K0)
    rhop = rhop + matmul(matmul(tp, rho), tp)
    tp = kron(K2, K1)
    rhop = rhop + matmul(matmul(tp, rho), tp)
    tp = kron(K2, K2)
    rhop = rhop + matmul(matmul(tp, rho), tp)
    return rhop

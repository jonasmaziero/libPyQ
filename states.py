

def werner(w):
    # Returns a two-qubit Werner state (1-w)id/4 + w*psi
    from numpy import array
    rho = array([[0.25*(1.0-w), 0.0, 0.0, 0.0], [0.0, 0.25*(1.0+w), 0.5*w, 0.0],
                 [0.0, 0.5*w, 0.25*(1.0+w), 0.0], [0.0, 0.0, 0.0, 0.25*(1.0-w)]])
    return rho


def rho_bds(c1, c2, c3):
    # Returns a two-qubit Bell-diagonal state (BDS)
    from numpy import zeros, kron
    s0 = zeros((2, 2))
    s1 = zeros((2, 2))
    s2 = zeros((2, 2), dtype=complex)
    s3 = zeros((2, 2))
    s0[0][0] = 1.0
    s0[1][1] = 1.0
    s1[0][1] = 1.0
    s1[1][0] = 1.0
    s2[0][1] = - 1j
    s2[1][0] = 1j
    s3[0][0] = 1.0
    s3[1][1] = -1.0
    rho = zeros((4, 4), dtype=complex)
    rho = (kron(s0, s0) + c1*kron(s1, s1) + c2*kron(s2, s2) + c3*kron(s3, s3))/4.0
    return rho


def psi_1qb(theta, phi):
    # Returns a generic-pure 1-qubit state
    from numpy import zeros
    psi = zeros(2, dtype=complex)
    from math import sin, cos, exp
    psi[0] = cos(theta/2.0)
    psi[1] = (cos(phi) + sin(phi)*1j)*sin(theta/2.0)
    return psi


def rho_1qb(r1, r2, r3):
    # Returns a generic-mixed 1-qubit state
    from numpy import zeros
    s0 = zeros((2, 2))
    s1 = zeros((2, 2))
    s2 = zeros((2, 2), dtype=complex)
    s3 = zeros((2, 2))
    s0[0][0] = 1.0
    s0[1][1] = 1.0
    s1[0][1] = 1.0
    s1[1][0] = 1.0
    s2[0][1] = - 1j
    s2[1][0] = 1j
    s3[0][0] = 1.0
    s3[1][1] = -1.0
    rho = zeros((2, 2), dtype=complex)
    rho = (s0 + r1*s1 + r2*s2 + r3*s3)/2.0
    return rho


def rho_2qb(CM):
    # Returns a generic-mixed 2-qubit state
    from numpy import zeros, kron
    s0 = zeros((2, 2))
    s1 = zeros((2, 2))
    s2 = zeros((2, 2), dtype=complex)
    s3 = zeros((2, 2))
    s0[0][0] = 1.0
    s0[1][1] = 1.0
    s1[0][1] = 1.0
    s1[1][0] = 1.0
    s2[0][1] = - 1j
    s2[1][0] = 1j
    s3[0][0] = 1.0
    s3[1][1] = -1.0
    rho = zeros((4, 4), dtype=complex)
    rho += CM[0][0]*kron(s0, s0) + CM[1][0]*kron(s1, s0)
    rho += CM[2][0]*kron(s2, s0) + CM[3][0]*kron(s3, s0)
    rho += CM[0][1]*kron(s0, s1) + CM[1][1]*kron(s1, s1)
    rho += CM[2][1]*kron(s2, s1) + CM[3][1]*kron(s3, s1)
    rho += CM[0][2]*kron(s0, s2) + CM[1][2]*kron(s1, s2)
    rho += CM[2][2]*kron(s2, s2) + CM[3][2]*kron(s3, s2)
    rho += CM[0][3]*kron(s0, s3) + CM[1][3]*kron(s1, s3)
    rho += CM[2][3]*kron(s2, s3) + CM[3][3]*kron(s3, s3)
    rho = 0.25*rho
    return rho

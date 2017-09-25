#-----------------------------------------------------------------------------------------------------------------------------------
# Returns a two-qubit Bell-diagonal state (BDS)
def rho_bds(c1, c2, c3):
  from numpy import zeros, kron
  s0 = zeros((2,2));  s1 = zeros((2,2));  s2 = zeros((2,2), dtype = complex);  s3 = zeros((2,2))
  s0[0][0] = 1.0; s0[1][1] = 1.0; s1[0][1] = 1.0; s1[1][0] = 1.0; s2[0][1] = - 1j; s2[1][0] = 1j; s3[0][0] = 1.0; s3[1][1] = -1.0
  rho = zeros((4,4), dtype = complex)
  rho = (kron(s0,s0) + c1*kron(s1,s1) + c2*kron(s2,s2) + c3*kron(s3,s3))/4.0
  return rho
#-----------------------------------------------------------------------------------------------------------------------------------
def psi_1qb(theta,phi):
  from numpy import zeros;  psi = zeros(2, dtype = complex)
  from math import sin, cos, exp
  psi[0] = cos(theta/2.0);  psi[1] = (cos(phi) + sin(phi)*1j)*sin(theta/2.0)
  return psi
#-----------------------------------------------------------------------------------------------------------------------------------
def rho_1qb(r1, r2, r3):
  from numpy import zeros
  s0 = zeros((2,2));  s1 = zeros((2,2));  s2 = zeros((2,2), dtype = complex);  s3 = zeros((2,2))
  s0[0][0] = 1.0; s0[1][1] = 1.0; s1[0][1] = 1.0; s1[1][0] = 1.0; s2[0][1] = - 1j; s2[1][0] = 1j; s3[0][0] = 1.0; s3[1][1] = -1.0
  rho = zeros((4,4), dtype = complex);  rho = (s0 + r1*s1 + r2*s2 + r3*s3)/2.0
  return rho
#-----------------------------------------------------------------------------------------------------------------------------------


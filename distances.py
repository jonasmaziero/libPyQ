#------------------------------------------------------------------------------------------------------------------------------------
def test_distances():
  from numpy import zeros
  from math import sqrt
  
  #fidelity pure - pure
  z = zeros(2, dtype = complex);  w = zeros(2, dtype = complex)
  z[0] = 1.0;  z[1] = 0.0
  w[0] = 1.0/sqrt(2.0);  w[1] = 1.0/sqrt(2.0)
  print(fidelity_pp(z, w))
  '''
  #fidelity pure - mixed
  x = zeros(2, dtype = complex);  y = zeros(4, dtype = complex)
  x = [1.0/sqrt(2.0),1.0/sqrt(2.0)*(1j)]
  y = [[1,0],[0,0]]
  print(fidelity_pm(x, y))
  '''
  #fidelity mixed - mixed
  x = zeros((2,2), dtype = complex);  y = zeros((2,2), dtype = complex)
  x = [[1.0,0.0],[0.0,0.0]]
  y = [[1.0/2.0,1.0/2.0],[1.0/2.0,1.0/2.0]]
  print(fidelity_mm(2, x, y))
  
#------------------------------------------------------------------------------------------------------------------------------------
# Returns the FIDELITY between 2 PURE states
def fidelity_pp(psi, phi):
  from numpy import conjugate, inner
  psiD = conjugate(psi)
  F = abs(inner(psiD, phi))
  return F
#------------------------------------------------------------------------------------------------------------------------------------
# Returns the FIDELITY between a PURE and a MIXED state
def fidelity_pm(psi, rho):
  from numpy import conjugate, inner, matmul
  psiD = conjugate(psi);  phi = matmul(rho, psi);  F = abs(inner(psiD, phi))
  return F
#------------------------------------------------------------------------------------------------------------------------------------
# Returns the fidelity between 2 MIXED states
def fidelity_mm(d, rho, zeta):  
  #import scipy.linalg.lapack as lapak;  eigr = lapak.zheevd(rho)
  from numpy import linalg as LA;  val, evec = LA.eigh(rho)
  from numpy import kron, zeros;  from math import sqrt;  import matrix_functions as mf
  A = zeros((d,d), dtype = complex);  phi = zeros(d, dtype = complex);  psi = zeros(d, dtype = complex)
  for j in range(0,d):
    for l in range(0,d):
      #phi[l] = eigr[1][l][j]
      phi[l] = evec[l][j]
    for k in range(0,d):
      for m in range(0,d):
        #psi[m] = eigr[1][m][k]
        psi[m] = evec[m][k]
      #A = A + (sqrt(eigr[0][j]*eigr[0][k])*mf.sandwich(d, phi, zeta, psi))*mf.outer(d, phi, psi)
      A = A + (sqrt(val[j]*val[k])*mf.sandwich(d, phi, zeta, psi))*mf.outer(d, phi, psi)
  #eigA = lapak.zheevd(A)
  eigA = LA.eigvalsh(A);  F = 0.0
  for n in range(0,d):
    F = F + sqrt(eigA[n])
  return F
#-----------------------------------------------------------------------------------------------------------------------------------
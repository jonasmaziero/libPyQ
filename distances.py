#------------------------------------------------------------------------------------------------------------------------------------
def test_distances():
  from numpy import zeros
  from math import sqrt
  '''#fidelity pure - pure
  x = np.zeros(3, dtype = complex);  y = np.zeros(3, dtype = complex)
  x[0] = 1.0/sqrt(2.0);  x[1] = (1.0/sqrt(2.0))*(1j)
  y[0] = 1.0/sqrt(2.0);  y[1] = 1.0/sqrt(2.0)
  print(fidelity_pp(x, y))'''
  
  #fidelity pure - mixed
  x = zeros(2, dtype = complex);  y = zeros(4, dtype = complex)
  x = [1.0/sqrt(2.0),1.0/sqrt(2.0)*(1j)]
  y = [[1,0],[0,0]]
  print(fidelity_pm(x, y))
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
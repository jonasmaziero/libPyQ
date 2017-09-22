#------------------------------------------------------------------------------------------------------------------------------------
import numpy as np
from math import sqrt
#------------------------------------------------------------------------------------------------------------------------------------
def test_distances():
  '''#fidelity pure - pure
  x = np.zeros(3, dtype = complex);  y = np.zeros(3, dtype = complex)
  x[0] = 1.0/sqrt(2.0);  x[1] = (1.0/sqrt(2.0))*(1j)
  y[0] = 1.0/sqrt(2.0);  y[1] = 1.0/sqrt(2.0)
  print(fidelity_pp(x, y))'''
  
  #fidelity pure - mixed
  x = np.zeros(2, dtype = complex);  y = np.zeros(4, dtype = complex)
  x = [1.0/sqrt(2.0),1.0/sqrt(2.0)*(1j)]
  y = [[1,0],[0,0]]
  print(fidelity_pm(x, y))
#------------------------------------------------------------------------------------------------------------------------------------
# Returns the FIDELITY between 2 PURE states
def fidelity_pp(psi, phi):
  psiD = np.conjugate(psi)
  F = abs(np.inner(psiD, phi))
  return F
#------------------------------------------------------------------------------------------------------------------------------------
# Returns the FIDELITY between a PURE and a MIXED state
def fidelity_pm(psi, rho):
  psiD = np.conjugate(psi)
  phi = np.matmul(rho, psi)
  F = abs(np.inner(psiD, phi))
  return F
#------------------------------------------------------------------------------------------------------------------------------------
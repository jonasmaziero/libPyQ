#------------------------------------------------------------------------------------------------------------------------------------
# Returns the purity of a dxd density matrix rho
def purity(d,rho):
  purity = 0.0
  j = -1
  while (j < d-1):
    j = j + 1
    k = -1
    while (k < d-1): 
      k = k + 1
      purity = purity + (rho[j][k].real)**2.0 + (rho[j][k].imag)**2.0
  return purity
#------------------------------------------------------------------------------------------------------------------------------------
# Returns the Shannon entropy of a probability vector pv, of dimension d
def shannon(d,pv):
  from math import log2
  shannon = 0.0
  j = -1
  while (j < d-1):
    j = j + 1
    if pv[j] > 1.e-15 and pv[j] < (1.0-1.e-15):
      shannon = shannon - pv[j]*log2(pv[j])
  return shannon
#------------------------------------------------------------------------------------------------------------------------------------
# Returns the von Neumann entropy of a dxd density matrix rho
def neumann(d,rho):
  import scipy.linalg.lapack as lapak
  b = lapak.zheevd(rho)
  neumann = shannon(d,b[0])
  return neumann
#------------------------------------------------------------------------------------------------------------------------------------

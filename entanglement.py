#-----------------------------------------------------------------------------------------------------------------------------------
def test_entanglement():
  from pTranspose import Ta
  from numpy import zeros;  Ec = zeros(100);  EF = zeros(100);  En = zeros(100);  Eln = zeros(100);  x = zeros(100)
  from states import werner;  dw = 1.0/100.0;  w = -dw
  for j in range(0,100):
    w = w + dw
    if w > 1.0:
      break
    rho = werner(w);  Ec[j] = concurrence(rho);  EF[j] = EoF(rho)
    rhoTp = Ta(2, 2, rho);  En[j] = negativity(4, rhoTp);  Eln[j] = log_negativity(4, rhoTp)
    x[j] = w
  import matplotlib.pyplot as plt
  plt.plot(x, Ec, label = 'Ec');  plt.plot(x, EF, label = 'EoF');  plt.plot(x, En, label = 'En');  plt.plot(x, Eln, label = 'Eln')
  plt.xlabel('x');  plt.ylabel('');  plt.legend(loc=4);  plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------
# Returns the entanglement CONCURRENCE, for two-qubit states
def concurrence(rho):
  from numpy import zeros, conjugate, dot
  s2s2 = zeros((4,4), dtype = complex);  rhoc = zeros((4,4), dtype = complex);  rhoa = zeros((4,4), dtype = complex)
  rhot = zeros((4,4), dtype = complex);  R = zeros((4,4), dtype = complex);  ev = zeros(4)
  s2s2 = [[0.0,0.0,0.0,-1.0],[0.0,0.0,1.0,0.0],[0.0,1.0,0.0,0.0],[-1.0,0.0,0.0,0.0]] 
  rhoc = conjugate(rho);  rhoa = dot(s2s2,rhoc);  rhot = dot(rhoa,s2s2);  R = dot(rho,rhot)
  from numpy import linalg as LA;  ev = LA.eigvalsh(R);  evm = max(ev[0], ev[1], ev[2], ev[3])
  from math import sqrt
  C = 2.0*sqrt(abs(evm)) - sqrt(abs(ev[0])) - sqrt(abs(ev[1])) - sqrt(abs(ev[2])) - sqrt(abs(ev[3]))
  if C < 0.0:
    C = 0.0
  return C
#-----------------------------------------------------------------------------------------------------------------------------------
# Returns the ENTANGLEMENT of FORMATION for two-qubit states
def EoF(rho):
  from numpy import zeros;  pv = zeros(2)
  from math import sqrt;  Ec = concurrence(rho);  pv[0] = (1.0 + sqrt(1.0 - Ec**2.0))/2.0;  pv[1] = 1.0 - pv[0]
  from entropy import shannon;  EF = shannon(2, pv)
  return EF
#-----------------------------------------------------------------------------------------------------------------------------------
# Returns the entanglement NEGATIVITY of a "bipartite" system
# En is equal to the sum of the absolute values of the negative eigenvalues of rho_pt
def negativity(d, rhoTp):
  from distances import normTr
  En = 0.5*(normTr(d, rhoTp) - 1.0)  
  return En
#-----------------------------------------------------------------------------------------------------------------------------------
# Returns the entanglement LOGARITMIC NEGATIVITY of a "bipartite" system
def log_negativity(d, rhoTp):
  from math import log;  En = negativity(d, rhoTp);  Eln = log(2.0*En+1.0, 2)
  return Eln
#-----------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------
def test():
  from distances import fidelity_pp;  from numpy import zeros;  from math import sqrt
  ns = 10**3 # number of samples for the average
  nqb = 5 # maximum number of qubits regarded
  Favg = zeros(nqb);  Fexa = zeros(nqb);  d = zeros(nqb, dtype = int)
  for j in range(0,nqb):
    d[j] = 2**(j+1);  psi = zeros(d[j], dtype = complex);  phi = zeros(d[j], dtype = complex)
    Fexa[j] = 1.0/d[j]
    Favg[j] = 0.0
    for k in range(0,ns):
      psi = rsvg(d[j]);  phi = rsvg(d[j]);  Favg[j] = Favg[j] + fidelity_pp(psi, phi)
    Favg[j] = Favg[j]/ns
  import matplotlib.pyplot as plt
  plt.plot(d, Favg, label = '<F>');  plt.plot(d, Fexa, label = 'F')
  plt.xlabel('d');  plt.ylabel('F');  plt.legend()
  plt.show()
#------------------------------------------------------------------------------------------------------------------------------------
def rsvg(d):
  from numpy import zeros;  rn = zeros(d);  rpv = zeros(d);  rsv = zeros(d, dtype = complex)
  from rpvg import rpv_zhsl;  rpv = rpv_zhsl(d)
  from math import pi, sqrt, cos, sin;  tpi = 2.0*pi
  from random import random
  for j in range(0,d):
    rn[j] = random();  arg = tpi*rn[j];  ph = cos(arg) + (1j)*sin(arg)
    rsv[j] = sqrt(rpv[j])*ph
  return rsv
#------------------------------------------------------------------------------------------------------------------------------------
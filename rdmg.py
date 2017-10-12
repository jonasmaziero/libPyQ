#------------------------------------------------------------------------------------------------------------------------------------
def test():
  from numpy import random, zeros;  random.seed()
  from coherence import coh_l1n, coh_re
  ns = 10**3 # number of samples for the average
  nqb = 5 # maximum number of qubits regarded
  Cavg = zeros(nqb);  d = zeros(nqb, dtype = int)
  for j in range(0,nqb):
    d[j] = 2**(j+1);  rdm = zeros((d[j],d[j]), dtype = complex)
    Cavg[j] = 0.0
    for k in range(0,ns):
      rdm = rdm_ginibre(d[j]);  Cavg[j] = Cavg[j] + coh_re(d[j], rdm)
    Cavg[j] = Cavg[j]/ns
  import matplotlib.pyplot as plt
  plt.plot(d, Cavg, label = '')
  plt.xlabel('d');  plt.ylabel('C');  plt.legend()
  plt.show()
#------------------------------------------------------------------------------------------------------------------------------------
def rdm_ginibre(d):
  from numpy import zeros;  rdm = zeros((d,d), dtype = complex);  G = ginibre(d)
  from distances import normHS;  N2 = (normHS(d, G))**2.0
  for j in range(0,d):
    for k in range(0,d):
      rdm[j][k] = 0.0
      for l in range(0,d):
        rdm[j][k] = rdm[j][k] + (G[j][l].real)*(G[k][l].real) + (G[j][l].imag)*(G[k][l].imag) \
                              - (1j)*((G[j][l].real)*(G[k][l].imag) - (G[j][l].imag)*(G[k][l].real))
      rdm[j][k] = rdm[j][k]/N2
  return rdm
#------------------------------------------------------------------------------------------------------------------------------------
def ginibre(d):
  from numpy import random, zeros;  G = zeros((d,d), dtype = complex);  mu, sigma = 0.0, 1.0
  for j in range(0,d):
    grn = random.normal(mu, sigma, 2*d)
    for k in range(0,d):
      G[j][k] = grn[k] + (1j)*grn[k+d]
  return G
#------------------------------------------------------------------------------------------------------------------------------------
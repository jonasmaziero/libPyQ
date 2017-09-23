#------------------------------------------------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
from math import log
import kraus
#import plots
#------------------------------------------------------------------------------------------------------------------------------------
def test_discord():
  vp = np.zeros(100);  
  cc1 = np.zeros(100);  mi1 = np.zeros(100);  dd1 = np.zeros(100)
  cc2 = np.zeros(100);  mi2 = np.zeros(100);  dd2 = np.zeros(100)
  cc3 = np.zeros(100);  mi3 = np.zeros(100);  dd3 = np.zeros(100)
  rho1 = np.zeros((4,4));  rho2 = np.zeros((4,4));  rho3 = np.zeros((4,4));  rhop = np.zeros((4,4))
  import states;  
  c1 = 0.5;  c2 = -c1;  c3 = 0.25;  rho1 = states.rho_bds(c1, c2, c3)
  c1 = 0.3;  c2 = -c1;  c3 = 0.25;  rho2 = states.rho_bds(c1, c2, c3)
  c1 = 0.1;  c2 = -c1;  c3 = 0.25;  rho3 = states.rho_bds(c1, c2, c3)
  dp = 0.01; p = -dp
  for j in range(0,100):
    p = p + dp;  vp[j] = p
    rhop = kraus.rho_pd(p, rho1)
    cc1[j] = ccorr_hv_bds(rhop);  dd1[j] = discord_oz_bds(rhop);  mi1[j] = mi_bds(rhop)
    rhop = kraus.rho_pd(p, rho2)
    cc2[j] = ccorr_hv_bds(rhop);  dd2[j] = discord_oz_bds(rhop);  mi2[j] = mi_bds(rhop)
    rhop = kraus.rho_pd(p, rho3)
    cc3[j] = ccorr_hv_bds(rhop);  dd3[j] = discord_oz_bds(rhop);  mi3[j] = mi_bds(rhop)
  plt.plot(vp,cc1,label='C1')
  plt.plot(vp,cc2,label='C2')
  plt.plot(vp,cc3,label='C3')
  plt.xlabel('p')
  plt.legend()
  plt.show()
#------------------------------------------------------------------------------------------------------------------------------------
# Returns the OLLIVIER-ZUREK discord for 2-qubit Bell-diagonal states
def discord_oz_bds(rho): 
  D = mi_bds(rho) - ccorr_hv_bds(rho)
  return D
#---------------------------------
# Returns the Henderson-Vedral classical correlation for 2-qubit Bell-diagonal states
def ccorr_hv_bds(rho):  
  c3 = 4.0*rho[0][0] - 1.0;  c1 = 2.0*(rho[0][3] + rho[1][2]);  c2 = 2.0*(rho[1][2] - rho[0][3])
  c = max(abs(c1), abs(c2), abs(c3))
  cc = 0.5*((1.0-c)*log(1.0-c,2) + (1.0+c)*log(1.0+c,2))
  return cc
#---------------------------------
# Returns the mutual information for 2-qubit Bell-diagonal states
def mi_bds(rho):
  c3 = 4.0*rho[0][0] - 1.0;  c1 = 2.0*(rho[0][3] + rho[1][2]);  c2 = 2.0*(rho[1][2] - rho[0][3])
  l00 = (1.0 + c1 - c2 + c3)/4.0;  l01 = (1.0 + c1 + c2 - c3)/4.0;  l10 = (1.0 - c1 + c2 + c3)/4.0;  l11 = (1.0 - c1 - c2 - c3)/4.0
  mi = l00*log(4.0*l00,2) + l01*log(4.0*l01,2) + l10*log(4.0*l10,2) + l11*log(4.0*l11,2)
  return mi
#------------------------------------------------------------------------------------------------------------------------------------

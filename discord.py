#------------------------------------------------------------------------------------------------------------------------------------
def test_discord():
  from numpy import zeros
  
  vp = zeros(100)
  CC1 = zeros(100);  MI1 = zeros(100);  DD1 = zeros(100)
  CC2 = zeros(100);  MI2 = zeros(100);  DD2 = zeros(100)
  CC3 = zeros(100);  MI3 = zeros(100);  DD3 = zeros(100)
  CC4 = zeros(100);  MI4 = zeros(100);  DD4 = zeros(100)
  rho1 = zeros((4,4));  rho2 = zeros((4,4));  rho3 = zeros((4,4));  rho4 = zeros((4,4));  rhop = zeros((4,4))
  import states
  import scipy.linalg.lapack as lapak
  c1 = 0.60;  c2 = -c1;  c3 = 0.25;  rho1 = states.rho_bds(c1, c2, c3);  b = lapak.zheevd(rho1);  print(b[0])
  c1 = 0.40;  c2 = -c1;  c3 = 0.25;  rho2 = states.rho_bds(c1, c2, c3);  b = lapak.zheevd(rho2);  print(b[0])
  c1 = 0.25;  c2 = -c1;  c3 = 0.25;  rho3 = states.rho_bds(c1, c2, c3);  b = lapak.zheevd(rho3);  print(b[0])
  c1 = 0.05;  c2 = -c1;  c3 = 0.25;  rho4 = states.rho_bds(c1, c2, c3);  b = lapak.zheevd(rho3);  print(b[0])
  import kraus
  dp = 0.01; p = -dp
  for j in range(0,100):
    p = p + dp;  vp[j] = p
    rhop = kraus.rho_pd(p, rho1);  CC1[j] = ccorr_hv_bds(rhop);  DD1[j] = discord_oz_bds(rhop);  MI1[j] = mi_bds(rhop)
    rhop = kraus.rho_pd(p, rho2);  CC2[j] = ccorr_hv_bds(rhop);  DD2[j] = discord_oz_bds(rhop);  MI2[j] = mi_bds(rhop)
    rhop = kraus.rho_pd(p, rho3);  CC3[j] = ccorr_hv_bds(rhop);  DD3[j] = discord_oz_bds(rhop);  MI3[j] = mi_bds(rhop)
    rhop = kraus.rho_pd(p, rho4);  CC4[j] = ccorr_hv_bds(rhop);  DD4[j] = discord_oz_bds(rhop);  MI4[j] = mi_bds(rhop)
  import matplotlib.pyplot as plt
  plt.plot(vp,CC1,label='C1(c1=0.6)')#;  plt.plot(vp,DD1,label='D1');  plt.plot(vp,MI1,label='I1')
  plt.plot(vp,CC2,label='C2(c1=0.4)')#;  plt.plot(vp,DD2,label='D2');  plt.plot(vp,MI2,label='I2')
  plt.plot(vp,CC3,label='C3(c1=0.25)')#;  plt.plot(vp,DD3,label='D3');  plt.plot(vp,MI3,label='I3')
  plt.plot(vp,CC4,label='C4(c1=0.05)')#;  plt.plot(vp,DD4,label='D4');  plt.plot(vp,MI4,label='I4')
  plt.xlabel('p')
  plt.legend()
  plt.show()
  '''
  from math import cos
  F1 = zeros(100);  F2 = zeros(100);  F3 = zeros(100);  F4 = zeros(100);  vtheta = zeros(100)
  dtheta = 0.01; theta = -dtheta
  for j in range(0,100):
    theta = theta + dtheta;  vtheta[j] = theta;  alpha = cos(theta/2.0)
    c1 = 0.60;  c2 = -c1;  c3 = 0.25;  F1[j] = 0.5*(1.0 + c3 + 4.0*(c1-c3)*(alpha**2.0)*(1.0-alpha**2.0))
    c1 = 0.40;  c2 = -c1;  c3 = 0.25;  F2[j] = 0.5*(1.0 + c3 + 4.0*(c1-c3)*(alpha**2.0)*(1.0-alpha**2.0))
    c1 = 0.25;  c2 = -c1;  c3 = 0.25;  F3[j] = 0.5*(1.0 + c3 + 4.0*(c1-c3)*(alpha**2.0)*(1.0-alpha**2.0))
    c1 = 0.05;  c2 = -c1;  c3 = 0.25;  F4[j] = 0.5*(1.0 + c3 + 4.0*(c1-c3)*(alpha**2.0)*(1.0-alpha**2.0))
  import matplotlib.pyplot as plt
  plt.plot(vtheta,F1,label='F1(c1=0.6)')
  plt.plot(vtheta,F2,label='F2(c1=0.4)')
  plt.plot(vtheta,F3,label='F3(c1=0.25)')
  plt.plot(vtheta,F4,label='F4(c1=0.05)')
  plt.xlabel('theta')
  plt.legend()
  plt.show()
  '''
#------------------------------------------------------------------------------------------------------------------------------------
# Returns the OLLIVIER-ZUREK discord for 2-qubit Bell-diagonal states
def discord_oz_bds(rho):
  D = mi_bds(rho) - ccorr_hv_bds(rho)
  return D
#---------------------------------
# Returns the Henderson-Vedral classical correlation for 2-qubit Bell-diagonal states
def ccorr_hv_bds(rho):
  from math import log
  c3 = 4.0*rho[0][0] - 1.0;  c1 = 2.0*(rho[0][3] + rho[1][2]);  c2 = 2.0*(rho[1][2] - rho[0][3])
  c = max(abs(c1), abs(c2), abs(c3))
  cc = 0.5*((1.0-c)*log(1.0-c,2) + (1.0+c)*log(1.0+c,2))
  return cc
#---------------------------------
# Returns the mutual information for 2-qubit Bell-diagonal states
def mi_bds(rho):
  from math import log
  c3 = 4.0*rho[0][0] - 1.0;  c1 = 2.0*(rho[0][3] + rho[1][2]);  c2 = 2.0*(rho[1][2] - rho[0][3])
  l00 = (1.0 + c1 - c2 + c3)/4.0;  l01 = (1.0 + c1 + c2 - c3)/4.0;  l10 = (1.0 - c1 + c2 + c3)/4.0;  l11 = (1.0 - c1 - c2 - c3)/4.0
  mi = l00*log(4.0*l00,2) + l01*log(4.0*l01,2) + l10*log(4.0*l10,2) + l11*log(4.0*l11,2)
  return mi
#------------------------------------------------------------------------------------------------------------------------------------



def discord_oz_bds(rho):
    # Returns the OLLIVIER-ZUREK discord for 2-qubit Bell-diagonal states
    D = mi_bds(rho) - ccorr_hv_bds(rho)
    return D


def ccorr_hv_bds(rho):
    # Returns the Henderson-Vedral classical correlation for 2-qubit
    # Bell-diagonal states
    from math import log
    c3 = 4.0*rho[0][0] - 1.0
    c1 = 2.0*(rho[0][3] + rho[1][2])
    c2 = 2.0*(rho[1][2] - rho[0][3])
    c = max(abs(c1), abs(c2), abs(c3))
    cc = float(0.5*((1.0-c)*log(1.0-c, 2) + (1.0+c)*log(1.0+c, 2)))
    return cc


def mi_bds(rho):
    # Returns the mutual information for 2-qubit Bell-diagonal states
    from math import log
    c3 = 4.0*rho[0][0] - 1.0
    c1 = 2.0*(rho[0][3] + rho[1][2])
    c2 = 2.0*(rho[1][2] - rho[0][3])
    l00 = (1.0 + c1 - c2 + c3)/4.0
    l01 = (1.0 + c1 + c2 - c3)/4.0
    l10 = (1.0 - c1 + c2 + c3)/4.0
    l11 = (1.0 - c1 - c2 - c3)/4.0
    mi = float(l00*log(4.0*l00, 2) + l01*log(4.0*l01, 2) +
               l10*log(4.0*l10, 2) + l11*log(4.0*l11, 2))
    return mi


def test_discord():
    from numpy import zeros
    '''
  vp = zeros(100)
  CC1 = zeros(100);  MI1 = zeros(100);  DD1 = zeros(100)
  CC2 = zeros(100);  MI2 = zeros(100);  DD2 = zeros(100)
  CC3 = zeros(100);  MI3 = zeros(100);  DD3 = zeros(100)
  CC4 = zeros(100);  MI4 = zeros(100);  DD4 = zeros(100)
  rho1 = zeros((4,4));  rho2 = zeros((4,4));  rho3 = zeros((4,4))
  rho4 = zeros((4,4));  rhop = zeros((4,4))
  import states
  import scipy.linalg.lapack as lapak
  c1 = 0.60;  c2 = -c1;  c3 = 0.25;  rho1 = states.rho_bds(c1, c2, c3)
  b = lapak.zheevd(rho1);  print(b[0])
  c1 = 0.40;  c2 = -c1;  c3 = 0.25;  rho2 = states.rho_bds(c1, c2, c3)
  b = lapak.zheevd(rho2);  print(b[0])
  c1 = 0.25;  c2 = -c1;  c3 = 0.25;  rho3 = states.rho_bds(c1, c2, c3)
  b = lapak.zheevd(rho3);  print(b[0])
  c1 = 0.05;  c2 = -c1;  c3 = 0.25;  rho4 = states.rho_bds(c1, c2, c3)
  b = lapak.zheevd(rho3);  print(b[0])
  import kraus
  dp = 1.0/100.0; p = -dp
  for j in range(0,100):
    p = p + dp;  vp[j] = p
    if p > 1.0:
      break
    rhop = kraus.rho_pd(p, rho1);  CC1[j] = ccorr_hv_bds(rhop)
    DD1[j] = discord_oz_bds(rhop);  MI1[j] = mi_bds(rhop)
    rhop = kraus.rho_pd(p, rho2);  CC2[j] = ccorr_hv_bds(rhop)
    DD2[j] = discord_oz_bds(rhop);  MI2[j] = mi_bds(rhop)
    rhop = kraus.rho_pd(p, rho3);  CC3[j] = ccorr_hv_bds(rhop)
    DD3[j] = discord_oz_bds(rhop);  MI3[j] = mi_bds(rhop)
    rhop = kraus.rho_pd(p, rho4);  CC4[j] = ccorr_hv_bds(rhop)
    DD4[j] = discord_oz_bds(rhop);  MI4[j] = mi_bds(rhop)
  import matplotlib.pyplot as plt
  plt.plot(vp,CC1,label='C1(c1=0.6)')#;  plt.plot(vp,DD1,label='D1')
  plt.plot(vp,MI1,label='I1')
  plt.plot(vp,CC2,label='C2(c1=0.4)')#;  plt.plot(vp,DD2,label='D2')
  plt.plot(vp,MI2,label='I2')
  plt.plot(vp,CC3,label='C3(c1=0.25)')#;  plt.plot(vp,DD3,label='D3')
  plt.plot(vp,MI3,label='I3')
  plt.plot(vp,CC4,label='C4(c1=0.05)')#;  plt.plot(vp,DD4,label='D4')
  plt.plot(vp,MI4,label='I4')
  plt.xlabel('p')
  plt.legend()
  plt.show()
  '''
    from math import cos
    F1 = zeros(100)
    F2 = zeros(100)
    F3 = zeros(100)
    F4 = zeros(100)
    vtheta = zeros(100)
    from math import pi
    dtheta = (pi/2.0)/100.0
    theta = -dtheta
    for j in range(0, 100):
        theta = theta + dtheta
        vtheta[j] = theta
        alpha = cos(theta/2.0)
        c1 = 0.60
        c2 = -c1
        c3 = 0.25
        F1[j] = 0.5*(1.0 + c3 + 4.0*(c1-c3)*(alpha**2.0)*(1.0-alpha**2.0))
        c1 = 0.40
        c2 = -c1
        c3 = 0.25
        F2[j] = 0.5*(1.0 + c3 + 4.0*(c1-c3)*(alpha**2.0)*(1.0-alpha**2.0))
        c1 = 0.25
        c2 = -c1
        c3 = 0.25
        F3[j] = 0.5*(1.0 + c3 + 4.0*(c1-c3)*(alpha**2.0)*(1.0-alpha**2.0))
        c1 = 0.05
        c2 = -c1
        c3 = 0.25
        F4[j] = 0.5*(1.0 + c3 + 4.0*(c1-c3)*(alpha**2.0)*(1.0-alpha**2.0))
    from distances import fidelity_pm
    from states import psi_1qb, rho_1qb
    FE1 = zeros(100)
    vt = zeros(100)
    vt[0] = 0.0
    psi = psi_1qb(0.0, 0.0)
    rho = rho_1qb(0.158, 0.16, 0.37)
    FE1[0] = fidelity_pm(psi, rho)
    print(vt[0], FE1[0])
    vt[1] = 0.2
    psi = psi_1qb(0.2, 0.0)
    rho = rho_1qb(0.202, 0.150, 0.332)
    FE1[1] = fidelity_pm(psi, rho)
    print(vt[1], FE1[1])
    vt[2] = 0.4
    psi = psi_1qb(0.4, 0.0)
    rho = rho_1qb(0.226, 0.172, 0.294)
    FE1[2] = fidelity_pm(psi, rho)
    print(vt[2], FE1[2])
    vt[3] = 0.6
    psi = psi_1qb(0.6, 0.0)
    rho = rho_1qb(0.226, 0.128, 0.266)
    FE1[3] = fidelity_pm(psi, rho)
    print(vt[3], FE1[3])
    vt[4] = 0.8
    psi = psi_1qb(0.8, 0.0)
    rho = rho_1qb(0.246, 0.144, 0.238)
    FE1[4] = fidelity_pm(psi, rho)
    print(vt[4], FE1[4])
    vt[5] = 1.0
    psi = psi_1qb(1.0, 0.0)
    rho = rho_1qb(0.212, 0.118, 0.178)
    FE1[5] = fidelity_pm(psi, rho)
    print(vt[5], FE1[5])
    vt[6] = 1.2
    psi = psi_1qb(1.2, 0.0)
    rho = rho_1qb(0.212, 0.118, 0.178)
    FE1[6] = fidelity_pm(psi, rho)
    print(vt[6], FE1[6])
    vt[7] = 1.4
    psi = psi_1qb(1.4, 0.0)
    rho = rho_1qb(0.212, 0.118, 0.178)
    FE1[7] = fidelity_pm(psi, rho)
    print(vt[7], FE1[7])
    vt[8] = pi/2.0
    psi = psi_1qb(pi/2.0, 0.0)
    rho = rho_1qb(0.212, 0.118, 0.178)
    FE1[8] = fidelity_pm(psi, rho)
    print(vt[8], FE1[8])

    import matplotlib.pyplot as plt
    plt.plot(vtheta, F1, label='F1(c1=0.6)')
    plt.plot(vt, FE1, label='FE(c1=0.6)')
    plt.plot(vtheta, F2, label='F2(c1=0.4)')
    plt.plot(vtheta, F3, label='F3(c1=0.25)')
    plt.plot(vtheta, F4, label='F4(c1=0.05)')
    plt.xlabel('theta')
    plt.legend()
    plt.show()

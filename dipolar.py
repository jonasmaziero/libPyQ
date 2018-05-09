import numpy as np
import matplotlib.pyplot as plt
from matFunc import proj, Dagger
from states import Bell, rho1qb
from entanglement import concurrence
import conc


def Epsi(tha, thb, D, t):
    aa = np.cos(tha/2)
    ba = np.sin(tha/2)
    ab = np.cos(thb/2)
    bb = np.sin(thb/2)
    f = 2*aa*ba*ab*bb*(np.cos(2*D*t) - np.cos(4*D*t))
    g = (aa**2*bb**2 + ba**2*ab**2)*np.sin(2*D*t) + 2*aa*ba*ab*bb*np.sin(4*D*t)
    return np.sqrt(f**2 + g**2)


# def Udd(D, t):
#    return proj(Bell(1, 1)) + np.exp(-1j*2*D*t)*proj(Bell(0, 1)) + np.exp(1j*D*t)*(proj(Bell(1, 0)) + proj(Bell(0, 0)))


# def rho3t(r3a, r3b, D, t):
#    return Udd(D, t)*np.kron(rho1qb(0, 0, r3a), rho1qb(0, 0, r3b))*Dagger(Udd(D, t))


# def rho1t(r1a, r1b, D, t):
#    return Udd(D, t)*np.kron(rho1qb(r1a, 0, 0), rho1qb(r1b, 0, 0))*Dagger(Udd(D, t))

#import concurrence


def Erho(a, b, D, t):
    rho3 = 0.25*np.array([[(1+a)*(1+b), 0, 0, 0],
                          [0, 1-a*b+(a-b)*np.cos(2*D*t), 1j*(a-b)*np.sin(2*D*t), 0],
                          [0, -1j*(a-b)*np.sin(2*D*t), 1-a*b-(a-b)*np.cos(2*D*t), 0],
                          [0, 0, 0, (1-a)*(1-b)]])
    cc = conc.concurrence_2qb(rho3)
    return cc


xll = -1
xul = 1
yll = 0
yul = np.pi
x = np.linspace(xll, xul, 80)
y = np.linspace(yll, yul, 80)
X, Y = np.meshgrid(x, y)
D = 1
r3b = 1
#Z = Epsi(X, 1, D, Y)
Z = Erho(X, r3b, D, Y)
contours = plt.contour(X, Y, Z, 6, colors='black')
plt.clabel(contours, inline=True, fontsize=7)
plt.imshow(Z, extent=[xll, xul, yll, yul], origin='lower', cmap=plt.cm.jet, alpha=0.9)
# plt.colorbar()
plt.xlabel(r'$r_{3a}$')
plt.ylabel(r'$t$')
plt.clim(0, 1)
#plt.savefig('ERhor3b1.eps', format='eps', dpi=100)
plt.show()

#print((rho3t(1, -1, 1, 1) - Udd(1, 1)*np.kron(rho1qb(0, 0, 1), rho1qb(0, 0, -1))*Dagger(Udd(1, 1))).real)

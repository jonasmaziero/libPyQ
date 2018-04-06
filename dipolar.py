import numpy as np
import matplotlib.pyplot as plt

def Epsi(tha, thb, D, t):
    aa = np.cos(tha/2);  ba = np.sin(tha/2)
    ab = np.cos(thb/2);  bb = np.sin(thb/2)
    f = 2*aa*ba*ab*bb*(np.cos(2*D*t) - np.cos(4*D*t))
    g = (aa**2*bb**2 + ba**2*ab**2)*np.sin(2*D*t) + 2*aa*ba*ab*bb*np.sin(4*D*t)
    return np.sqrt(f**2 + g**2)

xll=0; xul=2*np.pi; yll=0; yul=2*np.pi
x = np.linspace(xll, xul, 80)
y = np.linspace(yll, yul, 80)
X, Y = np.meshgrid(x, y)
#thb = 4*np.pi/8
D = 1
t = 2.5*np.pi/8
Z = Epsi(X, Y, D, t)
contours = plt.contour(X, Y, Z, 6, colors='black')
plt.clabel(contours, inline=True, fontsize=7)
plt.imshow(Z,extent=[xll,xul,yll,yul],origin='lower',cmap=plt.cm.jet,alpha=0.9)
#plt.colorbar()
plt.xlabel(r'$\theta_{a}$')
plt.ylabel(r'$\theta_{b}$')
plt.clim(0,1)
plt.savefig('EPsit25Pi8.eps', format='eps', dpi=100)
plt.show()

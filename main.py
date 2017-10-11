#-----------------------------------------------------------------------------------------------------------------------------------
#import discord
#discord.test_discord()

#import matrix_functions as mf
#mf.test_kp()

'''from states import rho_1qb
rho = rho_1qb(0.0,1.0,0.0)
print(rho)'''

#import lapack
#lapack.test_lapack()

#import distances
#distances.test_distances()

'''
import tomography;  path = "/home/jonasmaziero/Dropbox/Research/IBM_QC/tomography/BDS/c1_005/"
rhoE = tomography.tomo_2qb(path)
CM = [[1.0, 0.0, 0.0, 0.0],[0.0, 0.05, 0.0, 0.0],[0.0, 0.0, -0.05, 0.0], [0.0, 0.0, 0.0, 0.25]]
from states import rho_2qb;  rhoT = rho_2qb(CM)#;  print(rhoE-rhoT)
import distances; print('fidelity = ', distances.fidelity_mm(4, rhoT, rhoE))
from numpy import linalg as LA;  valT = LA.eigvalsh(rhoT);  valE = LA.eigvalsh(rhoE);  
print('eigens of rhoT =',valT);  print('eigens of rhoE =',valE)
#print(rhoT.real);  
import plots;  plots.plot_rho2qb(rhoE.real)
'''

#rho = [[1.0,0.0],[0.0,0.0]]
#import plots;  plots.plot_rho1qb(rho)

'''
from pTranspose import Ta, Tb;  from numpy import zeros;  rho = zeros((4,4), dtype = complex);  rhoTp = zeros((4,4), dtype = complex)
rho = [[0.5,0,0,0.5],[0,0,0,0],[0,0,0,0],[0.5,0,0,0.5]]; rhoTp = Tb(2, 2, rho)#;  print(rhoTp.real)
from entanglement import negativity;  print(negativity(4, rhoTp))
from states import werner;  rho = werner(1.0);  print(rho.real)
'''
#import entanglement;  entanglement.test_entanglement()

import rsvg;  rsvg.test_rsvg()
#-----------------------------------------------------------------------------------------------------------------------------------
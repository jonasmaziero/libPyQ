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

import tomography;  path = "/home/jonasmaziero/Dropbox/Research/IBM_QC/tomography/BDS/c1_06/"
rhoE = tomography.tomo_2qb(path)
CM = [[1.0, 0.0, 0.0, 0.0],[0.0, 0.6, 0.0, 0.0],[0.0, 0.0, -0.6, 0.0], [0.0, 0.0, 0.0, 0.25]]
from states import rho_2qb;  rhoT = rho_2qb(CM)#;  print(rhoE-rhoT)
import distances; print(distances.fidelity_mm(4, rhoT, rhoE))
#print(rhoE.real)
import plots;  plots.plot3dBar(rhoT.real)
#-----------------------------------------------------------------------------------------------------------------------------------
import states
import tomography
import numpy as np
import tomography
rho = np.zeros((4, 4), dtype=complex)

#rho = states.Werner(1)

sj = str(14)
path1 = '/home/jonasmaziero/Dropbox/Research/ibm/bds/'
path2 = 'werner_qx2/dados_plot/'
path = path1 + path2 + sj + '/'
rho = tomography.tomo_2qb(path)

tomography.plot_rho2qb(rho.real)
# tomography.plot_rho2qb(rho.imag)

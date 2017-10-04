#------------------------------------------------------------------------------------------------------------------------------------
def tomo_2qb(path):
  from numpy import genfromtxt, zeros
  ns = 8192.0
  CM = zeros((4,4));  CM[0][0] = 1.0
  #path = "/Users/jonasmaziero/Dropbox/Research/IBM_QC/scp_tel/experiment/tomography_BDS/c1_06/"
  fname = path + "XX.csv"
  pXX = genfromtxt(fname, delimiter = ",", skip_header = 1);  CM[1][1] = ((pXX[0][1] + pXX[3][1]) - (pXX[1][1] + pXX[2][1]))/ns
  CM[1][0] = ((pXX[0][1] + pXX[2][1]) - (pXX[1][1] + pXX[3][1]))/ns
  CM[0][1] = ((pXX[0][1] + pXX[1][1]) - (pXX[2][1] + pXX[3][1]))/ns
  fname = path + "XY.csv"
  pXY = genfromtxt(fname, delimiter = ",", skip_header = 1);  CM[1][2] = ((pXY[0][1] + pXY[3][1]) - (pXY[1][1] + pXY[2][1]))/ns
  fname = path + "XZ.csv"
  pXZ = genfromtxt(fname, delimiter = ",", skip_header = 1);  CM[1][3] = ((pXZ[0][1] + pXZ[3][1]) - (pXZ[1][1] + pXZ[2][1]))/ns
  fname = path + "YX.csv"
  pYX = genfromtxt(fname, delimiter = ",", skip_header = 1);  CM[2][1] = ((pYX[0][1] + pYX[3][1]) - (pYX[1][1] + pYX[2][1]))/ns
  fname = path + "YY.csv"
  pYY = genfromtxt(fname, delimiter = ",", skip_header = 1);  CM[2][2] = ((pYY[0][1] + pYY[3][1]) - (pYY[1][1] + pYY[2][1]))/ns
  CM[2][0] = ((pYY[0][1] + pYY[2][1]) - (pYY[1][1] + pYY[3][1]))/ns
  CM[0][2] = ((pYY[0][1] + pYY[1][1]) - (pYY[2][1] + pYY[3][1]))/ns
  fname = path + "YZ.csv"
  pYZ = genfromtxt(fname, delimiter = ",", skip_header = 1);  CM[2][3] = ((pYZ[0][1] + pYZ[3][1]) - (pYZ[1][1] + pYZ[2][1]))/ns
  fname = path + "ZX.csv"
  pZX = genfromtxt(fname, delimiter = ",", skip_header = 1);  CM[3][1] = ((pZX[0][1] + pZX[3][1]) - (pZX[1][1] + pZX[2][1]))/ns
  fname = path + "ZY.csv"
  pZY = genfromtxt(fname, delimiter = ",", skip_header = 1);  CM[3][2] = ((pZY[0][1] + pZY[3][1]) - (pZY[1][1] + pZY[2][1]))/ns
  fname = path + "ZZ.csv"
  pZZ = genfromtxt(fname, delimiter = ",", skip_header = 1);  CM[3][3] = ((pZZ[0][1] + pZZ[3][1]) - (pZZ[1][1] + pZZ[2][1]))/ns
  CM[3][0] = ((pZZ[0][1] + pZZ[2][1]) - (pZZ[1][1] + pZZ[3][1]))/ns
  CM[0][3] = ((pZZ[0][1] + pZZ[1][1]) - (pZZ[2][1] + pZZ[3][1]))/ns
  from states import rho_2qb;  rho = rho_2qb(CM)
  return rho
#------------------------------------------------------------------------------------------------------------------------------------

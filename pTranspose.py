#------------------------------------------------------------------------------------------------------------------------------------
# Returns its PARTIAL TRANSPOSE with relation to system A
def Ta(da, db, rho): 
  from numpy import zeros;  d = da*db;  rhoTa = zeros((d,d), dtype = complex)
  for ja in range(0,da):
    for ka in range(0,da):
      for jb in range(0,db):
        for kb in range(0,db):
          rhoTa[ka*db+jb][ja*db+kb] = rho[ja*db+jb][ka*db+kb]
  return rhoTa
#------------------------------------------------------------------------------------------------------------------------------------
# Returns its PARTIAL TRANSPOSE with relation to system A
def Tb(da, db, rho): 
  from numpy import zeros;  d = da*db;  rhoTb = zeros((d,d), dtype = complex)
  for ja in range(0,da):
    for ka in range(0,da):
      for jb in range(0,db):
        for kb in range(0,db):
          rhoTb[ja*db+kb][ka*db+jb] = rho[ja*db+jb][ka*db+kb]
  return rhoTb
#------------------------------------------------------------------------------------------------------------------------------------
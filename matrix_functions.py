#------------------------------------------------------------------------------------------------------------------------------------
# returns <phi|A|psi>
def sandwich(d, phi, A, psi):
  from numpy import zeros, conj
  csi = zeros(d, dtype = complex);  csi = conj(phi)
  sdw = 0.0
  for j in range(0,d):
    for k in range(0,d):
      sdw = sdw + csi[j]*A[j][k]*psi[k]
  return sdw
#------------------------------------------------------------------------------------------------------------------------------------
# returns |psi><phi|
def outer(d, psi, phi):
  from numpy import zeros, conj
  csi = zeros(d, dtype = complex);  csi = conj(phi);  ot = zeros((d,d), dtype = complex)
  for j in range(0,d):
    for k in range(0,d):
      ot[j][k] = psi[j]*csi[k]
  return ot
#------------------------------------------------------------------------------------------------------------------------------------
'''def test_kp():
  ra = 2; ca = 2; rb = 2; cb = 2; r = ra*rb; c = ca*cb;
  from numpy import zeros, kron, real
  A = zeros((ra,ca))
  B = zeros((rb,cb))
  kp = zeros((r,c))
  A[0][0] = 1.0; A[0][1] = 2.0; A[1][0] = 3.0; A[1][1] = 4.0
  B[0][0] = 1.0; B[0][1] = 1.0; B[1][0] = 1.0; B[1][1] = 1.0
  #kp = kronecker_product(ra, ca, rb, cb, A, B);
  kp = kron(A,B)
  for k in range(0,r):
    print(real(kp[k][0]), real(kp[k][1]), real(kp[k][2]), real(kp[k][3]))
#------------------------------------------------------------------------------------------------------------------------------------
# Returns the tensor product of two general complex matrices (this function is already implemented in numpy)
def kronecker_product(ra, ca, rb, cb, A, B):
# ra, ca, rb, cb  # Number of rows and columns of the two matrices
# M1(ra,ca), M2(rb,cb)  ! Matrices to take the tensor product of
# M1_kp_M2(ra*rb,ca*cb)  ! Matrix containing the tensor product of M1 and M2
  #int ja, ka, jb, kb, j, k, jl, ju, kl, ku;  // Auxiliary variables for counters
 from numpy import zeros
 kp = zeros((ra*rb,rb*cb))
  for ja in range(0,ra):
    jl = ja*rb
    ju = jl + rb
    for ka in range(0,ca):
      if (A[ja][ka] != 0.0):
        jb = -1;
        kl = ka*cb
        ku = kl + cb
        for j in range(jl,ju): 
          jb += 1
          kb = -1
          for k in range(kl,ku):
            kb += 1
            kp[j][k] = A[ja][ka]*B[jb][kb]
  return kp'''
#-----------------------------------------------------------------------------------------------------------------------------------
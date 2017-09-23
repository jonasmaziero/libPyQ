#-----------------------------------------------------------------------------------------------------------------------------------
import numpy as np
from math import sqrt
#-----------------------------------------------------------------------------------------------------------------------------------
# Returns the two-qubit evolved state for local PHASE DAMPING channels
def rho_pd(p, rho):
  K0 = np.zeros((2,2));  K0[0][0] = sqrt(1.0-p);  K0[1][1] = K0[0][0]
  K1 = np.zeros((2,2));  K1[0][0] = sqrt(p) 
  K2 = np.zeros((2,2));  K2[1][1] = K1[0][0]
  rhop = np.zeros((4,4));  tp = np.zeros((4,4))
  tp = np.kron(K0,K0);  rhop = rhop + np.matmul(np.matmul(tp,rho),tp)
  tp = np.kron(K0,K1);  rhop = rhop + np.matmul(np.matmul(tp,rho),tp)
  tp = np.kron(K0,K2);  rhop = rhop + np.matmul(np.matmul(tp,rho),tp)
  tp = np.kron(K1,K0);  rhop = rhop + np.matmul(np.matmul(tp,rho),tp)
  tp = np.kron(K1,K1);  rhop = rhop + np.matmul(np.matmul(tp,rho),tp)
  tp = np.kron(K1,K2);  rhop = rhop + np.matmul(np.matmul(tp,rho),tp)
  tp = np.kron(K2,K0);  rhop = rhop + np.matmul(np.matmul(tp,rho),tp)
  tp = np.kron(K2,K1);  rhop = rhop + np.matmul(np.matmul(tp,rho),tp)
  tp = np.kron(K2,K2);  rhop = rhop + np.matmul(np.matmul(tp,rho),tp)
  return rhop
#-----------------------------------------------------------------------------------------------------------------------------------

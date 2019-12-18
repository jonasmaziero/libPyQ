import numpy as np
import math
import mat_func as mf
import states as st

def pauli(j):
    if j == 0:
        return np.array([[1,0],[0,1]])
    elif j == 1:
        return np.array([[0,1],[1,0]])
    elif j == 2:
        return np.array([[0,-1j],[1j,0]])
    elif j == 3:
        return np.array([[1,0],[0,-1]])

def hadamard():
    return np.array([[1,1],[1,-1]])/math.sqrt(2)

def O2(x):
    return np.array([[math.cos(x),-math.sin(x)],[math.sin(x),math.cos(x)]])

def id(d):
    id = np.zeros((d,d))
    for j in range(0,d):
        id[j][j] = 1
    return id

def cnot(n,c,t):
    '''returns the control-NOT for n qubits, control c & target t'''
    list1 = []; list2 = []
    for j in range(0,n):
        if j == c:
            list1.append(mf.proj(2,st.cb(2,0)))
            list2.append(mf.proj(2,st.cb(2,1)))
        elif j == t:
            list1.append(pauli(0))
            list2.append(pauli(1))
        else:
            list1.append(pauli(0))
            list2.append(pauli(0))
    kp1 = np.kron(list1[0],list1[1]); kp2 = np.kron(list2[0],list2[1])
    for j in range(2,n):
        kp1 = np.kron(kp1,list1[j])
        kp2 = np.kron(kp2,list2[j])
    cn = kp1+kp2
    return cn

#print(cnot(2,1,0).real)

def test_lapack():
    import scipy.linalg.lapack as lapak
    import numpy as np
    A = np.zeros((2, 2), dtype=complex)
    A[0][0] = 0.0
    A[0][1] = 0.0 - 1j
    A[1][0] = 0.0 + 1j
    A[1][1] = 0.0
    eig = lapak.zheevd(A)
    # in the first array are the eigenvalues
    print(eig[0][0], eig[0][1])
    # and in the 2nd array are the eigenvector, as the columns
    for j in range(0, 2):
        print(eig[1][j][0])
    print(eig[1][0][0], eig[1][1][0])

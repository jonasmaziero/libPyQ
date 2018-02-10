def tomo_2qb(path):
    from numpy import genfromtxt, zeros
    ns = 8192.0
    CM = zeros((4, 4))
    CM[0][0] = 1.0
    #path = "/Users/jonasmaziero/Dropbox/Research/IBM_QC/scp_tel/experiment/tomography_BDS/c1_06/"
    fname = path + "XX.csv"
    pXX = genfromtxt(fname, delimiter=",", skip_header=1)
    CM[1][1] = ((pXX[0][1] + pXX[3][1]) - (pXX[1][1] + pXX[2][1]))/ns
    CM[1][0] = ((pXX[0][1] + pXX[2][1]) - (pXX[1][1] + pXX[3][1]))/ns
    CM[0][1] = ((pXX[0][1] + pXX[1][1]) - (pXX[2][1] + pXX[3][1]))/ns
    fname = path + "XY.csv"
    pXY = genfromtxt(fname, delimiter=",", skip_header=1)
    CM[1][2] = ((pXY[0][1] + pXY[3][1]) - (pXY[1][1] + pXY[2][1]))/ns
    fname = path + "XZ.csv"
    pXZ = genfromtxt(fname, delimiter=",", skip_header=1)
    CM[1][3] = ((pXZ[0][1] + pXZ[3][1]) - (pXZ[1][1] + pXZ[2][1]))/ns
    fname = path + "YX.csv"
    pYX = genfromtxt(fname, delimiter=",", skip_header=1)
    CM[2][1] = ((pYX[0][1] + pYX[3][1]) - (pYX[1][1] + pYX[2][1]))/ns
    fname = path + "YY.csv"
    pYY = genfromtxt(fname, delimiter=",", skip_header=1)
    CM[2][2] = ((pYY[0][1] + pYY[3][1]) - (pYY[1][1] + pYY[2][1]))/ns
    CM[2][0] = ((pYY[0][1] + pYY[2][1]) - (pYY[1][1] + pYY[3][1]))/ns
    CM[0][2] = ((pYY[0][1] + pYY[1][1]) - (pYY[2][1] + pYY[3][1]))/ns
    fname = path + "YZ.csv"
    pYZ = genfromtxt(fname, delimiter=",", skip_header=1)
    CM[2][3] = ((pYZ[0][1] + pYZ[3][1]) - (pYZ[1][1] + pYZ[2][1]))/ns
    fname = path + "ZX.csv"
    pZX = genfromtxt(fname, delimiter=",", skip_header=1)
    CM[3][1] = ((pZX[0][1] + pZX[3][1]) - (pZX[1][1] + pZX[2][1]))/ns
    fname = path + "ZY.csv"
    pZY = genfromtxt(fname, delimiter=",", skip_header=1)
    CM[3][2] = ((pZY[0][1] + pZY[3][1]) - (pZY[1][1] + pZY[2][1]))/ns
    fname = path + "ZZ.csv"
    pZZ = genfromtxt(fname, delimiter=",", skip_header=1)
    CM[3][3] = ((pZZ[0][1] + pZZ[3][1]) - (pZZ[1][1] + pZZ[2][1]))/ns
    CM[3][0] = ((pZZ[0][1] + pZZ[2][1]) - (pZZ[1][1] + pZZ[3][1]))/ns
    CM[0][3] = ((pZZ[0][1] + pZZ[1][1]) - (pZZ[2][1] + pZZ[3][1]))/ns
    # print(CM.real)
    from states import rho_2qb
    rho = rho_2qb(CM)
    return rho


def plot_rho2qb(rho):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    result = [rho[0][0], rho[0][1], rho[0][2], rho[0][3],
              rho[1][0], rho[1][1], rho[1][2], rho[1][3],
              rho[2][0], rho[2][1], rho[2][2], rho[2][3],
              rho[3][0], rho[3][1], rho[3][2], rho[3][3]]
    result = np.array(result, dtype=np.float)
    fig = plt.figure(figsize=(5, 5), dpi=150)
    ax1 = fig.add_subplot(111, projection='3d')
    xlabels = np.array([r'$|00\rangle$', r'$|01\rangle$', r'$|10\rangle$', r'$|11\rangle$'])
    xpos = np.arange(xlabels.shape[0])
    ylabels = np.array([r'$|00\rangle$', r'$|01\rangle$', r'$|10\rangle$', r'$|11\rangle$'])
    ypos = np.arange(ylabels.shape[0])
    xposM, yposM = np.meshgrid(xpos, ypos, copy=False)
    zpos = result
    dx = 0.5
    dy = 0.5
    dz = zpos
    ax1.w_xaxis.set_ticks(xpos + dx/2.0)
    ax1.w_xaxis.set_ticklabels(xlabels)
    ax1.w_yaxis.set_ticks(ypos + dy/2.0)
    ax1.w_yaxis.set_ticklabels(ylabels)
    values = np.linspace(0.2, 1.0, xposM.ravel().shape[0])
    colors = cm.rainbow(values)
    ax1.bar3d(xposM.ravel(), yposM.ravel(), dz*0, dx, dy, dz, color=colors)
    plt.show()


def tomo_1qb(path):
    from numpy import genfromtxt
    ns = 8192.0
    fname = path + "X.csv"
    pX = genfromtxt(fname, delimiter=",", skip_header=1)
    r1 = (pX[0][1] - pX[1][1])/ns
    fname = path + "Y.csv"
    pY = genfromtxt(fname, delimiter=",", skip_header=1)
    r2 = (pY[0][1] - pY[1][1])/ns
    fname = path + "Z.csv"
    pZ = genfromtxt(fname, delimiter=",", skip_header=1)
    r3 = (pZ[0][1] - pZ[1][1])/ns
    #print(r1, r2, r3)
    from states import rho_1qb
    rho = rho_1qb(r1, r2, r3)
    return rho


def plot_rho1qb(rho):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    result = [rho[0][0], rho[0][1],
              rho[1][0], rho[1][1]]
    result = np.array(result, dtype=np.float)
    fig = plt.figure(figsize=(5, 5), dpi=150)
    ax1 = fig.add_subplot(111, projection='3d')
    xlabels = np.array([r'$|0\rangle$', r'$|1\rangle$'])
    xpos = np.arange(xlabels.shape[0])
    ylabels = np.array([r'$|1\rangle$', r'$|0\rangle$'])
    ypos = np.arange(ylabels.shape[0])
    xposM, yposM = np.meshgrid(xpos, ypos, copy=False)
    zpos = result
    dx = 0.5
    dy = 0.5
    dz = zpos
    ax1.w_xaxis.set_ticks(xpos + dx/2.0)
    ax1.w_xaxis.set_ticklabels(xlabels)
    ax1.w_yaxis.set_ticks(ypos + dy/2.0)
    ax1.w_yaxis.set_ticklabels(ylabels)
    values = np.linspace(0.2, 1.0, xposM.ravel().shape[0])
    colors = cm.rainbow(values)
    ax1.bar3d(xposM.ravel(), yposM.ravel(), dz*0, dx, dy, dz, color=colors)
    plt.show()

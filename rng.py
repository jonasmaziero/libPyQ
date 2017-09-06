#------------------------------------------------------------------------------------------------------------------------------------
import random
import numpy as np
import plots

def test_rng():
    random.seed()
    d = 10**4
    x = np.zeros(d)
    y = np.zeros(d)
    for j in range(0,d):
        x[j] = random.random()
        y[j] = random.random()
        #print(j, x[j], y[j])
    plots.plotScatter(x,y)
#------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------
from numpy import random;  random.seed()
import numpy as np
import plots
#------------------------------------------------------------------------------------------------------------------------------------
def test():
  # test for 2d distribution
  ns = 10**4
  '''x = np.zeros(ns)
  y = np.zeros(ns)
  for j in range(0,ns):
    x[j] = random.random()
    y[j] = random.random()
  plots.plotScatter(x,y)'''
  
  # test for Gaussian distribution (code copied from numpy)
  mu, sigma = 0.0, 1.0 # mean and standard deviation
  s = np.random.normal(mu, sigma, 1000)
  import matplotlib.pyplot as plt
  count, bins, ignored = plt.hist(s, 30, normed = True)
  plt.plot(bins, (1.0/(sigma*np.sqrt(2*np.pi)))*np.exp(-((bins-mu)**2.0)/(2.0*sigma**2.0) ), linewidth = 2, color = 'r')
  plt.show()
#------------------------------------------------------------------------------------------------------------------------------------


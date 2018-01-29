#-----------------------------------------------------------------------------------------------------------------------------------
from numpy import genfromtxt, zeros
Ns = zeros(2);  #fname = "mnistNs";  Ns = genfromtxt(fname)
Ns[0] = 1000;  Ns[1] = 1000
import pickle, gzip
trainL = open('trainL', 'w');  trainD = open('trainD', 'w')
testL = open('testL', 'w');  testD = open('testD', 'w')
fmnist = gzip.open('mnist.pkl.gz', 'rb')
train_set, valid_set, test_set = pickle.load(fmnist, encoding = 'latin1')
fmnist.close()
for i in range(0,int(Ns[0])):
  trainL.write(str(train_set[1][i]))
  trainL.write('\n')
  for j in range(0,784):
    trainD.write(str(train_set[0][i][j]))
    trainD.write(' ')
  trainD.write('\n')
for i in range(0,int(Ns[1])):
  testL.write(str(test_set[1][i]))
  testL.write('\n')
  for j in range(0,784):
    testD.write(str(test_set[0][i][j]))
    testD.write(' ')
  testD.write('\n')
#-----------------------------------------------------------------------------------------------------------------------------------
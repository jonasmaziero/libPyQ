def googol():
    import random;  random.seed()
    nj = 10**1  # No. of games
    nrn = 10**4  # max int No.
    score = 0.0
    for j in range(0,nj):
        r0 = random.randint(0,nrn)/(1.0*nrn)
        r1 = random.randint(0,nrn)/(1.0*nrn)
        f0 = sigmoid(r0)
        rr = random.uniform(0,1)
        print(r0,r1,f0,rr)
        if f0 > rr:
            gtc = r0
        else:
            gtc = r0
        if r0 > r1:
            gt = r0
        else:
            gt = r1
        if gt == gtc:
            score +=1
    print((1.0*score)/nj)


def sigmoid(x):
    import numpy as np
    return 1.0/(1.0+np.exp(-x))

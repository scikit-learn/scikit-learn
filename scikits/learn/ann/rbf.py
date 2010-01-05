# rbf2.py
# tilde
# 2006/08/18
# - new attempt at RBF net to get my ideas straight...deadline is fast approaching!

from numpy import *

class rbf:

    _type = 'rbf'
    
    def __init__(nin,nhid,nout,trndata):
        # set easy params
        self.nin = nin
        self.nhid = nhid
        self.nout = nout
        # choose subset (1/5?) of training data for basis fxn centers and
        self.centers = []
        for i in trndata:
            if random.random < 0.2:
                self.centers.append(i)
        # set common variance proportional to max dist between centers
        d_max = 0.0
        for i in self.centers:
            for j in self.centers:
                tmp = sqrt((i-j)**2)
                if tmp > d_max:
                    d_max = tmp
        self.variance = d_max/2.0*size(trndata)
        
    
    def fwd(self,inputs):
        """ Propagate values forward through the net.
        Inputs:
                inputs      - vector of input values
        """
        z = exp((-1.0/(2*self.variance))*
        o = dot(z,self.w) + dot(ones((len(z),1)),self.b)


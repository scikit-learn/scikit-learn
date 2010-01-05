# srn.py
# by: Fred Mailhot
# last mod: 2006-06-21

from scipy import *

class srn:
    """Class to define, train and test a simple recurrent network
    a.k.a. 'Elman net' 
    (cf. Elman 1991's Machine Learning paper,inter alia)
    """

    _type = 'srn'
    _outfxns = ('linear','logistic','softmax')
    _algs = ('rtrl','bptt','srbptt')

    def __init__(self,nin,nhid,nout,fxn,alg='srbptt',w=None):
        """ Set up instance of srn. Initial weights are drawn from a 
        zero-mean Gaussian w/ variance is scaled by fan-in.
        (see Bishop 1995 for justification)
        Inputs:
            nin         - integer number of input units
            nhid        - integer number of hiden & context units
            nout        - integer number of output units, 
            fxn         - string description of output unit activation fxn;
                          one of {'linear','logistic','softmax'}
                          (n.b. hidden/context units use tanh)
        """
        if fxn not in self._outfxns:
            print "Undefined activation fxn. Using linear"
            self.outfxn = 'linear'
        else:
            self.outfxn = fxn
        self.nin = nin
        self.nhid = nhid
        self.ncon = nhid    # context units
        self.nout = nout
        self.alg = alg
        if w:
            self.nwts = size(w)
            self.w_packed = w
            self.w1 = zeros((nin,nhid),dtype=Float)
            self.b1 = zeros((1,nhid),dtype=Float)
            self.w_con = zeros((nhid,nhid),dtype=Float)     # context weights
            self.w2 = zeros((nhid,nout),dtype=Float)
            self.b2 = zeros((1,nout),dtype=Float)
            self.unpackwts()
        else:
            self.nwts = (nin+1)*nhid + (nhid+1)*nout
            self.w1 = randn(nin,nhid)/sqrt(nin+1)
            self.b1 = randn(1,nhid)/sqrt(nin+1)
            self.w_con = randn(nhid,nhid)/sqrt(nhid+1)      # context weights
            self.w2 = randn(nhid,nout)/sqrt(nhid+1)
            self.b2 = randn(1,nout)/sqrt(nhid+1)
            self.packwts()

    def unpackwts(self):
        """ Decompose 1-d vector of weights w into appropriate weight 
        matrices (w1,b1,w2,b2) and reinsert them into net
        """
        self.w1 = reshape(array(self.w_packed)[:self.nin*self.nhid],(self.nin,self.nhid))
        self.b1 = reshape(array(self.w_packed)[(self.nin*self.nhid):(self.nin*self.nhid)+self.nhid],(1,self.nhid))
        self.w2 = reshape(array(self.w_packed)[(self.nin*self.nhid)+self.nhid:\
                                (self.nin*self.nhid)+self.nhid+(self.nhid*self.nout)],(self.nhid,self.nout))
        self.b2 = reshape(array(self.w_packed)[(self.nin*self.nhid)+self.nhid+(self.nhid*self.nout):],(1,self.nout))

    def packwts(self):
        """ Compile weight matrices w1,b1,w2,b2 from net into a
        single vector, suitable for optimization routines.
        """
        self.w_packed = hstack([self.w1.reshape(size(self.w1)),
                                self.b1.reshape(size(self.b1)),
                                self.w2.reshape(size(self.w2)),
                                self.b2.reshape(size(self.b2))])

    def fwd(self,inputs,wts=None,hid=False):
        """ Propagate values forward through the net.
        Inputs:
            inputs  - self.nin*1 vector of inputs
            hid     - boolean specifying whether or not to return hidden
                      unit activations, False by default
        """
        if wts is not None:
            self.w_packed = wts
        self.unpackwts()
        
        z = tanh(dot(inputs,self.w1) + dot(ones((len(inputs),1)),self.b1))
        o = dot(z,self.w2) + dot(ones((len(z),1)),self.b2)
        
        if self.outfxn == 'linear':
            y = o
        elif self.outfxn == 'logistic':     # TODO: check for overflow here...
            y = 1/(1+exp(-o))
        elif self.outfxn == 'softmax':      # TODO: and here...
            tmp = exp(o)
            y = tmp/(sum(temp,1)*ones((1,self.nout)))
            
        if hid:
            return array(y),array(z)
        else:
            return array(y)

def main():

    print "this doesn't do much, yet"
    opt = raw_input("\nEnter desired training algorithm: ")
    print "Creating 1-2-1 SRN for 'temporal-XOR'"
    net = srn(1,2,1,'logistic',opt)

if __name__ == '__main__':
    main()


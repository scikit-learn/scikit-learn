# mlp.py
# by: Fred Mailhot

import sys
from scipy import *

class mlp:
    """Class to define, train and test a multilayer perceptron."""

    _type = "mlp"
    _outfxns = ("linear","logistic","softmax")
    _algs = ('leastsq','cg','quasinew')
    _options = zeros(18)

    def __init__(self,nin,nhid,nout,fxn,w=None):
        """Set up instance of mlp. Initial weights are drawn from a zero-mean,
        unit-variance Gaussian (variance is scaled by fan-in).
            Parameters:
                nin/nhid/nout - integer number of input/hidden/output units, 
                                respectively
                fxn           - string description of output unit activation 
                                fxn (hidden units use tanh); can be "linear", 
                                "logistic" or "softmax" 
        """
        if fxn not in self._outfxns:
            print "Undefined error fxn. Aborting."
            sys.exit(1)
        else:
            self.outfxn = fxn
        # units
        self.nin = nin
        self.nhid = nhid
        self.nout = nout
        # weights
        if w is not None:
            if size(w) != (nin+1)*nhid + (nhid+1)*nout:
                print "Incorrectly sized weight vector. Aborting."
                sys.exit(1)
            self.nwts = size(w)
            self.w_packed = w
            self.w1 = zeros((nin,nhid))
            self.b1 = zeros((1,nhid))
            self.w2 = zeros((nhid,nout))
            self.b2 = zeros((1,nout))
            self.unpackwts()
        else:
            self.nwts = (nin+1)*nhid + (nhid+1)*nout
            self.w1 = randn(nin,nhid)/sqrt(nin+1)
            self.b1 = randn(1,nhid)/sqrt(nin+1)
            self.w2 = randn(nhid,nout)/sqrt(nhid+1)
            self.b2 = randn(1,nout)/sqrt(nhid+1)
            self.packwts()
        self.options = self._options

    def unpackwts(self):
        """Decompose 1-d vector of weights w into appropriate weight 
            matrices (w1,b1,w2,b2) and reinsert them into net
        """
        self.w1 = reshape(self.w_packed[:self.nin*self.nhid],(self.nin,self.nhid))
        self.b1 = reshape(self.w_packed[(self.nin*self.nhid):(self.nin*self.nhid)+self.nhid],(1,self.nhid))
        self.w2 = reshape(self.w_packed[(self.nin*self.nhid)+self.nhid: \
                            (self.nin*self.nhid)+self.nhid+(self.nhid*self.nout)],(self.nhid,self.nout))
        self.b2 = reshape(self.w_packed[(self.nin*self.nhid)+self.nhid+(self.nhid*self.nout):],(1,self.nout))

    def packwts(self):
        """Compile weight matrices w1,b1,w2,b2 from net into a
            single vector, suitable for optimization functions.
        """
        self.w_packed = hstack([reshape(self.w1,(size(self.w1),)),
                                reshape(self.b1,(size(self.b1),)),
                                reshape(self.w2,(size(self.w2),)),
                                reshape(self.b2,(size(self.b2),))])

    def fwd(self,inputs,w):
        if self.w_packed != w:  # ADDED THIS TO MAKE MY MODEL EXPLICITLY
            self.w_packed = w   # DEPENDENT ON INPUTS AND PARAMETERS
        self.unpackwts()
        z = tanh(array(matrix(inputs)*matrix(self.w1) + matrix(ones((len(inputs),1)))*matrix(self.b1)))
        o = array(matrix(z)*matrix(self.w2) + matrix(ones((len(z),1)))*matrix(self.b2))
        if self.outfxn == 'linear':
            y = o
        elif self.outfxn == 'logistic':     # TODO: check for overflow here and...
            y = 1/(1+exp(-o))
        elif self.outfxn == 'softmax':      # TODO: here...
            tmp = exp(o)
            y = tmp/(sum(temp,1)*ones((1,self.nout)))
        else:   # this shouldn't happen
            print "Undefined output fxn. Uh oh."
            sys.exit(1)
        return y

#def residuals(w,x,t,i,h,o,f):
#        n = mlp(i,h,o,f,w)
#        err = t - n.fwd(x)
#        return sum(err,axis=1)

def residuals(w,t,x,n): 
    err = t - n.fwd(x)
    return sum(err,axis=1)

def mlptrain(net,x,t,alg):
    """ Train a multilayer perceptron with a user-specified algorithm.
        Parameters:
            x   - matrix of input data
            t   - matrix of target outputs
            alg - training algorithm, can be "cg" (conjugate gradients), 
                        "qn" (quasi-newton), or "gd" (batch gradient descent)
    """
    #if alg in self._algs:
    from scipy.optimize import leastsq
    w = leastsq(residuals,net.w_packed,args=(t,x,net),full_output=True)
    #else:
    #    print "Undefined optimization fxn. Aborting."
    #    sys.exit(1)
    return (w[0],w[1:])

def mlptest(net,x,t):
    
    e = t - net.fwd(x)
    return array(sum(matrix(e)*matrix(e).T))

def main():
    """ Approx test of module, using the oilTrn/oilTst data files that are distributed
        with Netlab (the Matlab ANN toolkit). Prints sum-squared-err before and after
        optimization.
        
        ** Something's not working the way I want it to here. The leastsq fxn is
        exceeding it's max number of iterations.**
    """
    from scipy.io import read_array, write_array
    print "Creating 12-12-2 MLP with linear outputs"
    net = mlp(12,12,2,'linear')
    net.packwts()
    w_old = net.w_packed
    print "Loading train and test sets."
    trn_input = read_array('oilTrn.dat',lines=(3,-1),columns=(0,(1,12)))
    trn_targs = read_array('oilTrn.dat',lines=(3,-1),columns=(12,-1))
    tst_input = read_array('oilTst.dat',lines=(3,-1),columns=(0,(1,12)))
    tst_targs = read_array('oilTst.dat',lines=(3,-1),columns=(12,-1))
    print "Testing initialized config...",
    sqrd_err = mlptest(net,tst_input,tst_targs)
    print "Done."
    print "Training net...",
    retval = mlptrain(net,trn_input,trn_targs,"leastsq")[0]
    net.w_packed = retval[0]
    print "Done."
    print retval[1]
    print "Testing trained net..."
    sqrd_err = mlptest(net,tst_input,tst_targs)
    print "Done."

if __name__ == '__main__':
    main()


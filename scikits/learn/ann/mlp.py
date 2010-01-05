# mlp.py
# by: Fred Mailhot
# last mod: 2006-08-18

from scipy import * # I'll want to change this for numpy eventually
from scipy.optimize import leastsq
import copy

class mlp:
    """Class to define, train and test a multilayer perceptron."""

    _type = 'mlp'
    _outfxns = ('linear','logistic','softmax')
    _algs = ('simplex','powell','bfgs','cg','ncg','leastsq')

    def __init__(self,nin,nhid,nout,fxn,alg='leastsq',w=None):
        """ Set up instance of mlp. Initial weights are drawn from a 
        zero-mean Gaussian w/ variance is scaled by fan-in.
        (see Bishop 1995 for justification)
        Inputs:
            nin/nhid/nout - integer number of input/hidden/output units, 
                                respectively
            fxn           - string description of output unit activation 
                                fxn (hidden units use tanh); can be 'linear', 
                                'logistic' or 'softmax' 
        """
        if fxn not in self._outfxns:
            print "Undefined activation fxn. Using linear"
            self.outfxn = 'linear'
        else:
            self.outfxn = fxn
        self.nin = nin
        self.nhid = nhid
        self.nout = nout
        self.alg = alg
        if w:
            self.nwts = size(w)
            self.w_packed = w
            self.w1 = zeros((nin,nhid),dtype=Float)
            self.b1 = zeros((1,nhid),dtype=Float)
            self.w2 = zeros((nhid,nout),dtype=Float)
            self.b2 = zeros((1,nout),dtype=Float)
            self.unpackwts()
        else:
            self.nwts = (nin+1)*nhid + (nhid+1)*nout
            self.w1 = randn(nin,nhid)/sqrt(nin+1)
            self.b1 = randn(1,nhid)/sqrt(nin+1)
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

    def errfxn(self,w,x,t):
        """ Return vector of squared-errors for the leastsq optimizer
        """
        y = self.fwd(x,w)
        return sum(array(y-t)**2,axis=1)

    def train(self,x,t):
        """ Train a multilayer perceptron using scipy's leastsq optimizer
        Input:
            x   - matrix of input data
            t   - matrix of target outputs
        Returns:
            post-optimization weight vector
        """
        # something's going wrong w/ the full_output option
        # return leastsq(self.errfxn,self.w_packed,args=(x,t),full_output=True)
        return leastsq(self.errfxn,self.w_packed,args=(x,t))

def main():
    """ Approx test of module, using the oilTrn/oilTst data files that are 
    """
    from scipy.io import read_array, write_array
    # build the net
    print "\nCreating 12-5-2 MLP with linear outputs"
    net = mlp(12,5,2,'linear')
    w_init = copy.copy(net.w_packed)
    # prep the train/test data
    print "\nLoading training and test sets...",
    trn_input = read_array('data/oilTrn.dat',lines=(3,-1),columns=(0,(1,12)))
    trn_targs = read_array('data/oilTrn.dat',lines=(3,-1),columns=(12,-1))
    tst_input = read_array('data/oilTst.dat',lines=(3,-1),columns=(0,(1,12)))
    tst_targs = read_array('data/oilTst.dat',lines=(3,-1),columns=(12,-1))
    print "done."
    # initial squared-error
    print "\nInitial SSE on training set: ",\
          sum(net.errfxn(net.w_packed,trn_input,trn_targs))
    print "\nInitial SSE on testing set: ",\
          sum(net.errfxn(net.w_packed,tst_input,tst_targs))
    # train the net
    net.w_packed = net.train(trn_input,trn_targs)[0]
    # final squared-error
    print "\nFinal SSE on training set: ",\
          sum(net.errfxn(net.w_packed,trn_input,trn_targs))
    print "\nFinal SSE on testing set: ",\
          sum(net.errfxn(net.w_packed,tst_input,tst_targs))
    # view extended output?
    # REMOVING THIS OPTION FOR NOW
    #if raw_input("Do you want to see the full training output? (y/n").lower() == 'y':
    #    print retval[1]
        
if __name__ == '__main__':
    main()


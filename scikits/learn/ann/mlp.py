# mlp.py
# by: Fred Mailhot
# last mod: 2006-08-19

import numpy as N
from scipy.optimize import leastsq

class mlp:
    """Class to define, train and test a multilayer perceptron.
    """

    _type = 'mlp'
    _outfxns = ('linear','logistic','softmax')

    def __init__(self,ni,nh,no,f='linear',w=None):
        """ Set up instance of mlp. Initial weights are drawn from a 
        zero-mean Gaussian w/ variance is scaled by fan-in.
        Input:
            ni  - <int> # of inputs
            nh  - <int> # of hidden & context units
            no  - <int> # of outputs
            f   - <str> output activation fxn
            w   - <array dtype=Float> weight vector
        """
        if f not in self._outfxns:
            print "Undefined activation fxn. Using linear"
            self.outfxn = 'linear'
        else:
            self.outfxn = f
        self.ni = ni
        self.nh = nh
        self.no = no
        #self.alg = alg
        if w:
            self.nw = N.size(w)
            self.wp = w
            self.w1 = N.zeros((ni,nh),dtype=Float)
            self.b1 = N.zeros((1,nh),dtype=Float)
            self.w2 = N.zeros((nh,no),dtype=Float)
            self.b2 = N.zeros((1,no),dtype=Float)
            self.unpack()
        else:
            self.nw = (ni+1)*nh + (nh+1)*no
            self.w1 = N.random.randn(ni,nh)/N.sqrt(ni+1)
            self.b1 = N.random.randn(1,nh)/N.sqrt(ni+1)
            self.w2 = N.random.randn(nh,no)/N.sqrt(nh+1)
            self.b2 = N.random.randn(1,no)/N.sqrt(nh+1)
            self.pack()

    def unpack(self):
        """ Decompose 1-d vector of weights w into appropriate weight 
        matrices (w1,b1,w2,b2) and reinsert them into net
        """
        self.w1 = N.array(self.wp)[:self.ni*self.nh].reshape(self.ni,self.nh)
        self.b1 = N.array(self.wp)[(self.ni*self.nh):(self.ni*self.nh)+self.nh].reshape(1,self.nh)
        self.w2 = N.array(self.wp)[(self.ni*self.nh)+self.nh:(self.ni*self.nh)+self.nh+(self.nh*self.no)].reshape(self.nh,self.no)
        self.b2 = N.array(self.wp)[(self.ni*self.nh)+self.nh+(self.nh*self.no):].reshape(1,self.no)

    def pack(self):
        """ Compile weight matrices w1,b1,w2,b2 from net into a
        single vector, suitable for optimization routines.
        """
        self.wp = N.hstack([self.w1.reshape(N.size(self.w1)),
                            self.b1.reshape(N.size(self.b1)),
                            self.w2.reshape(N.size(self.w2)),
                            self.b2.reshape(N.size(self.b2))])

    def fwd_all(self,x,w=None):
        """ Propagate values forward through the net. 
        Input:
            x   - array (size>1) of input patterns
            w   - optional 1-d vector of weights 
        Returns:
            y   - array of outputs for all input patterns
        """
        if w is not None:
            self.wp = w
        self.unpack()
        # compute vector of hidden unit values
        z = N.tanh(N.dot(x,self.w1) + N.dot(N.ones((len(x),1)),self.b1))
        # compute vector of net outputs
        o = N.dot(z,self.w2) + N.dot(N.ones((len(z),1)),self.b2)
        # compute final output activations
        if self.outfxn == 'linear':
            y = o
        elif self.outfxn == 'logistic':     # TODO: check for overflow here...
            y = 1/(1+N.exp(-o))
        elif self.outfxn == 'softmax':      # TODO: and here...
            tmp = N.exp(o)
            y = tmp/(N.sum(temp,1)*N.ones((1,self.no)))
            
        return N.array(y)

    def errfxn(self,w,x,t):
        """ Return vector of squared-errors for the leastsq optimizer
        """
        y = self.fwd_all(x,w)
        return N.sum(N.array(y-t)**2,axis=1)

    def train(self,x,t):
        """ Train network using scipy's leastsq optimizer
        Input:
            x   - array of input data 
            t   - array of targets
            
            N.B. x and t comprise the *entire* collection of training data
            
        Returns:
            post-optimization weight vector
        """
        return leastsq(self.errfxn,self.wp,args=(x,t))

    def test_all(self,x,t):
        """ Test network on an array (size>1) of patterns
        Input:
            x   - array of input data
            t   - array of targets
        Returns:
            sum-squared-error over all data
        """
        return N.sum(self.errfxn(self.wp,x,t))

def main():
    """ Build/train/test MLP 
    """
    from scipy.io import read_array, write_array
    print "\nCreating 2-2-1 MLP with logistic outputs"
    net = mlp(2,2,1,'logistic')
    print "\nLoading training and test sets...",
    trn_input = read_array('data/xor-trn.dat',lines=(3,-1),columns=(0,(1,2)))
    trn_targs = read_array('data/xor-trn.dat',lines=(3,-1),columns=(2,-1))
    trn_targs = trn_targs.reshape(N.size(trn_targs),1)
    tst_input = read_array('data/xor-tst.dat',lines=(3,-1),columns=(0,(1,2)))
    tst_targs = read_array('data/xor-tst.dat',lines=(3,-1),columns=(2,-1))
    tst_targs = tst_targs.reshape(N.size(tst_targs),1)
    print "done."
    print "\nInitial SSE:\n"
    print "\ttraining set: ",net.test_all(trn_input,trn_targs)
    print "\ttesting set: ",net.test_all(tst_input,tst_targs),"\n"
    net.wp = net.train(trn_input,trn_targs)[0]
    print "\nFinal SSE:\n"
    print "\ttraining set: ",net.test_all(trn_input,trn_targs)
    print "\ttesting set: ",net.test_all(tst_input,tst_targs),"\n"
        
if __name__ == '__main__':
    main()


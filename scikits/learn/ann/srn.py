# srn.py
# by: Fred Mailhot
# last mod: 2006-08-18

import numpy as N
from scipy.optimize import leastsq

class srn:
    """Class to define, train and test a simple recurrent network
    """

    _type = 'srn'
    _outfxns = ('linear','logistic','softmax')

    def __init__(self,ni,nh,no,f='linear',w=None):
        """ Set up instance of srn. Initial weights are drawn from a 
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
        self.nc = nh
        self.no = no
        if w:
            self.nw = N.size(w)
            self.wp = w
            self.w1 = N.zeros((ni,nh),dtype=Float)    # input-hidden wts
            self.b1 = N.zeros((1,nh),dtype=Float)     # input biases
            self.wc = N.zeros((nh,nh),dtype=Float)    # context wts
            self.w2 = N.zeros((nh,no),dtype=Float)    # hidden-output wts
            self.b2 = N.zeros((1,no),dtype=Float)     # hidden biases
            self.unpack()
        else:
            # N.B. I just understood something about the way reshape() works
            # that should simplify things, allowing me to just make changes
            # to the packed weight vector, and using "views" for the fwd
            # propagation.
            # I'll implement this next week.
            self.nw = (ni+1)*nh + (nh*nh) + (nh+1)*no
            self.w1 = N.random.randn(ni,nh)/N.sqrt(ni+1)
            self.b1 = N.random.randn(1,nh)/N.sqrt(ni+1)
            self.wc = N.random.randn(nh,nh)/N.sqrt(nh+1)
            self.w2 = N.random.randn(nh,no)/N.sqrt(nh+1)
            self.b2 = N.random.randn(1,no)/N.sqrt(nh+1)
            self.pack()

    def unpack(self):
        """ Decompose 1-d vector of weights w into appropriate weight 
        matrices (w1,b1,w2,b2) and reinsert them into net
        """
        self.w1 = N.array(self.wp)[:self.ni*self.nh].reshape(self.ni,self.nh)
        self.b1 = N.array(self.wp)[(self.ni*self.nh):(self.ni*self.nh)+self.nh].reshape(1,self.nh)
        self.wc = N.array(self.wp)[(self.ni*self.nh)+self.nh:(self.ni*self.nh)+self.nh+(self.nh*self.nh)].reshape(self.nh,self.nh)
        self.w2 = N.array(self.wp)[(self.ni*self.nh)+self.nh+(self.nh*self.nh):(self.ni*self.nh)+self.nh+(self.nh*self.nh)+(self.nh*self.no)].reshape(self.nh,self.no)
        self.b2 = N.array(self.wp)[(self.ni*self.nh)+self.nh+(self.nh*self.nh)+(self.nh*self.no):].reshape(1,self.no)

    def pack(self):
        """ Compile weight matrices w1,b1,wc,w2,b2 from net into a
        single vector, suitable for optimization routines.
        """
        self.wp = N.hstack([self.w1.reshape(N.size(self.w1)),
                            self.b1.reshape(N.size(self.b1)),
                            self.wc.reshape(N.size(self.wc)),
                            self.w2.reshape(N.size(self.w2)),
                            self.b2.reshape(N.size(self.b2))])

    def fwd_all(self,x,w=None):
        """ Propagate values forward through the net. 
        Input:
            x   - matrix of all input patterns
            w   - 1-d vector of weights
        Returns:
            y   - matrix of all outputs
        """
        if w is not None:
            self.wp = w
        self.unpack()
        # compute vector of context values for current weight matrix
        c = N.tanh(N.dot(x,self.w1) + N.dot(N.ones((len(x),1)),self.b1))
        c = N.vstack([c[1:],c[0]])
        # compute vector of hidden unit values
        z = N.tanh(N.dot(x,self.w1) + N.dot(c,self.wc) + N.dot(N.ones((len(x),1)),self.b1))
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
        """ Train a multilayer perceptron using scipy's leastsq optimizer
        Input:
            x   - matrix of input data
            t   - matrix of target outputs
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
    """ Set up a 1-2-1 SRN to solve the temporal-XOR problem from Elman 1990.
    """
    from scipy.io import read_array, write_array
    print "\nCreating 1-2-1 SRN for 'temporal-XOR'"
    net = srn(1,2,1,'logistic')
    print "\nLoading training and test sets...",
    trn_input = read_array('data/txor-trn.dat')
    trn_targs = N.hstack([trn_input[1:],trn_input[0]])
    trn_input = trn_input.reshape(N.size(trn_input),1)
    trn_targs = trn_targs.reshape(N.size(trn_targs),1)
    tst_input = read_array('data/txor-tst.dat')
    tst_targs = N.hstack([tst_input[1:],tst_input[0]])
    tst_input = tst_input.reshape(N.size(tst_input),1)
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


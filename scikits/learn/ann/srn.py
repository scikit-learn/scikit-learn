# srn.py
# by: Fred Mailhot
# last mod: 2006-06-22

from scipy import *
import copy

class srn:
    """Class to define, train and test a simple recurrent network
    a.k.a. 'Elman net' (cf. Elman 1991's Machine Learnig paper,inter alia)
    """

    _type = 'srn'
    _outfxns = ('linear','logistic','softmax')
    _alg = ('srn')

    def __init__(self,ni,nh,no,f,h=-1,w=None):
        """ Set up instance of srn. Initial weights are drawn from a 
        zero-mean Gaussian w/ variance is scaled by fan-in.
        (see Bishop 1995 for justification)
        Inputs:
            ni  - integer number of input units
            nh  - integer number of hiden & context units
            no  - integer number of output units, 
            f   - string description of output unit activation fxn;
                    one of {'linear','logistic','softmax'}
                    (n.b. hidden/context units use tanh)
            w   - initialized 1-d weight vector
        """
        if f not in self._outfxns:
            print "Undefined activation fxn. Using linear"
            self.outfxn = 'linear'
        else:
            self.outfxn = f
        # set up layers of units
        self.ni = ni
        self.nh = nh
        self.nc = nh
        self.no = no
        self.z = zeros((h,nh),dtype=Float)        # hidden activations for 1 epoch
        self.c = zeros((h,nh),dtype=Float)        # context activations for 1 epoch
        self.o = zeros((h,no),dtype=Float)       # output activiation for 1 epoch
        self.p = zeros((nh,nw,nw),dtype=Float)
        if w:
            self.nw = size(w)
            self.wp = w
            self.w1 = zeros((ni,nh),dtype=Float)    # input-hidden wts
            self.b1 = zeros((1,nh),dtype=Float)     # input biases
            self.wc = zeros((nh,nh),dtype=Float)    # context wts
            self.w2 = zeros((nh,no),dtype=Float)    # hidden-output wts
            self.b2 = zeros((1,no),dtype=Float)     # hidden biases
            self.unpack()
        else:
            # N.B. I just understood something about the way reshape() works
            # that should simplify things, allowing me to just make changes
            # to the packed weight vector, and using "views" for the fwd
            # propagation.
            # I'll implement this next week.
            self.nw = (ni+1)*nh + (nh*nh) + (nh+1)*no
            self.w1 = randn(ni,nh)/sqrt(ni+1)
            self.b1 = randn(1,nh)/sqrt(ni+1)
            self.wc = randn(nh,nh)/sqrt(nh+1)
            self.w2 = randn(nh,no)/sqrt(nh+1)
            self.b2 = randn(1,no)/sqrt(nh+1)
            self.pack()
            if size(self.wp) != self.nw:
               raise ValueError, "Unexpected number of weights" 

    def unpack(self):
        """ Decompose 1-d vector of weights w into appropriate weight 
        matrices (w1,b1,w2,b2) and reinsert them into net
        """
        self.w1 = reshape(array(self.wp)[:self.ni*self.nh],(self.ni,self.nh))
        self.b1 = reshape(array(self.wp)[(self.ni*self.nh):(self.ni*self.nh)+self.nh],(1,self.nh))
        self.wc = reshape(array(self.wp)[(self.ni*self.nh)+self.nh:\
                                (self.ni*self.nh)+self.nh+(self.nh*self.nh)],(self.nh,self.nh))
        self.w2 = reshape(array(self.wp)[(self.ni*self.nh)+self.nh+(self.nh*self.nh):\
                                (self.ni*self.nh)+self.nh+(self.nh*self.nh)+(self.nh*self.no)],(self.nh,self.no))
        self.b2 = reshape(array(self.wp)[(self.ni*self.nh)+self.nh+(self.nh*self.no):],(1,self.no))

    def pack(self):
        """ Compile weight matrices w1,b1,wc,w2,b2 from net into a
        single vector, suitable for optimization routines.
        """
        self.wp = hstack([self.w1.reshape(size(self.w1)),
                                self.b1.reshape(size(self.b1)),
                                self.wc.reshape(size(self.wc)),
                                self.w2.reshape(size(self.w2)),
                                self.b2.reshape(size(self.b2))])

    def fwd(self,x,w=None,hid=False):
        """ Propagate values forward through the net. 
        This involves the following steps:
        (i) feeds the current input and context values to the hidden layer, 
        (ii) hidden layer net input is transformed and then sent to the outputs
        (iii) output values are copied to the context layer
        Inputs:
            x   - matrix of all input patterns
            w   - 1-d vector of weights
            hid - boolean specifying whether or not to return hidden
                      unit activations, False by default
        Outputs:
            y   - matrix of all outputs
            z   - matrix of all hidden activations (if hid=True)
        """
        if wts is not None:
            self.wp = w
        self.unpack()
        # compute net input to hiddens and then squash it
        self.z = tanh(dot(x,self.w1) + dot(self.c,self.wc) + dot(ones((len(x),1)),self.b1))
        # send hidden vals to output and copy to context
        o = dot(self.z,self.w2) + dot(ones((len(self.z),1)),self.b2)
        self.c = copy.copy(self.z)
        # compute output activations
        if self.outfxn == 'linear':
            y = o
        elif self.outfxn == 'logistic':     # TODO: check for overflow here...
            y = 1/(1+exp(-o))
        elif self.outfxn == 'softmax':      # TODO: and here...
            tmp = exp(o)
            y = tmp/(sum(temp,1)*ones((1,self.no)))
            
        if hid:
            return array(y),array(z)
        else:
            return array(y)

    def train(self,x,t,N):
        """ Train net by standard backpropagation
        Inputs:
            x   - all input patterns
            t   - all target patterns
            N   - number of times to go over patterns
        Outputs:
            w   - new weight vector
        """
        for i in range(N):
            

    def errfxn(self,w,x,t):
        """ Error functions for each of the output-unit activation functions.
        Inputs:
            w   - current weight vector
            x   - current pattern input(s) (len(x) == self.h)
            t   - current pattern target(s)
        """
        y,z = self.fwd(w,x,True)
        if self.outfxn == 'linear':
            # calculate & return SSE
            err = 0.5*sum(sum(array(y-t)**2,axis=1))
        elif self.outfxn == 'logistic':
            # calculate & return x-entropy
            err = -1.0*sum(sum(t*log2(y)+(1-t)*log2(1-y),axis=1))
        elif self.outfxn == 'softmax':
            # calculate & return entropy
            err = -1.0*sum(sum(t*log2(y),axis=1))
        else:
            # this shouldn't happen, return SSE as safe default
            err = 0.5*sum(sum(array(y-t)**2,axis=1))
        
        # returning a tuple of info for now...not sure why
        return err,y,z

def main():
    """ Set up a 1-2-1 SRN to solve the temporal-XOR problem from Elman 1990.
    """
    from scipy.io import read_array, write_array
    print "Creating 1-2-1 SRN for 'temporal-XOR' (net.h = 2)"
    net = srn(1,2,1,'logistic',2)
    print net
    print "\nLoading training and test sets...",
    trn_input = read_array('data/t-xor1.dat')
    trn_targs = hstack([trn_input[1:],trn_input[0]])
    tst_input = read_array('data/t-xor2.dat')
    tst_targs = hstack([tst_input[1:],tst_input[0]])
    print "done."
    N = input("Number of iterations over training set: ")
    
    print "\nInitial error: ",net.errfxn(net.wp,tst_input,tst_targs)
    net.train(trn_input,trn_targs,N)
    print "\nFinal error: ",net.errfxn(net.wp,tst_input,tst_targs)

if __name__ == '__main__':
    main()


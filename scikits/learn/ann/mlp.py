# mlp.py
# by: Fred Mailhot
# last mod: 2006-06-21

from scipy import *

class mlp:
    """Class to define, train and test a multilayer perceptron."""

    _type = 'mlp'
    _outfxns = ('linear','logistic','softmax')
    _algs = ('simplex','powell','bfgs','cg','ncg','leastsq')
    _options = zeros(18)                # don't know if I'll use this or not

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
        self.options = self._options

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
        """ Implementing 'canonical' error fxns for each of the output
        activation fxns (see Nabney pp.123-128,156-158 for more info).
        Borrowing heavily from the Netlab implementations for now.
        """
        y = self.fwd(x,w)
        if self.alg == 'leastsq':
            return sum(array(y-t)**2,axis=1)
        if self.outfxn == 'linear':
            # calculate & return SSE
            return 0.5*sum(sum(array(y-t)**2,axis=1))
        elif self.outfxn == 'logistic':
            # calculate & return x-entropy
            return -1.0*sum(sum(t*log2(y)+(1-t)*log2(1-y),axis=1))
        elif self.outfxn == 'softmax':
            # calculate & return entropy
            return -1.0*sum(sum(t*log2(y),axis=1))
        else:
            # this shouldn't happen...return SSE as a reasonable default
            return 0.5*sum(sum((y - t)**2,axis=1))
        
    def errgrad(self,w,x,t):
        """ Error gradient fxns for canonical error fxns (see above, and
        Nabney pp.127-128,156-158)
        ** Includes error-backpropagation (Netlab splits these fxns)
        Inputs:
            w   - weight vector (don't really know why I pass this around...)
            x   - input patterns
            t   - targets
        Outputs:
            g   - gradient
            
            
            ***N.B.*********************************************************
             I'M DOING SOMETHING WRONG HERE, EVIDENTLY, AS THE OPTIMIZATION 
             FXNS THAT DEPEND ON A GRADIENT FXN AREN'T DOING WHAT THEY'RE 
             SUPPOSED TO (i.e. they're not optimizing anything)
            ****************************************************************
        """
        # get output and hidden activation patterns for a full forward pass
        y,z = self.fwd(x,w,True)
        outdeltas = y-t
        
        # compute second-layer weight and bias gradients
        # THIS IS AN AWFUL-LOOKING HACK, BUT I HAVEN'T FOUND A BETTER
        # WAY TO DO IT, YET...
        w2grad = zeros((shape(x)[0],shape(z)[1]*shape(outdeltas)[1]),dtype=Float)
        for i in range(shape(w2grad)[0]):
            w2grad[i] = outer(outdeltas[i],z[i]).reshape(size(outdeltas[i])*size(z[i]))
        w2grad = sum(w2grad)
        b2grad = sum(outdeltas)
        # backpropagate...AGAIN WITH THE FUGLY HACK...PLUS I HAVE TO EXPLICITLY
        # LOOP OVER ALL INPUT PATTERNS....*bleah*...
        hiddeltas = zeros((shape(x)[0],self.nhid),dtype=Float)
        for i in range(shape(hiddeltas)[0]):
            for j in range(shape(hiddeltas)[1]):
                hiddeltas[i][j] = (1-z[i][j]**2)*sum(diag(outer(self.w2[j],outdeltas[i])))
        # compute first-layer weight and bias gradients
        w1grad = zeros((shape(x)[0],shape(x)[1]*shape(hiddeltas)[1]),dtype=Float)
        for i in range(shape(w1grad)[0]):
            w1grad[i] = outer(hiddeltas[i],x[i]).reshape(size(hiddeltas[i])*size(x[i]))
        w1grad = sum(w1grad)
        b1grad = sum(hiddeltas)
        # pack into a single vector and return it
        g = hstack([w1grad.reshape(size(w1grad)),
                    b1grad.reshape(size(b1grad)),
                    w2grad.reshape(size(w2grad)),
                    b2grad.reshape(size(b2grad))])
        return g
        
    def train(self,x,t):
        """ Train a multilayer perceptron with a user-specified algorithm.
        Inputs:
            x   - matrix of input data
            t   - matrix of target outputs
            alg - training algorithm, one of {simplex,bfgs,ncg,leastsq}
        Outputs:
            w   - post-optimization weight vector
        """
        # N.B. doing nothing with the specified algorithm for now,
        # just optimizing with the leastsq fxn in scipy.optimize
        if self.alg == 'simplex':
            from scipy.optimize import fmin
            w = fmin(self.errfxn,self.w_packed,args=(x,t),full_output=True)
        elif self.alg == 'bfgs':
            from scipy.optimize import fmin_bfgs
            # version of this that uses errgrad doesn't converge
            #w = fmin_bfgs(self.errfxn,self.w_packed,fprime=self.errgrad,args=(x,t),full_output=True)
            w = fmin_bfgs(self.errfxn,self.w_packed,args=(x,t),full_output=True)
        elif self.alg == 'cg':
            from scipy.optimize import fmin_cg
            #w = fmin_cg(self.errfxn,self.w_packed,self.errgrad,args=(x,t),full_output=True)
            w = fmin_cg(self.errfxn,self.w_packed,args=(x,t),full_output=True)
        elif self.alg == 'ncg':
            from scipy.optimize import fmin_ncg
            w = fmin_ncg(self.errfxn,self.w_packed,self.errgrad,args=(x,t),\
                         full_output=True)
        else:
            # leastsq, or undef'd algorithm, in which case use leastsq as
            # a reasonable default
            if self.alg != 'leastsq':
                import sys
                print "Undefined algorithm, using least-squares"
                sys.stdout.flush()
            from scipy.optimize import leastsq
            w = leastsq(self.errfxn,self.w_packed,args=(x,t),\
                        full_output=True)    
        return w

def main():
    import os,sys,copy
    from scipy.io import read_array, write_array
    """ Approx test of module, using the oilTrn/oilTst data files that are 
    distributed with Netlab. Prints a bunch of info about weight vector and 
    error measures before and after optimization.
    """
    opt = raw_input("\nEnter desired optimizer (simplex,bfgs,cg,ncg,leastsq): ")
    print "\nCreating 12-5-2 MLP with linear outputs"
    net = mlp(12,5,2,'linear',opt)
    w_init = copy.copy(net.w_packed)
    print "\nLoading training and test sets...",
    trn_input = read_array('data/oilTrn.dat',lines=(3,-1),columns=(0,(1,12)))
    trn_targs = read_array('data/oilTrn.dat',lines=(3,-1),columns=(12,-1))
    tst_input = read_array('data/oilTst.dat',lines=(3,-1),columns=(0,(1,12)))
    tst_targs = read_array('data/oilTst.dat',lines=(3,-1),columns=(12,-1))
    print "done."
    sys.stdout.flush()
    
    print "\nInitial error: ",net.errfxn(net.w_packed,tst_input,tst_targs)
    retval = net.train(trn_input,trn_targs)
    net.w_packed = retval[0]
    print "\nFinal error: ",net.errfxn(net.w_packed,tst_input,tst_targs)

if __name__ == '__main__':
    main()


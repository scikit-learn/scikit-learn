# mlp.py
# by: Fred Mailhot
# last mod: 2006-06-16

from scipy import *

class mlp:
    """Class to define, train and test a multilayer perceptron."""

    _type = 'mlp'
    _outfxns = ('linear','logistic','softmax')
    _algs = ('simplex','bfgs','ncg','leastsq')
    _options = zeros(18)

    def __init__(self,nin,nhid,nout,fxn,w=None):
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
        self.outfxn = fxn
        self.nin = nin
        self.nhid = nhid
        self.nout = nout
        if w is not None:
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
        self.w1 = reshape(self.w_packed[:self.nin*self.nhid],\
                          (self.nin,self.nhid))
        self.b1 = reshape(self.w_packed[(self.nin*self.nhid):\
                                        (self.nin*self.nhid)+self.nhid],\
                                        (1,self.nhid))
        self.w2 = reshape(self.w_packed[(self.nin*self.nhid)+self.nhid: \
                            (self.nin*self.nhid)+self.nhid+\
                            (self.nhid*self.nout)],(self.nhid,self.nout))
        self.b2 = reshape(self.w_packed[(self.nin*self.nhid)+self.nhid+\
                                        (self.nhid*self.nout):],(1,self.nout))

    def packwts(self):
        """ Compile weight matrices w1,b1,w2,b2 from net into a
        single vector, suitable for optimization routines.
        """
        self.w_packed = hstack([reshape(self.w1,(size(self.w1),)),
                                reshape(self.b1,(size(self.b1),)),
                                reshape(self.w2,(size(self.w2),)),
                                reshape(self.b2,(size(self.b2),))])

    def fwd(self,inputs,hid=False):
        """ Propagate values forward through the net.
        Inputs:
            inputs  - self.nin*1 vector of inputs
            hid     - boolean specifying whether or not to return hidden
                      unit activations, False by default
        """
        self.unpackwts()
        
        z = tanh(array(matrix(inputs)*matrix(self.w1) + \
                       matrix(ones((len(inputs),1)))*matrix(self.b1)))
        o = array(matrix(z)*matrix(self.w2) + \
                  matrix(ones((len(z),1)))*matrix(self.b2))
        
        if self.outfxn == 'linear':
            y = o
        elif self.outfxn == 'logistic':     # TODO: check for overflow here...
            y = 1/(1+exp(-o))
        elif self.outfxn == 'softmax':      # TODO: and here...
            tmp = exp(o)
            y = tmp/(sum(temp,1)*ones((1,self.nout)))
            
        if hid:
            return y,z
        else:
            return y

    def residuals(self,w,x,t):
        """ Calculate output error vectors for all input data points
        Inputs:
            w   - weight vector
            x   - vector of inputs
            t   - vector of target vals
        Outputs:
            err - error vector for input patterns
            (i.e. err[0] == error vector for input pattern x[0])

        *N.B.*  This dupes a lot (all?) of the code from fwd(), since I
            think it's necessary for the optimize.leastsq routine. I may 
            eliminate the fwd() fxn eventually and return the residuals 
            and outputs in a tuple.
            
        **N.B.B.** I think I may be able to get rid of this, with the new
        implementation of the other error functions...I can tweak errfxn
        to return the residuals, I think
        """
        if self.w_packed[0] != w[0]:      # make dependence on w explicit
            w_old = self.w_packed
            self.w_packed = w 
        self.unpackwts()
        z = tanh(array(matrix(x)*matrix(self.w1) + \
                       matrix(ones((len(x),1)))*matrix(self.b1)))
        o = array(matrix(z)*matrix(self.w2) + \
                  matrix(ones((len(z),1)))*matrix(self.b2))
        if self.outfxn == 'linear':
            y = o
        elif self.outfxn == 'logistic':     # TODO: check for overflow here...
            y = 1/(1+exp(-o))
        elif self.outfxn == 'softmax':      # TODO: and here...
            tmp = exp(o)
            y = tmp/(sum(temp,1)*ones((1,self.nout)))
        err = t - y
        # returning err w/ squared entries summed along axis1
        return sum(err**2,axis=1)
    
    def errfxn(self,w,x,t):
        """ Implementing 'canonical' error fxns for each of the output
        activation fxns (see Nabney pp.123-128,156-158 for more info).
        Borrowing heavily from the Netlab implementations for now.
        """
        from math import log
        y = self.fwd(x)
        if self.outfxn == 'linear':
            # calculate & return SSE
            return 0.5*sum(sum((y - t)**2,axis=1))
        elif self.outfxn == 'logistic':
            # calculate & return x-entropy
            return -1.0*sum(sum(t*math.log(y,2)+(1-t)*math.log(1-y,2),axis=1))
        elif self.outfxn == 'softmax':
            # calculate & return entropy
            return -1.0*sum(sum(t*math.log(y,2),axis=1))
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
        """
        # get output and hidden activation patterns for a full forward pass
        y,z = self.fwd(x,True)
        outdeltas = y-t
        # compute second-layer weight and bias gradients
        w2grad = z.transpose()*outdeltas
        b2grad = sum(outdeltas,axis=1)
        # backpropagate
        hiddeltas = outdeltas*self.w2
        hiddeltas = hiddeltas*(1-z**2)
        # compute first-layer weight and bias gradients
        w1grad = x.transpose()*hiddeltas
        b1grad = sum(hiddeltas,axis=1)
        # pack into a single vector and return it
        g = hstack([reshape(w1grad,(size(w1grad),)),
                                reshape(b1grad,(size(b1grad),)),
                                reshape(w2grad,(size(w2grad),)),
                                reshape(b2,(size(b2grad),))])
        return g
        
        
    def train(self,x,t,alg):
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
        if alg == 'simplex':
            from scipy.optimize import fmin
            w = fmin(self.errfxn,self.w_packed,args=(x,t),full_output=True)
        elif alg == 'bfgs':
            from scipy.optimize import fmin_bfgs
            w = fmin_bfgs(self.errfxn,self.w_packed,fprime=self.errgrad,\
                          args=(x,t),full_output=True)
        elif alg == 'ncg':
            from scipy.optimize import fmin_ncg
            w = fmin_ncg(self.errfxn,self.w_packed,self.errgrad,args=(x,t),\
                         full_output=True)
        else:
            # leastsq, or undef'd algorithm, in which case use leastsq as
            # a reasonable default
            from scipy.optimize import leastsq
            w = leastsq(self.residuals,self.w_packed,args=(x,t),\
                        full_output=True)    
        return w

    def test(self,x,t):
        # ################################# #
        # ### THIS NEEDS EXTENSIVE WORK ### #
        # ### (i.e. do I even need it?) ### #
        # ################################# #
        e = t - self.fwd(x)
        return array(sum(matrix(e)*matrix(e).T))

def main():
    import os,sys,copy
    """ Approx test of module, using the oilTrn/oilTst data files that are 
    distributed with Netlab (the Matlab ANN toolkit). Prints a bunch of info
    about weight vector and error measures before and after optimization.
    """
    # ###################################### #
    # !!! TODO !!!                           #
    # ###################################### #
    # Implement "canonical" error fxns for   #
    # each of the output activation fxns, as #
    # well as their respective gradient fxns.#
    #                                        #
    # See Nabney pp. 123-128, 156-159        #
    # ###################################### #
    from scipy.io import read_array, write_array
    print "Creating 12-5-2 MLP with linear outputs...",
    net = mlp(12,5,2,'linear')
    net.packwts()
    w_init = copy.copy(net.w_packed)
    print "Done."
    sys.stdout.flush()
    
    print "Loading training and test sets...",
    trn_input = read_array('data/oilTrn.dat',lines=(3,-1),columns=(0,(1,12)))
    trn_targs = read_array('data/oilTrn.dat',lines=(3,-1),columns=(12,-1))
    tst_input = read_array('data/oilTst.dat',lines=(3,-1),columns=(0,(1,12)))
    tst_targs = read_array('data/oilTst.dat',lines=(3,-1),columns=(12,-1))
    print "Done."
    sys.stdout.flush()
    
    err_vec_init = net.residuals(net.w_packed,tst_input,tst_targs)
    
    opt = None
    while(opt not in net._algs):
        opt = raw_input("Enter desired optimizer (simplex,bfgs,ncg,leastsq): ")

    print "Training net with "+opt+" optimizer...",
    retval = net.train(trn_input,trn_targs,opt)
    net.w_packed = retval[0]
    print "Done."
    sys.stdout.flush()
    
    print "Weight vectors pre- and post-optimization:"
    print "Pre:"
    print "\tAverage: "+str(mean(w_init))
    print "\tVariance: "+str(var(w_init))
    print "Post:"
    print "\tAverage: "+str(mean(net.w_packed))    
    print "\tVariance: "+str(var(net.w_packed))
    sys.stdout.flush()

    print "Error vectors pre- and post-optimization:"
    err_vec_new = net.residuals(net.w_packed,tst_input,tst_targs)
    print "Pre:"
    print "\tMSE: "+str(sum(err_vec_init)/size(err_vec_init))
    print "\tSSE: "+str(sum(matrix(err_vec_init)*matrix(err_vec_init).T))
    print "Post:"
    print "\tMSE: "+str(sum(err_vec_new)/size(err_vec_new))
    print "\tSSE: "+str(sum(matrix(err_vec_new)*matrix(err_vec_new).T))
    sys.stdout.flush()

if __name__ == '__main__':
    main()


# /usr/bin/python
# Last Change: Fri Jul 14 05:00 PM 2006 J

# Module to implement GaussianMixture class.

import numpy as N
from numpy.random import randn, rand
import numpy.linalg as lin
import densities

MAX_DEV = 1e-10

# Right now, two main usages of a Gaussian Model are possible
#   - init a Gaussian Model with meta-parameters, and trains it
#   - set-up a Gaussian Model to sample it, draw ellipsoides 
#   of confidences. In this case, we would like to init it with
#   known values of parameters.
#
#   For now, we have to init with meta-parameters, and set 
#   the parameters afterward. There should be a better way ?
class GmParamError:
    """Exception raised for errors in gmm params

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message    = message
    
    def __str__(self):
        return self.message

class GM:
    """Gaussian Mixture class. This is a simple container class
    to hold Gaussian Mixture parameters (weights, mean, etc...).
    It can also draw itself (confidence ellipses) and samples itself.

    Is initiated by giving dimension, number of components and 
    covariance mode"""

    # I am not sure it is useful to have a spherical mode...
    _cov_mod    = ['diag', 'full']

    def __init__(self, d, k, mode = 'diag'):
        """Init a Gaussian model of k components, each component being a 
        d multi-variate Gaussian, with covariance matrix of style mode"""
        if mode not in self._cov_mod:
            raise GmmParamError("mode %s not recognized" + str(mode))

        self.d      = d
        self.k      = k
        self.mode   = mode

        # Init to 0 all parameters, with the right dimensions.
        # Not sure this is useful in python from an efficiency POV ?
        self.w   = N.zeros(k, float)
        self.mu  = N.zeros((k, d), float)
        if mode == 'diag':
            self.va  = N.zeros((k, d), float)
        elif mode == 'full':
            self.va  = N.zeros((k * d, d), float)

        self.is_valid   = False

    def set_param(self, weights, mu, sigma):
        """Set parameters of the model"""
        k, d, mode  = check_gmm_param(weights, mu, sigma)
        if not k == self.k:
            raise GmmParamError("Number of given components is %d, expected %d" 
                    % (shape(k), shape(self.k)))
        if not d == self.d:
            raise GmmParamError("Dimension of the given model is %d, expected %d" 
                    % (shape(d), shape(self.d)))
        if not mode == self.mode:
            raise GmmParamError("Given covariance mode is %s, expected %d"
                    % (mode, self.mode))
        self.w  = weights
        self.mu = mu
        self.va = sigma

        self.is_valid   = True

    def sample(self, nframes):
        """ Sample nframes frames from the model """
        if not self.is_valid:
            raise GmmParamError("""Parameters of the model has not been 
                set yet, please set them using self.set_param()""")

        # State index (ie hidden var)
        S   = gen_rand_index(self.w, nframes)
        # standard gaussian
        X   = randn(nframes, self.d)        

        if self.mode == 'diag':
            X   = self.mu[S, :]  + X * N.sqrt(self.va[S,:])
        elif self.mode == 'full':
            # Faster:
            cho = N.zeros((self.k, self.va.shape[1], self.va.shape[1]), float)
            for i in range(self.k):
                # Using cholesky looks more stable than sqrtm; sqrtm is not
                # available in numpy anyway, only in scipy...
                cho[i]  = lin.cholesky(self.va[i*self.d:i*self.d+self.d,:])

            for s in range(self.k):
                tmpind      = N.where(S == s)[0]
                X[tmpind]   = N.dot(X[tmpind], cho[s].transpose()) + self.mu[s]
        else:
            raise GmmParamError('cov matrix mode not recognized, this is a bug !')

        return X

    def conf_ellipses(self, c = [0, 1], npoints = 100):
        """Returns a list of confidence ellipsoids describing the Gmm
        defined by mu and va. c is the dimension we are projecting
        the variances on a 2d space. For now, the confidence level
        is fixed to 0.39.
        
        Returns:
            -Xe:    a list of x coordinates for the ellipses
            -Ye:    a list of y coordinates for the ellipses

        Example:
            Suppose we have w, mu and va as parameters for a mixture, then:
            
            gm      = GM(d, k)
            gm.set_param(w, mu, va)
            X       = gm.sample(1000)
            Xe, Ye  = gm.conf_ellipsoids()
            pylab.plot(X[:,0], X[:, 1], '.')
            for k in len(w):
                pylab.plot(Xe[k], Ye[k], 'r')
                
            Will plot samples X draw from the mixture model, and
            plot the ellipses of equi-probability from the mean with
            fixed level of confidence 0.39.  """
        # TODO: adjustable level (to do in gauss_ell). 
        # For now, a level of 0.39 means that we draw
        # ellipses for 1 standard deviation. 
        Xe  = []
        Ye  = []   
        if self.mode == 'diag':
            for i in range(self.k):
                xe, ye  = densities.gauss_ell(self.mu[i,:], self.va[i,:], dim = c, 
                        npoints = npoints)
                Xe.append(xe)
                Ye.append(ye)
        elif self.mode == 'full':
            for i in range(self.k):
                xe, ye  = densities.gauss_ell(self.mu[i,:], 
                        self.va[i*self.d:i*self.d+self.d,:], dim = c, 
                        npoints = npoints)
                Xe.append(xe)
                Ye.append(ye)

        return Xe, Ye
    
    def gen_param(self, d, nc, varmode = 'diag', spread = 1):
        """Generate valid parameters for a gaussian mixture model.
        d is the dimension, nc the number of components, and varmode
        the mode for cov matrices.

        This is a class method.

        Returns: w, mu, va
        """
        w   = abs(randn(nc))
        w   = w / sum(w)

        mu  = spread * randn(nc, d)
        if varmode == 'diag':
            va  = abs(randn(nc, d))
        elif varmode == 'full':
            va  = randn(nc * d, d)
            for k in range(nc):
                va[k*d:k*d+d]   = N.dot( va[k*d:k*d+d], 
                    va[k*d:k*d+d].transpose())
        else:
            raise GmmParamError('cov matrix mode not recognized')

        return w, mu, va

    gen_param = classmethod(gen_param)


# Function to generate a random index: this is kept outside any class,
# as the function can be useful for other
def gen_rand_index(p, n):
    """Generate a N samples vector containing random index between 1 
    and length(p), each index i with probability p(i)"""
    # TODO Check args here
    
    # TODO: check each value of inverse distribution is
    # different
    invcdf  = N.cumsum(p)
    uni     = rand(n)
    index   = N.zeros(n)

    # This one should be a bit faster
    for k in range(len(p)-1, 0, -1):
        blop        = N.where(N.logical_and(invcdf[k-1] <= uni, 
                    uni < invcdf[k]))
        index[blop] = k
        
    return index

def check_gmm_param(w, mu, va):
    """Check that w, mu and va are valid parameters for
    a mixture of gaussian: w should sum to 1, there should
    be the same number of component in each param, the variances
    should be positive definite, etc... 
    
    Params:
        w   = vector or list of weigths of the mixture (K elements)
        mu  = matrix: K * d
        va  = list of variances (vector K * d or square matrices Kd * d)

    returns:
        K   = number of components
        d   = dimension
        mode    = 'diag' if diagonal covariance, 'full' of full matrices
    """
        
    # Check that w is valid
    if N.fabs(N.sum(w)  - 1) > MAX_DEV:
        raise GmmParamError('weight does not sum to 1')
    
    if not len(w.shape) == 1:
        raise GmmParamError('weight is not a vector')

    # Check that mean and va have the same number of components
    K           = len(w)

    if N.ndim(mu) < 2:
        msg = "mu should be a K,d matrix, and a row vector if only 1 comp"
        raise GmmParamError(msg)
    if N.ndim(va) < 2:
        msg = """va should be a K,d / K *d, d matrix, and a row vector if
        only 1 diag comp"""
        raise GmmParamError(msg)

    (Km, d)     = mu.shape
    (Ka, da)    = va.shape

    if not K == Km:
        msg = "not same number of component in mean and weights"
        raise GmmParamError(msg)

    if not d == da:
        msg = "not same number of dimensions in mean and variances"
        raise GmmParamError(msg)

    if Km == Ka:
        mode = 'diag'
    else:
        mode = 'full'
        if not Ka == Km*d:
            msg = "not same number of dimensions in mean and variances"
            raise GmmParamError(msg)
        
    return K, d, mode
        
if __name__ == '__main__':
    # Meta parameters:
    #   - k = number of components
    #   - d = dimension
    #   - mode : mode of covariance matrices
    d       = 5
    k       = 5
    mode    = 'full'
    nframes = 1e3

    # Build a model with random parameters
    w, mu, va   = GM.gen_param(d, k, mode, spread = 3)
    gm          = GM(d, k, mode)
    gm.set_param(w, mu, va)

    # Sample nframes frames  from the model
    X   = gm.sample(nframes)

    # Plot
    import pylab as P

    P.plot(X[:, 0], X[:, 1], '.', label = '_nolegend_')

    # Real confidence ellipses with level 0.39
    Xre, Yre  = gm.conf_ellipses()
    P.plot(Xre[0], Yre[0], 'g', label = 'true confidence ellipsoides')
    for i in range(1,k):
        P.plot(Xre[i], Yre[i], 'g', label = '_nolegend_')

    P.legend(loc = 0)
    P.show()

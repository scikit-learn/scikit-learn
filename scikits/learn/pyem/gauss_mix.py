# /usr/bin/python
# Last Change: Sat Jun 09 03:00 PM 2007 J

# Module to implement GaussianMixture class.

import numpy as N
from numpy.random import randn, rand
import numpy.linalg as lin
import densities
import misc

# Right now, two main usages of a Gaussian Model are possible
#   - init a Gaussian Model with meta-parameters, and trains it
#   - set-up a Gaussian Model to sample it, draw ellipsoides 
#   of confidences. In this case, we would like to init it with
#   known values of parameters. This can be done with the class method 
#   fromval

# TODO:
#   - change bounds methods of GM class instanciations so that it cannot 
#   be used as long as w, mu and va are not set
#   - We have to use scipy now for chisquare pdf, so there may be other
#   methods to be used, ie for implementing random index.
#   - there is no check on internal state of the GM, that is does w, mu and va values
#   make sense (eg singular values)
#   - plot1d is still very rhough. There should be a sensible way to 
#   modify the result plot (maybe returns a dic with global pdf, component pdf and
#   fill matplotlib handles). Should be coherent with plot
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
    """

    # I am not sure it is useful to have a spherical mode...
    _cov_mod    = ['diag', 'full']

    #===============================
    # Methods to construct a mixture
    #===============================
    def __init__(self, d, k, mode = 'diag'):
        """Init a Gaussian Mixture of k components, each component being a 
        d multi-variate Gaussian, with covariance matrix of style mode.
        
        If you want to build a Gaussian Mixture with knowns weights, means
        and variances, you can use GM.fromvalues method directly"""
        if mode not in self._cov_mod:
            raise GmParamError("mode %s not recognized" + str(mode))

        self.d      = d
        self.k      = k
        self.mode   = mode

        # Init to 0 all parameters, with the right dimensions.
        # Not sure this is useful in python from an efficiency POV ?
        self.w   = N.zeros(k)
        self.mu  = N.zeros((k, d))
        if mode == 'diag':
            self.va  = N.zeros((k, d))
        elif mode == 'full':
            self.va  = N.zeros((k * d, d))

        self.is_valid   = False

    def set_param(self, weights, mu, sigma):
        """Set parameters of the model. Args should
        be conformant with metparameters d and k given during
        initialisation"""
        k, d, mode  = check_gmm_param(weights, mu, sigma)
        if not k == self.k:
            raise GmParamError("Number of given components is %d, expected %d" 
                    % (k, self.k))
        if not d == self.d:
            raise GmParamError("Dimension of the given model is %d, expected %d" 
                    % (d, self.d))
        if not mode == self.mode and not d == 1:
            raise GmParamError("Given covariance mode is %s, expected %s"
                    % (mode, self.mode))
        self.w  = weights
        self.mu = mu
        self.va = sigma

        self.is_valid   = True

    @classmethod
    def fromvalues(cls, weights, mu, sigma):
        """This class method can be used to create a GM model
        directly from its parameters weights, mean and variance
        
        w, mu, va   = GM.gen_param(d, k)
        gm  = GM(d, k)
        gm.set_param(w, mu, va)

        and
        
        w, mu, va   = GM.gen_param(d, k)
        gm  = GM.fromvalue(w, mu, va)

        Are equivalent """
        k, d, mode  = check_gmm_param(weights, mu, sigma)
        res = cls(d, k, mode)
        res.set_param(weights, mu, sigma)
        return res
        
    #=====================================================
    # Fundamental facilities (sampling, confidence, etc..)
    #=====================================================
    def sample(self, nframes):
        """ Sample nframes frames from the model """
        if not self.is_valid:
            raise GmParamError("""Parameters of the model has not been 
                set yet, please set them using self.set_param()""")

        # State index (ie hidden var)
        S   = gen_rand_index(self.w, nframes)
        # standard gaussian
        X   = randn(nframes, self.d)        

        if self.mode == 'diag':
            X   = self.mu[S, :]  + X * N.sqrt(self.va[S,:])
        elif self.mode == 'full':
            # Faster:
            cho = N.zeros((self.k, self.va.shape[1], self.va.shape[1]))
            for i in range(self.k):
                # Using cholesky looks more stable than sqrtm; sqrtm is not
                # available in numpy anyway, only in scipy...
                cho[i]  = lin.cholesky(self.va[i*self.d:i*self.d+self.d,:])

            for s in range(self.k):
                tmpind      = N.where(S == s)[0]
                X[tmpind]   = N.dot(X[tmpind], cho[s].transpose()) + self.mu[s]
        else:
            raise GmParamError('cov matrix mode not recognized, this is a bug !')

        return X

    def conf_ellipses(self, dim = misc._DEF_VIS_DIM, npoints = misc._DEF_ELL_NP, \
        level = misc._DEF_LEVEL):
        """Returns a list of confidence ellipsoids describing the Gmm
        defined by mu and va. Check densities.gauss_ell for details

        Returns:
            -Xe:    a list of x coordinates for the ellipses (Xe[i] is
            the array containing x coordinates of the ith Gaussian)
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
        if not self.is_valid:
            raise GmParamError("""Parameters of the model has not been 
                set yet, please set them using self.set_param()""")

        Xe  = []
        Ye  = []   
        if self.mode == 'diag':
            for i in range(self.k):
                xe, ye  = densities.gauss_ell(self.mu[i,:], self.va[i,:], 
                        dim, npoints, level)
                Xe.append(xe)
                Ye.append(ye)
        elif self.mode == 'full':
            for i in range(self.k):
                xe, ye  = densities.gauss_ell(self.mu[i,:], 
                        self.va[i*self.d:i*self.d+self.d,:], 
                        dim, npoints, level)
                Xe.append(xe)
                Ye.append(ye)

        return Xe, Ye
    
    def check_state(self):
        """
        """
        if not self.is_valid:
            raise GmParamError("""Parameters of the model has not been 
                set yet, please set them using self.set_param()""")

        if self.mode == 'full':
            raise NotImplementedError, "not implemented for full mode yet"
        
        # # How to check w: if one component is negligeable, what shall
        # # we do ?
        # M   = N.max(self.w)
        # m   = N.min(self.w)

        # maxc    = m / M

        # Check condition number for cov matrix
        cond    = N.zeros(self.k)
        ava     = N.absolute(self.va)
        for c in range(self.k):
            cond[c] = N.amax(ava[c,:]) / N.amin(ava[c,:])

        print cond

    def gen_param(self, d, nc, varmode = 'diag', spread = 1):
        """Generate valid parameters for a gaussian mixture model.
        d is the dimension, nc the number of components, and varmode
        the mode for cov matrices.

        This is a class method.

        Returns: w, mu, va
        """
        w   = abs(randn(nc))
        w   = w / sum(w, 0)

        mu  = spread * randn(nc, d)
        if varmode == 'diag':
            va  = abs(randn(nc, d))
        elif varmode == 'full':
            va  = randn(nc * d, d)
            for k in range(nc):
                va[k*d:k*d+d]   = N.dot( va[k*d:k*d+d], 
                    va[k*d:k*d+d].transpose())
        else:
            raise GmParamError('cov matrix mode not recognized')

        return w, mu, va

    gen_param = classmethod(gen_param)

    #=======================
    # Regularization methods
    #=======================
    def _regularize(self):
        raise NotImplemented("No regularization")

    #=================
    # Plotting methods
    #=================
    def plot(self, dim = misc._DEF_VIS_DIM, npoints = misc._DEF_ELL_NP, 
            level = misc._DEF_LEVEL):
        """Plot the ellipsoides directly for the model
        
        Returns a list of lines, so that their style can be modified. By default,
        the style is red color, and nolegend for all of them.
        
        Does not work for 1d"""
        if not self.is_valid:
            raise GmParamError("""Parameters of the model has not been 
                set yet, please set them using self.set_param()""")

        assert self.d > 1
        k       = self.k
        Xe, Ye  = self.conf_ellipses(dim, npoints, level)
        try:
            import pylab as P
            return [P.plot(Xe[i], Ye[i], 'r', label='_nolegend_')[0] for i in range(k)]
            #for i in range(k):
            #    P.plot(Xe[i], Ye[i], 'r')
        except ImportError:
            raise GmParamError("matplotlib not found, cannot plot...")

    def plot1d(self, level = 0.5, fill = 0, gpdf = 0):
        """This function plots the pdfs of each component of the model. 
        If gpdf is 1, also plots the global pdf. If fill is 1, fill confidence
        areas using level argument as a level value
        
        Returns a dictionary h of plot handles so that their properties can
        be modified (eg color, label, etc...):
            - h['pdf'] is a list of lines, one line per component pdf
            - h['gpdf'] is the line for the global pdf
            - h['conf'] is a list of filling area
        """
        # This is not optimized at all, may be slow. Should not be
        # difficult to make much faster, but it is late, and I am lazy
        # XXX separete the computation from the plotting
        if not self.d == 1:
            raise GmParamError("the model is not one dimensional model")
        from scipy.stats import norm
        nrm     = norm(0, 1)
        pval    = N.sqrt(self.va[:,0]) * nrm.ppf((1+level)/2)

        # Compute reasonable min/max for the normal pdf: [-mc * std, mc * std]
        # gives the range we are taking in account for each gaussian
        mc  = 3
        std = N.sqrt(self.va[:,0])
        m   = N.amin(self.mu[:, 0] - mc * std)
        M   = N.amax(self.mu[:, 0] + mc * std)

        np  = 500
        x   = N.linspace(m, M, np)
        Yf  = N.zeros(np)
        Yt  = N.zeros(np)

        # Prepare the dic of plot handles to return
        ks  = ['pdf', 'conf', 'gpdf']
        hp  = dict((i,[]) for i in ks)
        try:
            import pylab as P
            for c in range(self.k):
                y   = self.w[c]/(N.sqrt(2*N.pi) * std[c]) * \
                        N.exp(-(x-self.mu[c][0])**2/(2*std[c]**2))
                Yt  += y
                h   = P.plot(x, y, 'r', label ='_nolegend_')
                hp['pdf'].extend(h)
                if fill:
                    #P.axvspan(-pval[c] + self.mu[c][0], pval[c] + self.mu[c][0], 
                    #        facecolor = 'b', alpha = 0.2)
                    id1 = -pval[c] + self.mu[c]
                    id2 = pval[c] + self.mu[c]
                    xc  = x[:, N.where(x>id1)[0]]
                    xc  = xc[:, N.where(xc<id2)[0]]
                    Yf  = self.w[c]/(N.sqrt(2*N.pi) * std[c]) * \
                            N.exp(-(xc-self.mu[c][0])**2/(2*std[c]**2))
                    xc  = N.concatenate(([xc[0]], xc, [xc[-1]]))
                    Yf  = N.concatenate(([0], Yf, [0]))
                    h   = P.fill(xc, Yf, 
                            facecolor = 'b', alpha = 0.1, label='_nolegend_')
                    hp['conf'].extend(h)
                    #P.fill([xc[0], xc[0], xc[-1], xc[-1]], 
                    #        [0, Yf[0], Yf[-1], 0], facecolor = 'b', alpha = 0.2)
            if gpdf:
                h           = P.plot(x, Yt, 'r:', label='_nolegend_')
                hp['gpdf']  = h
            return hp
        except ImportError:
            raise GmParamError("matplotlib not found, cannot plot...")

    def _get_component_pdf(self, x):
        """Returns a list of pdf, one for each component. Summing them gives
        the pdf of the mixture."""
        # XXX: have a public function to compute the pdf at given points
        # instead...
        std = N.sqrt(self.va[:,0])
        retval = N.empty((x.size, self.k))
        for c in range(self.k):
            retval[:, c] = self.w[c]/(N.sqrt(2*N.pi) * std[c]) * \
                    N.exp(-(x-self.mu[c][0])**2/(2*std[c]**2))

        return retval

    def density_on_grid(self, dim = misc._DEF_VIS_DIM, nx = 50, ny = 50,
            maxlevel = 0.95):
        """Do all the necessary computation for contour plot of mixture's density.
        
        Returns X, Y, Z and V as expected by mpl contour function."""

        # Ok, it is a bit gory. Basically, we want to compute the size of the
        # grid. We use conf_ellipse, which will return a couple of points for
        # each component, and we can find a grid size which then is just big
        # enough to contain all ellipses. This won't work well if two
        # ellipsoids are crossing each other a lot (because this assumes that
        # at a given point, one component is largely dominant for its
        # contribution to the pdf).

        # XXX: we need log pdf, not the pdf... this can save some computing
        Xe, Ye = self.conf_ellipses(level = maxlevel, dim = dim)
        ax = [N.min(Xe), N.max(Xe), N.min(Ye), N.max(Ye)]

        w = ax[1] - ax[0]
        h = ax[3] - ax[2]
        X, Y, den = self._densityctr(N.linspace(ax[0]-0.2*w, ax[1]+0.2*w, nx), \
                N.linspace(ax[2]-0.2*h, ax[3]+0.2*h, ny), dim = dim)
        lden = N.log(den)
        V = [-5, -3, -1, -0.5, ]
        V.extend(N.linspace(0, N.max(lden), 4).tolist())
        return X, Y, lden, N.array(V)

    def _densityctr(self, xrange, yrange, dim = misc._DEF_VIS_DIM):
        """Helper function to compute density contours on a grid."""
        gr = N.meshgrid(xrange, yrange)
        X = gr[0].flatten()
        Y = gr[1].flatten()
        xdata = N.concatenate((X[:, N.newaxis], Y[:, N.newaxis]), axis = 1)
        # XXX refactor computing pdf
        dmu = self.mu[:, dim]
        dva = self._get_va(dim)
        den = densities.multiple_gauss_den(xdata, dmu, dva) * self.w
        den = N.sum(den, 1)
        den = den.reshape(len(yrange), len(xrange))

        X = gr[0]
        Y = gr[1]
        return X, Y, den

    def _get_va(self, dim):
        """Returns variance limited do dimension in dim."""
        dim = N.array(dim)
        if dim.any() < 0 or dim.any() >= self.d:
            raise ValueError("dim elements should be between 0 and dimension"\
                    " of the mixture.")
        if self.mode == 'diag':
            return self.va[:, dim]
        elif self.mode == 'full':
            tidx = N.array([N.array(dim) + i * self.d for i in range(self.k)])
            tidx.flatten()
            return self.va[tidx, dim]
        else:
            raise ValueError("Unkown mode")
    # Syntactic sugar
    def __repr__(self):
        repr    = ""
        repr    += "Gaussian Mixture:\n"
        repr    += " -> %d dimensions\n" % self.d
        repr    += " -> %d components\n" % self.k
        repr    += " -> %s covariance \n" % self.mode
        if self.is_valid:
            repr    += "Has initial values"""
        else:
            repr    += "Has no initial values yet"""
        return repr

    def __str__(self):
        return self.__repr__()

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
    index   = N.zeros(n, dtype=int)

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
    if N.fabs(N.sum(w, 0)  - 1) > misc._MAX_DBL_DEV:
        raise GmParamError('weight does not sum to 1')
    
    if not len(w.shape) == 1:
        raise GmParamError('weight is not a vector')

    # Check that mean and va have the same number of components
    K           = len(w)

    if N.ndim(mu) < 2:
        msg = "mu should be a K,d matrix, and a row vector if only 1 comp"
        raise GmParamError(msg)
    if N.ndim(va) < 2:
        msg = """va should be a K,d / K *d, d matrix, and a row vector if
        only 1 diag comp"""
        raise GmParamError(msg)

    (Km, d)     = mu.shape
    (Ka, da)    = va.shape

    if not K == Km:
        msg = "not same number of component in mean and weights"
        raise GmParamError(msg)

    if not d == da:
        msg = "not same number of dimensions in mean and variances"
        raise GmParamError(msg)

    if Km == Ka:
        mode = 'diag'
    else:
        mode = 'full'
        if not Ka == Km*d:
            msg = "not same number of dimensions in mean and variances"
            raise GmParamError(msg)
        
    return K, d, mode
        
if __name__ == '__main__':
    # Meta parameters:
    #   - k = number of components
    #   - d = dimension
    #   - mode : mode of covariance matrices
    d       = 5
    k       = 4

    # Now, drawing a model
    mode    = 'full'
    nframes = 1e3

    # Build a model with random parameters
    w, mu, va   = GM.gen_param(d, k, mode, spread = 3)
    gm          = GM.fromvalues(w, mu, va)

    # Sample nframes frames  from the model
    X   = gm.sample(nframes)

    # Plot the data
    import pylab as P
    P.plot(X[:, 0], X[:, 1], '.', label = '_nolegend_')

    # Real confidence ellipses with confidence level 
    level       = 0.50
    h           = gm.plot(level=level)

    # set the first ellipse label, which will appear in the legend
    h[0].set_label('confidence ell at level ' + str(level))

    P.legend(loc = 0)
    P.show()

# /usr/bin/python
# Last Change: Mon Jul 02 05:00 PM 2007 J

"""Module implementing GM, a class which represents Gaussian mixtures.

GM instances can be used to create, sample mixtures. They also provide
different plotting facilities, such as isodensity contour for multi dimensional
models, ellipses of confidence."""

__docformat__ = 'restructuredtext'

import numpy as N
from numpy.random import randn, rand
import numpy.linalg as lin
import densities as D
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
#   - there is no check on internal state of the GM, that is does w, mu and va
#   values make sense (eg singular values) - plot1d is still very rhough. There
#   should be a sensible way to modify the result plot (maybe returns a dic
#   with global pdf, component pdf and fill matplotlib handles). Should be
#   coherent with plot
class GmParamError(Exception):
    """Exception raised for errors in gmm params

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        Exception.__init__(self)
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
        """Init a Gaussian Mixture.

        :Parameters:
            d : int
                dimension of the mixture.
            k : int
                number of component in the mixture.
            mode : string
                mode of covariance

        :Returns:
            an instance of GM.

        Note
        ----

        Only full and diag mode are supported for now.

        :SeeAlso:
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
        if d > 1:
            self.is1d = False
        else:
            self.is1d = True

    def set_param(self, weights, mu, sigma):
        """Set parameters of the model. 
        
        Args should be conformant with metparameters d and k given during
        initialisation.
        
        :Parameters:
            weights : ndarray
                weights of the mixture (k elements)
            mu : ndarray
                means of the mixture. One component's mean per row, k row for k
                components.
            sigma : ndarray
                variances of the mixture. For diagonal models, one row contains
                the diagonal elements of the covariance matrix. For full
                covariance, d rows for one variance.

        Examples
        --------
        Create a 3 component, 2 dimension mixture with full covariance matrices

        >>> w = numpy.array([0.2, 0.5, 0.3])
        >>> mu = numpy.array([[0., 0.], [1., 1.]])
        >>> va = numpy.array([[1., 0.], [0., 1.], [2., 0.5], [0.5, 1]])
        >>> gm = GM(2, 3, 'full')
        >>> gm.set_param(w, mu, va)

        :SeeAlso:
            If you know already the parameters when creating the model, you can
            simply use the method class GM.fromvalues."""
        #XXX: when fromvalues is called, parameters are called twice...
        k, d, mode  = check_gmm_param(weights, mu, sigma)
        if not k == self.k:
            raise GmParamError("Number of given components is %d, expected %d" 
                    % (k, self.k))
        if not d == self.d:
            raise GmParamError("Dimension of the given model is %d, "\
                "expected %d" % (d, self.d))
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
        
        :Parameters:
            weights : ndarray
                weights of the mixture (k elements)
            mu : ndarray
                means of the mixture. One component's mean per row, k row for k
                components.
            sigma : ndarray
                variances of the mixture. For diagonal models, one row contains
                the diagonal elements of the covariance matrix. For full
                covariance, d rows for one variance.

        :Returns:
            gm : GM
                an instance of GM.

        Examples
        --------

        >>> w, mu, va   = GM.gen_param(d, k)
        >>> gm  = GM(d, k)
        >>> gm.set_param(w, mu, va)

        and
        
        >>> w, mu, va   = GM.gen_param(d, k)
        >>> gm  = GM.fromvalue(w, mu, va)

        are strictly equivalent."""
        k, d, mode  = check_gmm_param(weights, mu, sigma)
        res = cls(d, k, mode)
        res.set_param(weights, mu, sigma)
        return res
        
    #=====================================================
    # Fundamental facilities (sampling, confidence, etc..)
    #=====================================================
    def sample(self, nframes):
        """ Sample nframes frames from the model.
        
        :Parameters:
            nframes : int
                number of samples to draw.
        
        :Returns:
            samples : ndarray
                samples in the format one sample per row (nframes, d)."""
        if not self.is_valid:
            raise GmParamError("""Parameters of the model has not been 
                set yet, please set them using self.set_param()""")

        # State index (ie hidden var)
        S   = gen_rand_index(self.w, nframes)
        # standard gaussian
        X   = randn(nframes, self.d)        

        if self.mode == 'diag':
            X   = self.mu[S, :]  + X * N.sqrt(self.va[S, :])
        elif self.mode == 'full':
            # Faster:
            cho = N.zeros((self.k, self.va.shape[1], self.va.shape[1]))
            for i in range(self.k):
                # Using cholesky looks more stable than sqrtm; sqrtm is not
                # available in numpy anyway, only in scipy...
                cho[i]  = lin.cholesky(self.va[i*self.d:i*self.d+self.d, :])

            for s in range(self.k):
                tmpind      = N.where(S == s)[0]
                X[tmpind]   = N.dot(X[tmpind], cho[s].transpose()) + self.mu[s]
        else:
            raise GmParamError("cov matrix mode not recognized, "\
                    "this is a bug !")

        return X

    def conf_ellipses(self, dim = misc.DEF_VIS_DIM, npoints = misc.DEF_ELL_NP, 
            level = misc.DEF_LEVEL):
        """Returns a list of confidence ellipsoids describing the Gmm
        defined by mu and va. Check densities.gauss_ell for details

        :Parameters:
            dim : sequence
                sequences of two integers which represent the dimensions where to
                project the ellipsoid.
            npoints : int
                number of points to generate for the ellipse.
            level : float
                level of confidence (between 0 and 1).

        :Returns:
            Xe : sequence
                a list of x coordinates for the ellipses (Xe[i] is the array
                containing x coordinates of the ith Gaussian)
            Ye : sequence
                a list of y coordinates for the ellipses.

        Examples
        --------
            Suppose we have w, mu and va as parameters for a mixture, then:
            
            >>> gm      = GM(d, k)
            >>> gm.set_param(w, mu, va)
            >>> X       = gm.sample(1000)
            >>> Xe, Ye  = gm.conf_ellipsoids()
            >>> pylab.plot(X[:,0], X[:, 1], '.')
            >>> for k in len(w):
            ...    pylab.plot(Xe[k], Ye[k], 'r')
                
            Will plot samples X draw from the mixture model, and
            plot the ellipses of equi-probability from the mean with
            default level of confidence."""
        if self.is1d:
            raise ValueError("This function does not make sense for 1d "
                "mixtures.")

        if not self.is_valid:
            raise GmParamError("""Parameters of the model has not been 
                set yet, please set them using self.set_param()""")

        Xe  = []
        Ye  = []   
        if self.mode == 'diag':
            for i in range(self.k):
                xe, ye  = D.gauss_ell(self.mu[i, :], self.va[i, :], 
                        dim, npoints, level)
                Xe.append(xe)
                Ye.append(ye)
        elif self.mode == 'full':
            for i in range(self.k):
                xe, ye  = D.gauss_ell(self.mu[i, :], 
                        self.va[i*self.d:i*self.d+self.d, :], 
                        dim, npoints, level)
                Xe.append(xe)
                Ye.append(ye)

        return Xe, Ye
    
    def check_state(self):
        """Returns true if the parameters of the model are valid. 

        For Gaussian mixtures, this means weights summing to 1, and variances
        to be positive definite.
        """
        if not self.is_valid:
            raise GmParamError("Parameters of the model has not been"\
                "set yet, please set them using self.set_param()")

        # Check condition number for cov matrix
        if self.mode == 'diag':
            tinfo = N.finfo(self.va.dtype)
            if N.any(self.va < tinfo.eps):
                raise GmParamError("variances are singular")
        elif self.mode == 'full':
            try:
                d = self.d
                for i in range(self.k):
                    N.linalg.cholesky(self.va[i*d:i*d+d, :])
            except N.linalg.LinAlgError:
                raise GmParamError("matrix %d is singular " % i)

        else:
            raise GmParamError("Unknown mode")

        return True

    @classmethod
    def gen_param(cls, d, nc, mode = 'diag', spread = 1):
        """Generate random, valid parameters for a gaussian mixture model.

        :Parameters:
            d : int
                the dimension
            nc : int
                the number of components
            mode : string
                covariance matrix mode ('full' or 'diag').

        :Returns:
            w : ndarray
                weights of the mixture
            mu : ndarray
                means of the mixture
            w : ndarray
                variances of the mixture

        Notes
        -----
        This is a class method.
        """
        w   = N.abs(randn(nc))
        w   = w / sum(w, 0)

        mu  = spread * N.sqrt(d) * randn(nc, d)
        if mode == 'diag':
            va  = N.abs(randn(nc, d))
        elif mode == 'full':
            # If A is invertible, A'A is positive definite
            va  = randn(nc * d, d)
            for k in range(nc):
                va[k*d:k*d+d]   = N.dot( va[k*d:k*d+d], 
                    va[k*d:k*d+d].transpose())
        else:
            raise GmParamError('cov matrix mode not recognized')

        return w, mu, va

    #gen_param = classmethod(gen_param)

    def pdf(self, x, log = False):
        """Computes the pdf of the model at given points.

        :Parameters:
            x : ndarray
                points where to estimate the pdf. One row for one
                multi-dimensional sample (eg to estimate the pdf at 100
                different points in 10 dimension, data's shape should be (100,
                20)).
            log : bool
                If true, returns the log pdf instead of the pdf.

        :Returns:
            y : ndarray
                the pdf at points x."""
        if log:
            return D.logsumexp(
                    D.multiple_gauss_den(x, self.mu, self.va, log = True)
                        + N.log(self.w))
        else:
            return N.sum(D.multiple_gauss_den(x, self.mu, self.va) * self.w, 1)

    #=================
    # Plotting methods
    #=================
    def plot(self, dim = misc.DEF_VIS_DIM, npoints = misc.DEF_ELL_NP, 
            level = misc.DEF_LEVEL):
        """Plot the ellipsoides directly for the model
        
        Returns a list of lines handle, so that their style can be modified. By
        default, the style is red color, and nolegend for all of them.
        
        :Parameters:
            dim : sequence
                sequence of two integers, the dimensions of interest.
            npoints : int
                Number of points to use for the ellipsoids.
            level : int
                level of confidence (to use with fill argument)
        
        :Returns:
            h : sequence
                Returns a list of lines handle so that their properties
                can be modified (eg color, label, etc...):

        Note
        ----
        Does not work for 1d. Requires matplotlib
        
        :SeeAlso:
            conf_ellipses is used to compute the ellipses. Use this if you want
            to plot with something else than matplotlib."""
        if self.is1d:
            raise ValueError("This function does not make sense for 1d "
                "mixtures.")

        if not self.is_valid:
            raise GmParamError("""Parameters of the model has not been 
                set yet, please set them using self.set_param()""")

        k       = self.k
        Xe, Ye  = self.conf_ellipses(dim, npoints, level)
        try:
            import pylab as P
            return [P.plot(Xe[i], Ye[i], 'r', label='_nolegend_')[0] for i in
                    range(k)]
            #for i in range(k):
            #    P.plot(Xe[i], Ye[i], 'r')
        except ImportError:
            raise GmParamError("matplotlib not found, cannot plot...")

    def plot1d(self, level = misc.DEF_LEVEL, fill = False, gpdf = False):
        """Plots the pdf of each component of the 1d mixture.
        
        :Parameters:
            level : int
                level of confidence (to use with fill argument)
            fill : bool
                if True, the area of the pdf corresponding to the given
                confidence intervales is filled.
            gpdf : bool
                if True, the global pdf is plot.
        
        :Returns:
            h : dict
                Returns a dictionary h of plot handles so that their properties
                can be modified (eg color, label, etc...):
                - h['pdf'] is a list of lines, one line per component pdf
                - h['gpdf'] is the line for the global pdf
                - h['conf'] is a list of filling area
        """
        if not self.is1d:
            raise ValueError("This function does not make sense for "
                "mixtures which are not unidimensional")

        # This is not optimized at all, may be slow. Should not be
        # difficult to make much faster, but it is late, and I am lazy
        # XXX separete the computation from the plotting
        if not self.d == 1:
            raise GmParamError("the model is not one dimensional model")
        from scipy.stats import norm
        nrm     = norm(0, 1)
        pval    = N.sqrt(self.va[:, 0]) * nrm.ppf((1+level)/2)

        # Compute reasonable min/max for the normal pdf: [-mc * std, mc * std]
        # gives the range we are taking in account for each gaussian
        mc  = 3
        std = N.sqrt(self.va[:, 0])
        m   = N.amin(self.mu[:, 0] - mc * std)
        M   = N.amax(self.mu[:, 0] + mc * std)

        np  = 500
        x   = N.linspace(m, M, np)
        Yf  = N.zeros(np)
        Yt  = N.zeros(np)

        # Prepare the dic of plot handles to return
        ks  = ['pdf', 'conf', 'gpdf']
        hp  = dict((i, []) for i in ks)
        try:
            import pylab as P
            for c in range(self.k):
                y   = self.w[c]/(N.sqrt(2*N.pi) * std[c]) * \
                        N.exp(-(x-self.mu[c][0])**2/(2*std[c]**2))
                Yt  += y
                h   = P.plot(x, y, 'r', label ='_nolegend_')
                hp['pdf'].extend(h)
                if fill:
                    #P.axvspan(-pval[c] + self.mu[c][0], pval[c] +
                    #self.mu[c][0], 
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
                    #        [0, Yf[0], Yf[-1], 0], facecolor = 'b', alpha =
                    #        0.2)
            if gpdf:
                h           = P.plot(x, Yt, 'r:', label='_nolegend_')
                hp['gpdf']  = h
            return hp
        except ImportError:
            raise GmParamError("matplotlib not found, cannot plot...")

    def density_on_grid(self, dim = misc.DEF_VIS_DIM, nx = 50, ny = 50,
            nl = 20, maxlevel = 0.95, V = None):
        """Do all the necessary computation for contour plot of mixture's
        density.
        
        :Parameters:
            dim : sequence
                sequence of two integers, the dimensions of interest.
            nx : int
                Number of points to use for the x axis of the grid
            ny : int
                Number of points to use for the y axis of the grid
            nl : int
                Number of contour to plot.
        
        :Returns:
            X : ndarray
                points of the x axis of the grid
            Y : ndarray
                points of the y axis of the grid
            Z : ndarray
                values of the density on X and Y
            V : ndarray
                Contour values to display.
            
        Note
        ----
        X, Y, Z and V are as expected by matplotlib contour function."""
        if self.is1d:
            raise ValueError("This function does not make sense for 1d "
                "mixtures.")

        # Ok, it is a bit gory. Basically, we want to compute the size of the
        # grid. We use conf_ellipse, which will return a couple of points for
        # each component, and we can find a grid size which then is just big
        # enough to contain all ellipses. This won't work well if two
        # ellipsoids are crossing each other a lot (because this assumes that
        # at a given point, one component is largely dominant for its
        # contribution to the pdf).

        Xe, Ye = self.conf_ellipses(level = maxlevel, dim = dim)
        ax = [N.min(Xe), N.max(Xe), N.min(Ye), N.max(Ye)]

        w = ax[1] - ax[0]
        h = ax[3] - ax[2]
        X, Y, lden = self._densityctr(N.linspace(ax[0]-0.2*w, ax[1]+0.2*w, nx), \
                N.linspace(ax[2]-0.2*h, ax[3]+0.2*h, ny), dim = dim)
        # XXX: how to find "good" values for level ?
        if V is None:
            V = N.linspace(-5, N.max(lden), nl)
        return X, Y, lden, N.array(V)

    def _densityctr(self, rangex, rangey, dim = misc.DEF_VIS_DIM):
        """Helper function to compute density contours on a grid."""
        gr = N.meshgrid(rangex, rangey)
        X = gr[0].flatten()
        Y = gr[1].flatten()
        xdata = N.concatenate((X[:, N.newaxis], Y[:, N.newaxis]), axis = 1)
        dmu = self.mu[:, dim]
        dva = self._get_va(dim)
        den = GM.fromvalues(self.w, dmu, dva).pdf(xdata, log = True)
        den = den.reshape(len(rangey), len(rangex))

        X = gr[0]
        Y = gr[1]
        return X, Y, den

    def _get_va(self, dim):
        """Returns variance limited do 2 dimension in tuple dim."""
        assert len(dim) == 2
        dim = N.array(dim)
        if dim.any() < 0 or dim.any() >= self.d:
            raise ValueError("dim elements should be between 0 and dimension"\
                    " of the mixture.")
        if self.mode == 'diag':
            return self.va[:, dim]
        elif self.mode == 'full':
            ld = dim.size
            vaselid = N.empty((ld * self.k, ld), N.int)
            for i in range(self.k):
                vaselid[ld*i] = dim[0] + i * self.d
                vaselid[ld*i+1] = dim[1] + i * self.d
            vadid = N.empty((ld * self.k, ld), N.int)
            for i in range(self.k):
                vadid[ld*i] = dim
                vadid[ld*i+1] = dim
            return self.va[vaselid, vadid]
        else:
            raise ValueError("Unkown mode")

    # Syntactic sugar
    def __repr__(self):
        msg = ""
        msg += "Gaussian Mixture:\n"
        msg += " -> %d dimensions\n" % self.d
        msg += " -> %d components\n" % self.k
        msg += " -> %s covariance \n" % self.mode
        if self.is_valid:
            msg += "Has initial values"""
        else:
            msg += "Has no initial values yet"""
        return msg

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
    a mixture of gaussian.
    
    w should sum to 1, there should be the same number of component in each
    param, the variances should be positive definite, etc... 
    
    :Parameters:
        w : ndarray
            vector or list of weigths of the mixture (K elements)
        mu : ndarray
            matrix: K * d
        va : ndarray
            list of variances (vector K * d or square matrices Kd * d)

    :Returns:
        k : int
            number of components
        d : int
            dimension
        mode : string
            'diag' if diagonal covariance, 'full' of full matrices
    """
        
    # Check that w is valid
    if not len(w.shape) == 1:
        raise GmParamError('weight should be a rank 1 array')

    if N.fabs(N.sum(w)  - 1) > misc._MAX_DBL_DEV:
        raise GmParamError('weight does not sum to 1')
    
    # Check that mean and va have the same number of components
    K = len(w)

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
    pass

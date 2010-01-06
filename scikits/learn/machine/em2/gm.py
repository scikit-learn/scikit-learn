# /usr/bin/python
# Last Change: Thu Jan 22 10:00 PM 2009 J

"""Module implementing GM, a class which represents Gaussian mixtures.

GM instances can be used to create, sample mixtures. They also provide
different plotting facilities, such as isodensity contour for multi dimensional
models, ellipses of confidence, etc..."""

import numpy as np
from numpy.random import randn, rand

# TODO:
#   - We have to use scipy now for chisquare pdf, so there may be other
#   methods to be used, ie for implementing random index.
#   - there is no check on internal state of the GM, that is does w, mu and va
#   values make sense (eg singular values) - plot1d is still very rhough. There
#   should be a sensible way to modify the result plot (maybe returns a dic
#   with global pdf, component pdf and fill matplotlib handles). Should be
#   coherent with plot
class GM:
    """Gaussian Mixture class. This is a simple container class to hold
    Gaussian Mixture parameters (weights, mean, etc...). It can estimate
    (log)-likelihood of mixtures, or one component. It can also be used to
    sample data from mixtures."""

    # I am not sure it is useful to have a spherical mode...
    _cov_mod    = ['diag', 'full']

    #===============================
    # Methods to construct a mixture
    #===============================
    def __init__(self, d, k, mode='diag'):
        """Create a Gaussian Mixture instance.

        Parameters
        ----------
        d : int
            dimension of the mixture.
        k : int
            number of component in the mixture.
        mode : {'diag', 'full', 'spherical'}
            mode of covariance

        Returns
        -------
        an instance of GM

        Notes
        -----
        Only full and diag mode are supported for now.

        SeeAlso
        -------
        If you want to build a Gaussian Mixture with knowns weights, means
        and variances, you can use GM.fromvalues method directly"""
        if mode not in self._cov_mod:
            raise ValueError("mode %s not recognized" + str(mode))

        self.d = d
        self.k = k
        self.mode = mode

        # Init to 0 all parameters, with the right dimensions.
        # Not sure this is useful in python from an efficiency POV ?
        self.w   = np.zeros(k)
        self.mu  = np.zeros((k, d))
        if mode == 'diag':
            self.va  = np.zeros((k, d))
        elif mode == 'full':
            self.va  = np.zeros((k * d, d))

        self.__is_valid   = False
        if d > 1:
            self.__is1d = False
        else:
            self.__is1d = True

    def setparams(self, w, mu, va):
        """Set parameters of the model.

        Args should be conformant with meta-parameters d and k given during
        initialisation.

        Parameters
        ----------
        w: ndarray
            weights of the mixture (k items)
        mu : ndarray
            means of the mixture. One component's mean per row, k row for k
            components.
        va : ndarray
            variances of the mixture. For diagonal models, one row contains
            the diagonal elements of the covariance matrix. For full
            covariance, d rows for one variance.

        Examples
        --------
        Create a 3 component, 2 dimension mixture with full covariance matrices

        >>> w = np.array([0.2, 0.5, 0.3])
        >>> mu = np.array([[0., 0.], [1., 1.]])
        >>> va = np.array([[1., 0.], [0., 1.], [2., 0.5], [0.5, 1]])
        >>> gm = GM(2, 3, 'full')
        >>> gm.setparams(w, mu, va)

        SeeAlso
        -------
        If you know already the parameters when creating the model, you can
        simply use the method class GM.fromvalues.
        """
        k, d, mode = _check_gmm_param(w, mu, va)
        if not k == self.k:
            raise ValueError("Number of given components is %d, expected %d"
                             % (k, self.k))
        if not d == self.d:
            raise ValueError("Dimension of the given model is %d, "\
                             "expected %d" % (d, self.d))
        if not mode == self.mode and not d == 1:
            raise ValueError("Given covariance mode is %s, expected %s"
                             % (mode, self.mode))
        self.w = w
        self.mu = mu
        self.va = va

        self.__is_valid   = True

    @classmethod
    def fromvalues(cls, w, mu, va):
        """This class method can be used to create a GM model directly from its
        parameters weights, mean and variance

        Parameters
        ----------
        w: ndarray
            weights of the mixture (k elements)
        mu : ndarray
            means of the mixture. One component's mean per row, k row for k
            components.
        va : ndarray
            variances of the mixture. For diagonal models, one row contains
            the diagonal elements of the covariance matrix. For full
            covariance, d rows for one variance.

        Returns
        -------
        gm : GM
            an instance of GM.

        Examples
        --------
        >>> w, mu, va = GM.genparams(d, k)
        >>> gm = GM(d, k)
        >>> gm.setparams(w, mu, va)

        and

        >>> w, mu, va = GM.genparams(d, k)
        >>> gm = GM.fromvalues(w, mu, va)

        are strictly equivalent."""
        k, d, mode  = _check_gmm_param(w, mu, va)
        res = cls(d, k, mode)
        res.setparams(w, mu, va)
        return res

    @classmethod
    def genparams(cls, d, nc, mode='diag', spread=1):
        """Generate random but valid parameters for a gaussian mixture model.

        Parameters
        ----------
        d : int
            the dimension
        nc : int
            the number of components
        mode : string
            covariance matrix mode ('full' or 'diag').

        Returns
        -------
        w : ndarray
            weights of the mixture
        mu : ndarray
            means of the mixture
        va : ndarray
            variances of the mixture

        Notes
        -----
        This is a class method.
        """
        w = np.abs(rand(nc))
        w /= np.sum(w)

        mu = spread * np.sqrt(d) * randn(nc, d)
        if mode == 'diag':
            va = np.abs(randn(nc, d))
        elif mode == 'full':
            # If A is invertible, A'A is positive definite
            va = randn(nc * d, d)
            for k in range(nc):
                va[k*d:k*d+d] = np.dot(va[k*d:k*d+d], va[k*d:k*d+d].T)
        else:
            raise ValueError('cov matrix mode not recognized')

        return w, mu, va

    def pdf(self, x, log=False):
        """Computes the pdf of the model at given points.

        Parameters
        ----------
        x: ndarray
            points where to estimate the pdf. One row for one multi-dimensional
            sample (eg to estimate the pdf at 100 different points in 10
            dimension, data's shape should be (100, 20)).
        log: bool
            If true, returns the log pdf instead of the pdf.

        Returns
        -------
        out: ndarray
            the pdf at points x."""
        if log:
            return logsumexp(mnormalik(x, self.mu, self.va, log=True)
                             + np.log(self.w))
        else:
            raise ValueError("Not implemented yet")

    def sample(self, nframes):
        """Sample nframes frames from the model.
        
        Parameters
        ----------
        nframes: int
            number of samples to draw.
        
        Returns
        -------
        samples: ndarray
            samples in the format one sample per row (nframes, d)."""
        if not self.__is_valid:
            raise ValueError("Parameters of the model has not been set yet" \
                             ", please set them using self.setparams()")

        # State index (ie hidden var)
        sti = randindex(self.w, nframes)
        # standard gaussian samples
        x = randn(nframes, self.d)        

        if self.mode == 'diag':
            x = self.mu[sti, :] + x * np.sqrt(self.va[sti, :])
        elif self.mode == 'full':
            # Faster:
            cho = np.zeros((self.k, self.va.shape[1], self.va.shape[1]))
            for i in range(self.k):
                # Using cholesky looks more stable than sqrtm; sqrtm is not
                # available in numpy anyway, only in scipy...
                cho[i] = np.linalg.cholesky(self.va[i*self.d:i*self.d+self.d, :])

            for s in range(self.k):
                tmpind = np.where(sti == s)[0]
                x[tmpind] = np.dot(x[tmpind], cho[s].T) + self.mu[s]
        else:
            raise ValueError("%s mode not recognized" % self.mode)

        return x

# Function to generate a random index: this is kept outside any class,
# as the function can be useful for other
def randindex(p, n):
    """Generate a N samples vector containing random index between 1
    and length(p), each index i with probability p(i)"""
    # TODO Check args here

    # TODO: check each value of inverse distribution is different
    invcdf = np.cumsum(p)
    uni = rand(n)
    index = np.zeros(n, dtype=int)

    # This one should be a bit faster
    # XXX: make this faster ?
    for k in range(len(p)-1, 0, -1):
        blop = np.where(np.logical_and(invcdf[k-1] <= uni, uni < invcdf[k]))
        index[blop] = k

    return index

def _check_gmm_param(w, mu, va):
    """Check that w, mu and va are valid parameters for a mixture of gaussian.

    w should sum to 1, there should be the same number of component in each
    param, the variances should be positive definite, etc...

    Raise a ValueError if arguments are not valid

    Parameters
    ----------
    w : ndarray
        vector or list of weigths of the mixture (K elements)
    mu : ndarray
        matrix: K * d
    va : ndarray
        list of variances (vector K * d or square matrices Kd * d)

    Returns
    -------
    k : int
        number of components
    d : int
        dimension
    mode : string
        Mode for covariances (diagonal, full, etc...)
    """

    # Check that w is valid
    if not len(w.shape) == 1:
        raise ValueError('weight should be a rank 1 array')

    if np.fabs(np.sum(w) - 1.) > 1e-12:
        raise ValueError('weight does not sum to 1')

    # Check that mean and va have the same number of components
    k = len(w)

    if np.ndim(mu) < 2:
        msg = "mu should be a K,d matrix, and a row vector if only 1 comp"
        raise ValueError(msg)
    if np.ndim(va) < 2:
        msg = """\
va should be a K,d / K *d, d matrix, and a row vector if only 1 diag comp"""
        raise ValueError(msg)

    (km, d) = mu.shape
    (ka, da) = va.shape

    if not k == km:
        msg = "not same number of component in mean and weights"
        raise ValueError(msg)

    if not d == da:
        msg = "not same number of dimensions in mean and variances"
        raise ValueError(msg)

    if km == ka:
        mode = 'diag'
    else:
        mode = 'full'
        if not ka == km*d:
            msg = "not same number of dimensions in mean and variances"
            raise ValueError(msg)

    return k, d, mode

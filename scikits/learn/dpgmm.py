import numpy as np

from . import mixture

from scipy.special import digamma

# Author: Alexandre Passos (alexandre.tp@gmail.com)
#
# Based on mixture.py by: 
#         Ron Weiss <ronweiss@gmail.com>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#

def norm(v):
    """The squared norm of vector v, as a scalar. As dot is not always
    scalar, using np.sum."""
    return np.sum(v*v)

def diagnorm(x, Sigma):
    """The x^T Sigma x norm when x is a vector representing a diagonal matrix."""
    return np.sum(x*Sigma*x)

def squarenorm(x, Sigma):
    """The x^T Sigma x norm when Sigma is a matrix."""
    return np.sum(np.dot(np.dot(x, Sigma), x))

def lognormalize(v):
    """Given a vector of unnormalized log-probabilites v returns a
 vector of normalized probabilities"""
    return np.exp(v-np.logaddexp.reduce(v))


class DPGMM(mixture.GMM):
    """Variational Inference for the Infinite Gaussian Mixture Model

    Stick-breaking Representation of a Gaussian mixture model
    probability distribution. This class allows for easy and efficient
    inference of an approximate posterior distribution over the
    parameters of a gaussian mixture model with a variable number of
    components (smaller than the truncation parameter n_states).

    Initialization is with normally-distributed means and identity
    covariance, for proper convergence.

    Parameters
    ----------
    n_states : int, optional
        Number of mixture components. Defaults to 1.

    cvtype : string (read-only), optional
        String describing the type of covariance parameters to
        use.  Must be one of 'spherical', 'tied', 'diag', 'full'.
        Defaults to 'diag'.


    Attributes
    ----------
    cvtype : string (read-only)
        String describing the type of covariance parameters used by
        the DP-GMM.  Must be one of 'spherical', 'tied', 'diag', 'full'.
    n_features : int
        Dimensionality of the Gaussians.
    n_states : int (read-only)
        Number of mixture components.
    weights : array, shape (`n_states`,)
        Mixing weights for each mixture component.
    means : array, shape (`n_states`, `n_features`)
        Mean parameters for each mixture component.
    covars : array
        Covariance parameters for each mixture component.  The shape
        depends on `cvtype`:
            (`n_states`,)                             if 'spherical',
            (`n_features`, `n_features`)              if 'tied',
            (`n_states`, `n_features`)                if 'diag',
            (`n_states`, `n_features`, `n_features`)  if 'full'
    converged_ : bool
        True when convergence was reached in fit(), False
        otherwise.

    Methods
    -------
    decode(X)
        Find most likely mixture components for each point in `X`.
    eval(X)
        Compute a lower-bound of the log likelihood of `X` under the model 
        and an approximate posterior distribution over mixture components.
    fit(X)
        Estimate the posterior of themodel parameters from `X` using the
        variational mean-field algorithm.
    predict(X)
        Like decode, find most likely mixtures components for each
        observation in `X`.
    rvs(n=1)
        Generate `n` samples from the posterior for the model.
    score(X)
        Compute the log likelihood of `X` under the model.

    Examples
    --------
    >>> import numpy as np
    >>> from scikits.learn import mixture
    >>> g = mixture.DPGMM(n_states=2)

    >>> # Generate random observations with two modes centered on 0
    >>> # and 10 to use for training.
    >>> np.random.seed(0)
    >>> obs = np.concatenate((np.random.randn(100, 1),
    ...                       10 + np.random.randn(300, 1)))
    >>> g.fit(obs)
    GMM(cvtype='diag', n_states=2)
    >>> g.weights
    array([ 0.25,  0.75])
    >>> g.means
    array([[ 0.05980802],
           [ 9.94199467]])
    >>> g.covars
    [array([[ 1.01682662]]), array([[ 0.96080513]])]
    >>> np.round(g.weights, 2)
    array([ 0.25,  0.75])
    >>> np.round(g.means, 2)
    array([[ 0.06],
           [ 9.94]])
    >>> np.round(g.covars, 2)
    ... #doctest: +NORMALIZE_WHITESPACE
    array([[[ 1.02]],
           [[ 0.96]]])
    >>> g.predict([[0], [2], [9], [10]])
    array([0, 0, 1, 1])
    >>> np.round(g.score([[0], [2], [9], [10]]), 2)
    array([-2.32, -4.16, -1.65, -1.19])

    >>> # Refit the model on new data (initial parameters remain the
    >>> # same), this time with an even split between the two modes.
    >>> g.fit(20 * [[0]] +  20 * [[10]])
    GMM(cvtype='diag', n_states=2)
    >>> np.round(g.weights, 2)
    array([ 0.5,  0.5])
    """

    def _bound_pxgivenz(self, x, k):
        bound = -0.5*self.n_features*np.log(2*np.pi)
        bound -= np.log(2*np.pi*np.e)
        if self.ctype == 'spherical':
            bound -= 0.5*self.n_features*(digamma(self._a[k])-np.log(self._b[k]))
            bound -= 0.5*(self._a[k]/self._b[k])*(norm(x-self._means[k])+D)
        elif self.ctype == 'diag':
            for i in xrange(self.n_features):
                bound -= 0.5*(digamma(self._a[k])-np.log(self._b[k]))
            bound -= 0.5*diagnorm(x-self._means[k], self._ab[k])
        elif self.ctype == 'tied':
            for i in xrange(self.n_features):
                bound -= 0.5*digamma((self._a-1-i)/2.) 
            bound -= 0.5*self.n_features*np.log(2)
            bound -= 0.5*np.log(self._detB)
            bound -= 0.5*squarenorm(x-self._means[k], self._ab)
        elif self.ctype == 'full':
            for i in xrange(self.n_features):
                bound -= 0.5*digamma((self._a[k]-1-i)/2.) 
            bound -= 0.5*self.n_features*np.log(2)
            bound -= 0.5*np.log(self._detB[k])
            bound -= 0.5*squarenorm(x-self._means[k], self._ab[k])
        else: 
            raise NotImplementedError("This ctype is not implemented: "+self.ctype)
        return bound
        
                                    
        

    def eval(self, obs):
        """Evaluate the model on data

        Compute the bound on log probability of `obs` under the model
        and return the posterior distribution (responsibilities) of
        each mixture component for each element of `obs`.

        This is done by computing the parameters for the mean-field of
        z for each observation.

        Parameters
        ----------
        obs : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.

        Returns
        -------
        logprob : array_like, shape (n_samples,)
            Log probabilities of each data point in `obs`
        posteriors: array_like, shape (n_samples, n_states)
            Posterior probabilities of each mixture component for each
            observation
        """
        obs = np.asanyarray(obs)
        z = np.zeros((obs.shape[0],self.n_states))
        p = np.zeros(self.n_states)
        bound = 0.
        for i in xrange(obs.shape[0]):
            for k in xrange(self.n_states):
                p[k] = z[i,k] = self._bound_pxgivenz(obs[i], k)
                z[i,k] += digamma(self._gamma[k,1]) 
                z[i,k] -= digamma(self._gamma[k,1]+self._gamma[k,2])
                for j in xrange(k):
                    z[i,k] += digamma(self._gamma[j,2])
                    z[i,k] -= digamma(self._gamma[j,1]+self._gamma[j,2])
                z[i] = lognormalize(z[i])
                bound += np.sum(z[i]*p)
        return bound, z


    def fit(self, X, n_iter=10, thresh=1e-2, params='wmc',
            init_params='wmc'):
        """Estimate model parameters with the variational
        algorithm.

        For a full derivation and description of the algorithm see 
        doc/dp-derivation/dp-derivation.tex

        A initialization step is performed before entering the em
        algorithm. If you want to avoid this step, set the keyword
        argument init_params to the empty string ''. Likewise, if you
        would like just to do an initialization, call this method with
        n_iter=0.

        Parameters
        ----------
        X : array_like, shape (n, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.
        n_iter : int, optional
            Number of EM iterations to perform.
        thresh : float, optional
            Convergence threshold.
        params : string, optional
            Controls which parameters are updated in the training
            process.  Can contain any combination of 'w' for weights,
            'm' for means, and 'c' for covars.  Defaults to 'wmc'.
        init_params : string, optional
            Controls which parameters are updated in the initialization
            process.  Can contain any combination of 'w' for weights,
            'm' for means, and 'c' for covars.  Defaults to 'wmc'.
        """

        ## initialization step

        X = np.asanyarray(X)

        if hasattr(self, 'n_features') and self.n_features != X.shape[1]:
            raise ValueError('Unexpected number of dimensions, got %s but '
                             'expected %s' % (X.shape[1], self.n_features))

        self.n_features = X.shape[1]

        if 'm' in init_params:
            self._means = np.random.normal(0, 1, size=(self.n_states, self.n_features))
        elif not hasattr(self, 'means'): 
            self._means = np.random.normal(0, 1, size=(self.n_states, self.n_features))

        if 'w' in init_params or not hasattr(self, 'weights'):
            self.weights = np.tile(1.0 / self._n_states, self._n_states)

        if 'c' in init_params:
            cv = np.identity(self.n_features)
            if not cv.shape:
                cv.shape = (1, 1)
            self._covars = _distribute_covar_matrix_to_match_cvtype(cv, self._cvtype, 
                                                                    self._n_states)
        elif not hasattr(self, 'covars'):
                self.covars = _distribute_covar_matrix_to_match_cvtype(
                    np.eye(self.n_features), self.cvtype, self.n_states)

        logprob = []
        # reset self.converged_ to False
        self.converged_ = False
        for i in xrange(n_iter):
            # Expectation step
            curr_logprob, posteriors = self.eval(X)
            logprob.append(curr_logprob.sum())

            # Check for convergence.
            if i > 0 and abs(logprob[-1] - logprob[-2]) < thresh:
                self.converged_ = True
                break

            # Maximization step
            self._do_mstep(X, posteriors, params, min_covar)

        return self

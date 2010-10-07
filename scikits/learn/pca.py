# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

import numpy as np
from scipy import linalg

from .base import BaseEstimator
from .utils.extmath import fast_logdet

def _assess_dimension_(spect, rk, n_samples, dim):
    """
    Compute the likelihood of a rank rk dataset
    embedded in gaussian noise of shape(n, dimf) having spectrum spect

    Parameters
    ----------
    spect: array of shape (n)
           data spectrum
    rk: int,  tested rank value
    n_samples: int, number of samples
    dim: int, embedding/emprical dimension

    Returns
    -------
    ll, float, The log-likelihood

    Notes
    -----
    This implements the method of Thomas P. Minka:
    Automatic Choice of Dimensionality for PCA. NIPS 2000: 598-604
    """
    if rk>dim:
        raise ValueError, "the dimension cannot exceed dim"
    from scipy.special import gammaln

    pu = -rk*np.log(2)
    for i in range(rk):
        pu += gammaln((dim - i)/2) - np.log(np.pi) * (dim - i)/2

    pl = np.sum(np.log(spect[:rk]))
    pl = -pl * n_samples / 2

    if rk == dim:
        pv = 0
        v = 1
    else:
        v = np.sum(spect[rk:dim]) / (dim - rk)
        pv = -np.log(v) * n_samples * (dim - rk)/2

    m = dim * rk - rk * (rk + 1) / 2
    pp = np.log(2 * np.pi) * (m + rk + 1)/2

    pa = 0
    spectrum_ = spect.copy()
    spectrum_[rk:dim] = v
    for i in range(rk):
        for j in range (i + 1, dim):
            pa += np.log((spect[i] - spect[j]) * \
                     (1./spectrum_[j] - 1./spectrum_[i])) + np.log(n_samples)

    ll = pu + pl + pv + pp -pa/2 - rk*np.log(n_samples)/2

    return ll

def _infer_dimension_(spect, n, p):
    """
    This method infers the dimension of a dataset of shape (n,p)
    with spectrum spect
    """
    ll = []
    for rk in range(min(n, p, len(spect))):
        ll.append(_assess_dimension_(spect, rk, n, p))
    ll = np.array(ll)
    return ll.argmax()


################################################################################
class PCA(BaseEstimator):
    """Principal component analysis (PCA)

    Parameters
    ----------
    X : array-like, shape = [n_samples, n_features]
        Training vector, where n_samples in the number of samples and
        n_features is the number of features.

    Attributes
    ----------
    n_comp : int, none or string,
             Number of components
             if k is not set all components are kept
             if k=='mle', Minka's mle is used to guess the dimension

    copy : bool
        If False, data passed to fit are overwritten

    components_ : array, [n_features, n_comp]
        Components with maximum variance

    explained_variance_ : array, [n_comp]
        Percentage of variance explained by each of the selected components.
        k is not set then all components are stored and the sum of
        explained variances is equal to 1.0

    Methods
    -------
    fit(X) : self
        Fit the model

    transform(X) : array
        Apply dimension reduction to k components


    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> from scikits.learn.pca import PCA
    >>> pca = PCA(n_comp=2)
    >>> pca.fit(X)
    PCA(n_comp=2, copy=True)
    >>> print pca.explained_variance_ratio_
    [ 0.99244289  0.00755711]

    See also
    --------

    """
    def __init__(self, n_comp=None, copy=True):
        self.n_comp = n_comp
        self.copy = copy

    def fit(self, X, **params):
        self._set_params(**params)
        X = np.atleast_2d(X)
        n_samples = X.shape[0]
        if self.copy:
            X = X.copy()
        # Center data
        self.mean_ = np.mean(X, axis=0)
        X -= self.mean_
        U, S, V = linalg.svd(X, full_matrices=False)
        self.explained_variance_ = (S**2)/n_samples
        self.explained_variance_ratio_ = self.explained_variance_ / \
                                        self.explained_variance_.sum()
        self.components_ = V.T
        if self.n_comp=='mle':
            self.n_comp = _infer_dimension_(self.explained_variance_,
                                            n_samples, X.shape[1])
            self.components_ = self.components_[:, :self.n_comp]
            self.explained_variance_ = self.explained_variance_[:self.n_comp]

        elif self.n_comp is not None:
            self.components_ = self.components_[:, :self.n_comp]
            self.explained_variance_ = self.explained_variance_[:self.n_comp]

        return self

    def transform(self, X):
        Xr = X - self.mean_
        Xr = np.dot(Xr, self.components_)
        return Xr


################################################################################
class ProbabilisticPCA(PCA):
    """ Additional layer on top of PCA that add a probabilistic evaluation
    """

    def fit(self, X, homoscedastic=True):
        """Additionally to PCA.fit, learns a covariance model

        Parameters
        ----------
        X: array of shape(n_samples, n_dim), test data
        homoscedastic: bool, optional,
                       if True, average variance across remaining dimensions
        """
        PCA.fit(self, X)
        self.dim = X.shape[1]
        Xr = X - self.mean_
        Xr -= np.dot(np.dot(Xr, self.components_), self.components_.T)
        n_samples = X.shape[0]
        if self.dim <= self.n_comp:
            delta = np.zeros(self.dim)
        elif homoscedastic:
            #delta = (Xr**2).sum()/(n_samples*(self.dim-self.n_comp)) *\
            #        np.ones(self.dim)
            delta = (Xr**2).sum() / (n_samples * self.dim) * np.ones(self.dim)
        else:
            delta = (Xr**2).mean(0) / (self.dim-self.n_comp)
        self.covariance_ = np.diag(delta)
        for k in range(self.n_comp):
            add_cov =  np.dot(
                self.components_[:, k:k+1], self.components_[:, k:k+1].T)
            self.covariance_ += self.explained_variance_[k] * add_cov
        return self

    def score(self, X):
        """Return a scoreassociated to new data

        Parameters
        ----------
        X: array of shape(n_samples, n_dim), test data

        Returns
        -------
        ll: array of shape (n_samples),
            log-likelihood of each row of X under the current model
        """
        Xr = X - self.mean_
        log_like = np.zeros(X.shape[0])
        self.precision_ = np.linalg.inv(self.covariance_)
        for i in range(X.shape[0]):
            log_like[i] = -.5 * np.dot(np.dot(self.precision_, Xr[i]), Xr[i])
        log_like += fast_logdet(self.precision_) - \
                                    self.dim / 2 * np.log(2 * np.pi)
        return log_like

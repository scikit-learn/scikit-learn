import numpy as np
from scipy.cluster.vq import kmeans2

from likelihoods import mnormalik, logsumexp
from gm import _check_gmm_param

def initkmeans(data, k):
    # XXX: This is bogus initialization should do better (in kmean with CV)
    (code, label) = kmeans2(data, data[:k], 5, minit='matrix')

    w = np.ones(k) / k
    mu = code.copy()
    va = np.zeros((k, d))
    for c in range(k):
        for i in range(d):
            va[c, i] = np.cov(data[np.where(label==c), i], rowvar=0)

    return w, mu, va

class Parameters:
    @classmethod
    def fromvalues(cls, w, mu, va):
        k, d, mode  = _check_gmm_param(w, mu, va)
        res = cls(d, k, mode)
        res.setparams(w, mu, va)
        return res

    def __init__(self, d, k, mode='diag'):
        self.mode = mode
        self.w = np.zeros(k)
        self.mu = np.zeros((k, d))
        if mode == 'diag':
            self.va = np.zeros((k, d))
        else:
            raise ValueError("Mode != diag not yet implemented")

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
        self.w[:] = w
        self.mu[:] = mu
        self.va[:] = va

    def update(self, ss):
        if not self.mode == ss.mode:
            raise ValueError("ss and parameters mode mismatch")

        if self.mode == 'diag':
            return self._diag_update(ss.w, ss.wx, ss.wxx)
        raise ValueError("Mode != diag not yet implemented")

    def _diag_update(self, w, wx, wxx):
        self.w[:] = w
        k = self.w.size

        for c in range(k):
            self.mu[c] = wx[c] / self.w[c]
            self.va[c] = wxx[c] / self.w[c] - self.mu[c] ** 2

class SStats:
    def __init__(self, d, k, mode='diag'):
        self.mode = mode
        self.k = k
        self.d = d
        self.w = np.zeros(k)
        self.wx = np.zeros((k, d))
        if mode == 'diag':
            self.wxx = np.zeros((k, d))
        else:
            raise ValueError("Mode != diag not yet implemented")
        self.np = 0.

    def reset(self):
        self.w[:] = 0
        self.wx[:] = 0
        self.wxx[:] = 0
        self.np = 0.

    def update(self, data, w, mu, va):
        if self.mode == 'diag':
            return self._diag_ss(data, w, mu, va)
        raise ValueError("Mode != diag not yet implemented")

    def compute(self, data, w, mu, va):
        self.reset()
        self.update(data, w, mu, va)

    def _diag_ss(self, data, w, mu, va):
        n = data.shape[0]
        # XXX: we should use clip here before exp, to remove too small values
        # (could be dealt with in logresp)
        nresp = np.exp(logresp(data, w, mu, va))
        if self.np == 0:
            self.w[:] = np.sum(nresp, axis=0) / n
            self.wx[:] = np.dot(nresp.T, data) / n
            self.wxx[:] = np.dot(nresp.T, data ** 2) / n
        else:
            self.w += np.sum(nresp, axis=0) / self.np
            self.w /= (1 + n / self.np)

            self.wx += np.dot(nresp.T, data) / self.np
            self.wx /= (1 + n / self.np)

            self.wxx += np.dot(nresp.T, data ** 2) / self.np
            self.wxx /= (1 + n / self.np)

        self.np += n

def logresp(data, w, mu, va):
    """Compute log responsabilities for a GMM, given data and parameters
        
    Note
    ----
    Computes the latent variable distribution (a posteriori probability)
    knowing the explicit data for the Gaussian model (w, mu, var): gamma(t,
    i) = P[state = i | observation = data(t); w, mu, va]
    """
    # compute the gaussian pdf
    resp = mnormalik(data, mu, va, log=True)
    # multiply by the weight
    resp += np.log(w)
    # Normalize to get a (log) pdf
    nresp = resp - logsumexp(resp)[:, np.newaxis]

    return nresp

class GMM:
    """A class to model a Gaussian Mixture Model (GMM)."""
    def __init__(self, gm):
        """Create a trainable mixture model.
        
        Initialize the model from a GM instance. This class implements all the
        necessary functionalities for EM.

        Parameters
        ----------
        gm: GM
            the mixture model to train.
        """
        self.gm = gm

    def logresp(self, data):
        """Compute log responsabilities.
        
        Note
        ----
        Computes the latent variable distribution (a posteriori probability)
        knowing the explicit data for the Gaussian model (w, mu, var): gamma(t,
        i) = P[state = i | observation = data(t); w, mu, va]
        """
        # compute the gaussian pdf
        resp = mnormalik(data, self.gm.mu, self.gm.va, log=True)
        # multiply by the weight
        resp += np.log(self.gm.w)
        # Normalize to get a (log) pdf
        nresp = resp - logsumexp(resp)[:, np.newaxis]

        return nresp

    def ss(self, data):
        if self.gm.mode == 'diag':
            return self._diag_ss(data)
        raise ValueError("Mode != diag not yet implemented")

    def _diag_ss(self, data):
        n = data.shape[0]
        nresp = np.exp(self.logresp(data))
        w = np.sum(nresp, axis=0) / n

        wx = np.dot(nresp.T, data) / n
        wxx = np.dot(nresp.T, data ** 2) / n

        return w, wx, wxx

    def update(self, w, wx, wxx):
        """Update the model parameters from SS."""
        if self.gm.mode == 'diag':
            return self._diag_update(w, wx, wxx)
        raise ValueError("Mode != diag not yet implemented")

    def _diag_update(self, w, wx, wxx):
        k, d = self.gm.k, self.gm.d

        mu = np.zeros((k, d))
        va = np.zeros((k, d))

        for c in range(k):
            mu[c] = wx[c] / w[c]
            va[c] = wxx[c] / w[c] - mu[c] ** 2

        return w, mu, va

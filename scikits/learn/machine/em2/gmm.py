import numpy as np

from likelihoods import mnormalik, logsumexp

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

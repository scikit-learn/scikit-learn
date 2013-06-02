import numpy as np
from scipy.optimize import fmin_l_bfgs_b

from sklearn.base import BaseEstimator


class SparseFiltering(BaseEstimator):

    def __init__(self, n_filters, maxfun=500, iprint=-1):
        self.n_filters = n_filters
        self.iprint = iprint
        self.maxfun = maxfun

    def fit(self, X, y=None, **params):
        """Fit the model with X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self.w_ = self._fit(X, **params)
        return self

    def fit_transform(self, X, y=None, **params):
        """Fit the model with X and apply the dimensionality reduction on X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)

        """
        self.w_ = self._fit(X, **params)
        return self._transform(X)

    def _fit(self, X):
        X = X.T  # transpose data in order to be consistent with the Matlab code
        X -= X.mean(0) # substract the mean from each image patch

        def objective_fct(w):
            # View 1d weight vector as a 2d matrix
            W = w.reshape(self.n_filters, X.shape[0])

            # Compute unnormalized features by multiplying weight matrix with
            # data
            F = W.dot(X)  # Linear Activation
            Fs = np.sqrt(F ** 2 + 1e-8)  # Soft-Absolute Activation

            # Normalize each feature to be equally active by dividing each
            # feature by its l2-norm across all examples
            L2Fs = np.apply_along_axis(np.linalg.norm, 1, Fs)
            NFs = Fs / L2Fs[:, None]
            # Normalize features per example, so that they lie on the unit
            # l2 -ball
            L2Fn = np.apply_along_axis(np.linalg.norm, 1, NFs.T)
            Fhat = NFs.T / L2Fn[:, None]
            # Compute sparsity of each feature over all example, i.e., compute
            # its l1-norm; the objective function is the sum over these
            # sparsities
            obj = np.apply_along_axis(np.linalg.norm, 1, Fhat, 1).sum()
            # Backprop through each feedforward step
            deltaW = l2grad(NFs.T, Fhat, L2Fn, np.ones_like(Fhat))
            deltaW = l2grad(Fs, NFs, L2Fs, deltaW.T)
            deltaW = (deltaW * (F / Fs)).dot(X.T)

            # Return objective value and gradient
            return obj, deltaW.flatten()

        def l2grad(X, Y, N, D):
            # Backpropagate through normalization
            return D / N[:, None] - Y * (D * X).sum(1)[:, None] / (N ** 2)[:, None]

        # Choose initial weights randomly
        w0 = np.random.random(X.shape[0] * self.n_filters) * 2 - 1
        # Use L-BFGS to find weights which correspond to a (local) minimum of
        # the objective function
        w, s, d = fmin_l_bfgs_b(objective_fct, w0, iprint=self.iprint,
                                maxfun=self.maxfun)

        return w.reshape(self.n_filters, X.shape[0],)

    def _transform(self, X):
        X = X.T  # transpose data in order to be consistent with the Matlab code
        X -= X.mean(0) # substract the mean from each image patch

        W = self.w_.reshape(self.n_filters, X.shape[0])

        # Compute unnormalized features by multiplying weight matrix with
        # data
        F = W.dot(X)  # Linear Activation
        Fs = np.sqrt(F ** 2 + 1e-8)  # Soft-Absolute Activation
        # Normalize each feature to be equally active by dividing each
        # feature by its l2-norm across all examples
        L2Fs = np.apply_along_axis(np.linalg.norm, 1, Fs)
        NFs = Fs / L2Fs[:, None]
        # Normalize features per example, so that they lie on the unit
        # l2 -ball
        L2Fn = np.apply_along_axis(np.linalg.norm, 1, NFs.T)
        Fhat = NFs.T / L2Fn[:, None]

        return Fhat.T

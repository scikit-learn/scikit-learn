"""Feature learning based on sparse filtering"""
# Author: Jan Hendrik Metzen
# License: BSD 3 clause

import numpy as np
from scipy.optimize import fmin_l_bfgs_b

from ..base import BaseEstimator


class SparseFiltering(BaseEstimator):
    """Sparse filtering

    Unsupervised learning of features using the sparse filtering algorithm.
    Features are linear in the inputs, i.e., f_j(x) = \sum_i w_{ij}x_i
    This algorithm does not try to model the
    data's distribution but rather to learn features which are sparsely
    activated in the sense of
        * Population Sparsity: for each image, only a small subset of features
                               is activated
        * Lifetime Sparsity: each feature is only activated on a small
                             subset of the examples
        * High Dispersal: Uniform activity distribution of features.

    This is encoded as an objective function which maps the the weight vector w
    onto a scalar value which is the smaller the more sparse the features are.
    L-BFGS is used to minimize this objective function.

    Parameters
    ----------
    n_features : int,
        Number of features to be learned.

    maxfun : int,
        Maximum number of evaluations of the objective function in L-BFGS-B.
        Defaults to 500.

    iprint : int,
        Verbosity of the L-BFGS-B. Prints information regarding the objective
        function every iprint iterations. Does not print any information if set
        to -1. Defaults to -1.

    Attributes
    ----------
    `w_` : array, [n_features, n_inputs]
        Sparse components extracted from the data.

    Notes
    -----
    This implements the method described `Jiquan Ngiam, Pang Wei Koh,
    Zhenghao Chen, Sonia Bhaskar, Andrew Y. Ng:
    Sparse Filtering. NIPS 2011: 1125-1133`
    and is based on the Matlab code provided in the supplementary material
    """

    def __init__(self, n_features, maxfun=500, iprint=-1):
        self.n_features = n_features
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
            W = w.reshape(self.n_features, X.shape[0])

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
        w0 = np.random.random(X.shape[0] * self.n_features) * 2 - 1
        # Use L-BFGS to find weights which correspond to a (local) minimum of
        # the objective function
        w, s, d = fmin_l_bfgs_b(objective_fct, w0, iprint=self.iprint,
                                maxfun=self.maxfun)

        return w.reshape(self.n_features, X.shape[0])

    def _transform(self, X):
        X = X.T  # transpose data in order to be consistent with the Matlab code
        X -= X.mean(0) # substract the mean from each image patch

        W = self.w_.reshape(self.n_features, X.shape[0])

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

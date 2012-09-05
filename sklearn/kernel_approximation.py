# -*- coding: utf-8 -*-
"""
The :mod:`sklearn.kernel_approximation` module implements several
approximate kernel feature maps base on Fourier transforms.
"""

# Author: Andreas Mueller <amueller@ais.uni-bonn.de>
#
# License: BSD Style.

import numpy as np
import scipy.sparse as sp

from .base import BaseEstimator
from .base import TransformerMixin
from .utils import array2d, atleast2d_or_csr, check_random_state
from .utils.extmath import safe_sparse_dot


class RBFSampler(BaseEstimator, TransformerMixin):
    """Approximates feature map of an RBF kernel by Monte Carlo approximation
    of its Fourier transform.

    Parameters
    ----------
    gamma: float
        parameter of RBF kernel: exp(-gamma * x**2)

    n_components: int
        number of Monte Carlo samples per original feature.
        Equals the dimensionality of the computed feature space.

    random_state : {int, RandomState}, optional
        If int, random_state is the seed used by the random number generator;
        if RandomState instance, random_state is the random number generator.

    Notes
    -----
    See "Random Features for Large-Scale Kernel Machines" by A. Rahimi and
    Benjamin Recht.
    """

    def __init__(self, gamma=1., n_components=100., random_state=None):
        self.gamma = gamma
        self.n_components = n_components
        self.random_state = random_state

    def fit(self, X, y=None):
        """Fit the model with X.

        Samples random projection according to n_features.

        Parameters
        ----------
        X: {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self : object
            Returns the transformer.
        """

        X = atleast2d_or_csr(X)
        self.random_state = check_random_state(self.random_state)
        n_features = X.shape[1]

        self.random_weights_ = (np.sqrt(self.gamma) *
                self.random_state.normal(size=(n_features, self.n_components)))
        self.random_offset_ = self.random_state.uniform(0,
                2 * np.pi, size=self.n_components)
        return self

    def transform(self, X, y=None):
        """Apply the approximate feature map to X.

        Parameters
        ----------
        X: {array-like, sparse matrix}, shape (n_samples, n_features)
            New data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new: array-like, shape (n_samples, n_components)
        """
        X = atleast2d_or_csr(X)
        projection = safe_sparse_dot(X, self.random_weights_)
        return (np.sqrt(2.) / np.sqrt(self.n_components)
                * np.cos(projection + self.random_offset_))


class SkewedChi2Sampler(BaseEstimator, TransformerMixin):
    """Approximates feature map of the "skewed chi-squared" kernel by Monte
    Carlo approximation of its Fourier transform.

    Parameters
    ----------
    skewedness: float
        "skewedness" parameter of the kernel. Needs to be cross-validated.

    n_components: int
        number of Monte Carlo samples per original feature.
        Equals the dimensionality of the computed feature space.

    random_state : {int, RandomState}, optional
        If int, random_state is the seed used by the random number generator;
        if RandomState instance, random_state is the random number generator.

    Notes
    -----
    See "Random Fourier Approximations for Skewed Multiplicative Histogram
    Kernels" by Fuxin Li, Catalin Ionescu and Cristian Sminchisescu.
    """

    def __init__(self, skewedness=1., n_components=100, random_state=None):
        self.skewedness = skewedness
        self.n_components = n_components
        self.random_state = random_state

    def fit(self, X, y=None):
        """Fit the model with X.

        Samples random projection according to n_features.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self : object
            Returns the transformer.
        """

        X = array2d(X)
        self.random_state = check_random_state(self.random_state)
        n_features = X.shape[1]
        uniform = self.random_state.uniform(size=(n_features,
            self.n_components))
        # transform by inverse CDF of sech
        self.random_weights_ = (1. / np.pi
                * np.log(np.tan(np.pi / 2. * uniform)))
        self.random_offset_ = self.random_state.uniform(0,
                2 * np.pi, size=self.n_components)
        return self

    def transform(self, X, y=None):
        """Apply the approximate feature map to X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            New data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new: array-like, shape (n_samples, n_components)
        """
        X = array2d(X)
        if (X < 0).any():
            raise ValueError("X may not contain entries smaller than zero.")

        projection = safe_sparse_dot(np.log(X + self.skewedness),
                self.random_weights_)

        return (np.sqrt(2.) / np.sqrt(self.n_components)
                * np.cos(projection + self.random_offset_))


class AdditiveChi2Sampler(BaseEstimator, TransformerMixin):
    """Approximate feature map for additive chi² kernel.

    Uses sampling the fourier transform of the kernel characteristic
    at regular intervals.

    Since the kernel that is to be approximated is additive, the components of
    the input vectors can be treated separately.  Each entry in the original
    space is transformed into 2×sample_steps+1 features, where sample_steps is
    a parameter of the method. Typical values of sample_steps include 1, 2 and
    3.

    Optimal choices for the sampling interval for certain data ranges can be
    computed (see the reference). The default values should be reasonable.

    Parameters
    ----------
    sample_steps: int, optional
        Gives the number of (complex) sampling points.
    sample_interval: float, optional
        Sampling interval. Must be specified when sample_steps not in {1,2,3}.

    Notes
    -----
    See `"Efficient additive kernels via explicit feature maps"
    <http://eprints.pascal-network.org/archive/00006964/01/vedaldi10.pdf>`_
    Vedaldi, A. and Zisserman, A., Computer Vision and Pattern Recognition 2010

    """

    def __init__(self, sample_steps=2, sample_interval=None):
        self.sample_steps = sample_steps
        self.sample_interval = sample_interval

    def fit(self, X, y=None):
        """Set parameters."""
        if self.sample_interval == None:
            # See reference, figure 2 c)
            if self.sample_steps == 1:
                self.sample_interval = 0.8
            elif self.sample_steps == 2:
                self.sample_interval = 0.5
            elif self.sample_steps == 3:
                self.sample_interval = 0.4
            else:
                raise ValueError("If sample_steps is not in [1, 2, 3],"
                    " you need to provide sample_interval")
        return self

    def transform(self, X, y=None):
        """Apply approximate feature map to X.

        Parameters
        ----------
        X: {array-like, sparse matrix}, shape = (n_samples, n_features)

        Returns
        -------
        X_new: {array, sparse matrix}, \
               shape = (n_samples, n_features × (2×sample_steps + 1))
            Whether the return value is an array of sparse matrix depends on
            the type of the input X.
        """

        X = atleast2d_or_csr(X)
        sparse = sp.issparse(X)

        # check if X has negative values. Doesn't play well with np.log.
        if ((X.data if sparse else X) < 0).any():
            raise ValueError("Entries of X must be non-negative.")
        # zeroth component
        # 1/cosh = sech
        # cosh(0) = 1.0

        transf = self._transform_sparse if sparse else self._transform_dense
        return transf(X)

    def _transform_dense(self, X):
        non_zero = (X != 0.0)
        X_nz = X[non_zero]

        X_step = np.zeros_like(X)
        X_step[non_zero] = np.sqrt(X_nz * self.sample_interval)

        X_new = [X_step]

        log_step_nz = self.sample_interval * np.log(X_nz)
        step_nz = 2 * X_nz * self.sample_interval

        for j in xrange(1, self.sample_steps):
            factor_nz = np.sqrt(step_nz /
                             np.cosh(np.pi * j * self.sample_interval))

            X_step = np.zeros_like(X)
            X_step[non_zero] = factor_nz * np.cos(j * log_step_nz)
            X_new.append(X_step)

            X_step = np.zeros_like(X)
            X_step[non_zero] = factor_nz * np.sin(j * log_step_nz)
            X_new.append(X_step)

        return np.hstack(X_new)

    def _transform_sparse(self, X):
        indices = X.indices.copy()
        indptr = X.indptr.copy()

        data_step = np.sqrt(X.data * self.sample_interval)
        X_step = sp.csr_matrix((data_step, indices, indptr),
                               shape=X.shape, dtype=X.dtype, copy=False)
        X_new = [X_step]

        log_step_nz = self.sample_interval * np.log(X.data)
        step_nz = 2 * X.data * self.sample_interval

        for j in xrange(1, self.sample_steps):
            factor_nz = np.sqrt(step_nz /
                             np.cosh(np.pi * j * self.sample_interval))

            data_step = factor_nz * np.cos(j * log_step_nz)
            X_step = sp.csr_matrix((data_step, indices, indptr),
                                   shape=X.shape, dtype=X.dtype, copy=False)
            X_new.append(X_step)

            data_step = factor_nz * np.sin(j * log_step_nz)
            X_step = sp.csr_matrix((data_step, indices, indptr),
                                   shape=X.shape, dtype=X.dtype, copy=False)
            X_new.append(X_step)

        return sp.hstack(X_new)

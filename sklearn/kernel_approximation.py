"""
The :mod:`sklearn.kernel_approximation` module implements several
approximate kernel feature maps base on Fourier transforms.
"""

# Author: Andreas Mueller <amueller@ais.uni-bonn.de>
#
# License: BSD Style.

import numpy as np
from .base import BaseEstimator
from .base import TransformerMixin
from .utils import check_random_state


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

    Notes
    -----
    For details see "Random Features for Large-Scale Kernel Machines" by A,
    Rahimi and Benjamin Recht for details."""

    def __init__(self, gamma=1., n_components=100., random_state=None):
        self.gamma = gamma
        self.n_components = n_components
        self.random_state = random_state

    def fit(self, X, y=None):
        """Fit the model with X.
        Samples random projection according to n_features

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.
        Returns
        -------
        self : object
            Returns the instance itself
        """

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
        X: array-like, shape (n_samples, n_features)
            New data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new array-like, shape (n_samples, n_components)
        """
        projection = np.dot(X, self.random_weights_)
        return (np.sqrt(2.) / np.sqrt(self.n_components)
                * np.cos(projection + self.random_offset_))


class SkewedChi2Sampler(BaseEstimator, TransformerMixin):
    """Approximates feature map of the "skewed chi-squared" kernel by Monte
    Carlo approximation of its Fourier transform.

    Parameters
    ----------
    c: float
        "skewedness" parameter of the kernel. Needs to be cross-validated.
    n_components: int
        number of Monte Carlo samples per original feature.
        Equals the dimensionality of the computed feature space.

    Notes
    -----
    See "Random Fourier Approximations for Skewed Multiplicative
    Histogram Kernels" by Fuxin Li, Catalin Ionescu and Cristian
    Sminchisescu for details.
    """

    def __init__(self, skewedness=1., n_components=100, random_state=None):
        self.skewedness = skewedness
        self.n_components = n_components
        self.random_state = random_state

    def fit(self, X, y=None):
        """Fit the model with X.
        Samples random projection according to n_features

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.
        Returns
        -------
        self : object
            Returns the instance itself
        """

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
        X_new array-like, shape (n_samples, n_components)
        """
        if (X < 0).any():
            raise ValueError("X may not contain entries smaller than zero.")

        projection = np.dot(np.log(X + self.skewedness), self.random_weights_)
        return (np.sqrt(2.) / np.sqrt(self.n_components)
                * np.cos(projection + self.random_offset_))


class AdditiveChi2Sampler(BaseEstimator, TransformerMixin):
    """Approximate feature map for additive chi^2 kernel.
    uses sampling the fourier transform of the kernel characteristic
    at regular intervals L.

    Since the kernel that is to be approximated is additive, the
    components of the input vectors can be treated separately.
    Each entry in the original space is transformed into 2n+1
    features, where n is a parameter of the method.
    Usually, n is 1, 2 or 3.

    Optimal choices for the sampling interval L for certain
    data ranges can be computed (see the reference).
    The default values should be reasonable.

    Parameters
    ----------
    n: int,     one of 1, 2 or 3.
                Gives the number of (complex) sampling points.
    L: float,   sampling interval

    Notes
    -----
    For details on the algorithm see `"Efficient additive kernels via explicit
    feature maps"
    <http://eprints.pascal-network.org/archive/00006964/01/vedaldi10.pdf>`_
    Vedaldi, A. and Zisserman, A.
    - Computer Vision and Pattern Recognition 2010
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
                raise ValueError("If n is not in [1, 2, 3], you need"
                    "to provide L")
        return self

    def transform(self, X, y=None):
        """Apply approximate feature map to X.

        Parameters
        ----------
        X array-like, shape (n_samples, n_features)

        Returns
        -------
        X_new array-like, shape (n_samples, n_features * (2n + 1))
        """

        # check if X has zeros. Doesn't play well with np.log.
        if (X <= 0).any():
            raise ValueError("Entries of X must be strictly positive.")
        X_new = []
        # zeroth component
        # 1/cosh = sech
        X_new.append(np.sqrt(X * self.sample_interval / np.cosh(0)))

        log_step = self.sample_interval * np.log(X)
        step = 2 * X * self.sample_interval

        for j in xrange(1, self.sample_steps):
            factor = np.sqrt(step / np.cosh(np.pi * j * self.sample_interval))
            X_new.append(factor * np.cos(j * log_step))
            X_new.append(factor * np.sin(j * log_step))
        return np.hstack(X_new)

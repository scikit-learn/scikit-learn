import numpy as np
#import cuv_python as cp
from sklearn.base import BaseEstimator


class FourierSampler(BaseEstimator):
    """ Base class for Fourier approximations to kernel functions"""
    def fit_transform(self, X, y=None):
        """Fit the model with X and apply the approximated feature map to X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new array-like, shape (n_samples, D)
        """
        return self.fit(X).transform(X)


class RBFSampler(FourierSampler):
    """Approximates feature map of an RBF kernel by Monte Carlo approximation
    of its Fourier transform.

    Parameters
    ----------
    gamma: float
        parameter of RBF kernel: exp(-gamma * x**2)
    D: int
        number of Monte Carlo samples per original feature.
        Equals the dimensionality of the computed feature space.

    See "Random Features for Large-Scale Kernel Machines" by A, Rahimi and
    Benjamin Recht for details."""

    def __init__(self, gamma=1., D=1.):
        self.gamma = gamma
        self.D = D

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

        n_features = X.shape[1]
        self.omega_ = (np.sqrt(self.gamma)
                * np.random.normal(size=(n_features, self.D)))
        self.b = np.random.uniform(0, 2 * np.pi, size=self.D)
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
        X_new array-like, shape (n_samples, D)
        """
        projection = np.dot(X, self.omega_)
        return np.sqrt(2.) / np.sqrt(self.D) * np.cos(projection + self.b)


class SkewedChi2Sampler(FourierSampler):
    """Approximates feature map of the "skewed chi-squared" kernel by Monte
    Carlo approximation of its Fourier transform.
    The "skewed chi-squared" kernel applies only to non-negative data and is
    particularly suited for histograms. It is designed to share certain
    beneficial properties of the chi-squared kernel.

    Parameters
    ----------
    c: float
        "skewedness" parameter of the kernel. Needs to be cross-validated.
    D: int
        number of Monte Carlo samples per original feature.
        Equals the dimensionality of the computed feature space.

    See "Random Fourier Approximations for Skewed Multiplicative
    Histogram Kernels" by Fuxin Li, Catalin Ionescu and Cristian
    Sminchisescu for details.
    """

    def __init__(self, c, D):
        self.c = c
        self.D = D

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
        n_features = X.shape[1]
        uniform = np.random.uniform(size=(n_features, self.D))
        # transform by inverse CDF of sech
        self.omega_ = 1. / np.pi * np.log(np.tan(np.pi / 2. * uniform))
        self.b = np.random.uniform(0, 2 * np.pi, size=self.D)
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
        X_new array-like, shape (n_samples, D)
        """
        projection = np.dot(np.log(X + self.c), self.omega_)
        return np.sqrt(2.) / np.sqrt(self.D) * np.cos(projection + self.b)

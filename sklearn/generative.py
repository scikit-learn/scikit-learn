"""
Bayesian Generative Classification
==================================
This module contains routines for general Bayesian generative classification.
Perhaps the best-known instance of generative classification is the Naive
Bayes Classifier, in which the distribution of each training class is
approximated by an axis-aligned multi-dimensional normal distribution, and
unknown points are evaluated by comparing their posterior probability under
each model.

This idea can be straightforwardly extended by using a more sophisticated
distribution model for each class: i.e. instead of simply modeling it as a
normal distribution, you might use a mixture of gaussians or a kernel density
estimate.

Mathematical Background
-----------------------
Bayesian Generative classification relies on Bayes' formula,

.. math::

    P(C|D) = \frac{P(D|C)P(C)}{P(D)}

where here, :math:`D` refers to the observed data of an unknown sample, and
:math:`C` refers to its classification.  The :math:`P(D|C)` term represents
the per-class likelihood (e.g. using the normal approximation to the given
class in the case of naive Bayes), and :math:`P(C)` gives the class prior,
often inferred from the training sample.  :math:`P(D)` acts as a normalization
parameter.  The final classification is the class :math:`C` which gives the
largest posterior probability :math:`P(D|C)`.
"""
# Author: Jake Vanderplas <jakevdp@cs.washington.edu>

__all__ = ['GenerativeBayes']
 
import numpy as np
from .neighbors import KernelDensity
from .mixture import GMM
from .base import BaseEstimator, clone
from .utils import array2d, check_random_state
from .utils.extmath import logsumexp
from .naive_bayes import BaseNB
 
 
class _NormalApproximation(BaseEstimator):
    """Normal Approximation Density Estimator"""
    def __init__(self):
        pass
 
    def fit(self, X):
        """Fit the Normal Approximation to data
 
        Parameters
        ----------
        X: array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.
        """
        X = array2d(X)
        epsilon = 1e-9
        self.mean = X.mean(0)
        self.var = X.var(0) + epsilon
        return self
 
    def score_samples(self, X):
        """Evaluate the model on the data
 
        Parameters
        ----------
        X : array_like
            An array of points to query.  Last dimension should match dimension
            of training data (n_features)
 
        Returns
        -------
        density : ndarray
            The array of density evaluations.  This has shape X.shape[:-1]
        """
        X = array2d(X)
        if X.shape[-1] != self.mean.shape[0]:
            raise ValueError("dimension of X must match that of training data")
        norm = 1. / np.sqrt(2 ** X.shape[-1] * np.sum(self.var))
        res = np.log(norm * np.exp(-0.5 * ((X - self.mean) ** 2
                                                 / self.var).sum(1)))
        return res
 
    def score(self, X):
        """Compute the log probability under the model.
 
        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.
 
        Returns
        -------
        logprob : array_like, shape (n_samples,)
            Log probabilities of each data point in X
        """
        return np.sum(np.log(self.eval(X)))
 
    def sample(self, n_samples=1, random_state=None):
        """Generate random samples from the model.
 
        Parameters
        ----------
        n_samples : int, optional
            Number of samples to generate. Defaults to 1.
 
        random_state: RandomState or an int seed (0 by default)
            A random number generator instance
 
        Returns
        -------
        X : array_like, shape (n_samples, n_features)
            List of samples
        """
        rng = check_random_state(random_state)
        return rng.normal(self.mean, np.sqrt(self.var),
                          size=(n_samples, len(self.mean)))        
 
 
MODEL_TYPES = {'norm_approx': _NormalApproximation,
               'gmm': GMM,
               'kde': KernelDensity}
 
 
class GenerativeBayes(BaseNB):
    """
    Generative Bayes Classifier

    This is a meta-estimator which performs generative Bayesian classification.

    Parameters
    ----------
    density_estimator : str, class, or instance
        The density estimator to use for each class.  Options are
            'norm_approx' : Normal Approximation
            'gmm' : Gaussian Mixture Model
            'kde' : Kernel Density Estimate
        Alternatively, a class or class instance can be specified.  The
        instantiated class should be a sklearn estimator, and contain a
        ``score_samples`` method with semantics similar to that in
        :class:`sklearn.neighbors.KDE` or :class:`sklearn.mixture.GMM`.
    **kwargs :
        additional keyword arguments to be passed to the constructor
        specified by density_estimator.
    """
    def __init__(self, density_estimator, **kwargs):
        if isinstance(density_estimator, str):
            dclass = MODEL_TYPES.get(density_estimator)
            self.density_estimator = dclass(**kwargs)
        elif isinstance(density_estimator, type):
            self.density_estimator = density_estimator(**kwargs)
        else:
            self.density_estimator = density_estimator
 
    def fit(self, X, y):
        """Fit the model"""
        X = array2d(X)
        y = np.asarray(y)
        self.classes_ = np.sort(np.unique(y))
        n_classes = len(self.classes_)
        n_samples, n_features = X.shape
        
        self.n_features_ = X.shape[1]
        self.class_prior_ = np.array([np.float(np.sum(y == y_i)) / n_samples
                                      for y_i in self.classes_])
        self.estimators_ = [clone(self.density_estimator).fit(X[y == c])
                            for c in self.classes_]
        return self
 
    def _joint_log_likelihood(self, X):
        """Compute the per-class log likelihood of each sample

        Parameters
        ----------
        X : array_like
            Array of samples on which to compute likelihoods.  Shape is
            (n_samples, n_features)

        Returns
        -------
        logL : array_like
            The log likelihood under each class.
            Shape is (n_samples, n_classes).  logL[i, j] gives the log
            likelihood of X[i] within the model representing the class
            self.classes_[j].
        """
        X = array2d(X)

        # GMM API, in particular score() and score_samples(), is
        # not consistent with the rest of the package.  This needs
        # to be addressed eventually...
        if isinstance(self.density_estimator, GMM):
            return np.array([np.log(prior) + dens.score(X)
                             for (prior, dens)
                             in zip(self.class_prior_,
                                    self.estimators_)]).T
        else:
            return np.array([np.log(prior) + dens.score_samples(X)
                             for (prior, dens)
                             in zip(self.class_prior_,
                                    self.estimators_)]).T

    def sample(self, n_samples=1, random_state=None):
        """Generate random samples from the model.

        Parameters
        ----------
        n_samples : int, optional
            Number of samples to generate. Defaults to 1.

        Returns
        -------
        X : array_like, shape (n_samples, n_features)
            List of samples
        y : array_like, shape (n_samples,)
            List of class labels for the generated samples
        """
        random_state = check_random_state(random_state)
        X = np.empty((n_samples, self.n_features_))
        rand = random_state.rand(n_samples)

        # split samples by class
        prior_cdf = np.cumsum(self.class_prior_)
        labels = prior_cdf.searchsorted(rand)

        print prior_cdf

        # for each class, generate all needed samples
        for i, model in enumerate(self.estimators_):
            model_mask = (labels == i)
            N_model = model_mask.sum()
            if N_model > 0:
                X[model_mask] = model.sample(N_model,
                                             random_state=random_state)

        return X, self.classes_[labels]

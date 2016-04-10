"""
Smoothness maximization is a graph-based method for semi-supervised
learning. It is a simple closed-form algorithm which yields encouraging
experimental results.

The basic idea is "to construct the classifying function which is
sufficiently smooth with respect to the intrinsic global structure
collectively revealed by known labeled and unlabeled points".

Notes
-----
In accordance with sklearn common practice, we use labels 0 and 1 for
negative and positive samples respectively, instead of -1 and 1
described in the original literature. So the marker value for unlabeled
samples are -1 rather than 0.

References:
[1] Zhou, Dengyong, et al. "Learning with local and global
consistency." Advances in neural information processing systems
16.16 (2004): 321-328.
[2] Zhou, Dengyong, et al. "Semi-supervised learning by maximizing
smoothness." J. of Mach. Learn. Research (2004).

"""

# Authors: Boyuan Deng <contact@boyuandeng.me>
# Licence: BSD

import numpy as np
from scipy import linalg

from ..base import BaseEstimator, ClassifierMixin
from ..metrics.pairwise import rbf_kernel
from ..utils.graph import graph_laplacian
from ..utils.multiclass import unique_labels
from ..utils.validation import check_array, check_is_fitted, check_X_y


class SmoothnessMaximization(BaseEstimator, ClassifierMixin):
    """Smoothness maximization classifier.

    Parameters
    ----------
    gamma : float
        Parameter for rbf kernel.

    alpha : float
        Specify the relative amount of the information from its
        neighbors and its initial label information.

    Attributes
    ----------
    indirect_sort : ndarray, shape = (n_samples,)
        The indices that would partition the samples into unlabeled
        ones following labeled ones (the original input is not
        necessarily partitioned).

    f : ndarray, shape = (n_samples,)
        Values of the real-valued classifying function given
        partitioned samples as input.

    """

    def __init__(self, gamma=20, alpha=0.99):
        self.gamma = gamma
        self.alpha = alpha

        self.indirect_sort = None
        self.f = None

    def fit(self, X, y):
        """Fit a smoothness maximization model.

        The input data is provided as matrix X (labeled and unlabeled)
        and corresponding label matrix y with a dedicated marker value
        for unlabeled samples.

        Parameters
        ----------
        X : ndarray, shape = (n_samples, n_features)
            Samples.

        y : ndarray, shape = (n_samples,)
            Labels.
            Negative and positive samples have labels 0 and 1
            respectively. Unlabeled samples are marked as -1.

        Returns
        -------
        self : returns an instance of self.

        """
        X, y = check_X_y(X, y)
        self.classes_ = unique_labels(y)
        self.X_ = X
        self.y_ = y

        # move all unlabeled samples behind labeled ones
        self.indirect_sort = np.argsort(y)[::-1]
        X_sorted = self.X_[self.indirect_sort]

        # convert labels according to literature's convention, which is
        # crucial because the values are involved in matrix computation
        y_sorted = self.y_[self.indirect_sort]
        negative_samples_mask = y_sorted == 0
        unlabeled_samples_mask = y_sorted == -1
        y_sorted[negative_samples_mask] = -1
        y_sorted[unlabeled_samples_mask] = 0

        # affinity matrix construction
        affinity_matrix = rbf_kernel(X_sorted, X_sorted, gamma=self.gamma)
        np.fill_diagonal(affinity_matrix, 0.0)

        # matrix S and computation of f
        laplacian = graph_laplacian(affinity_matrix, normed=True)
        I = np.ones_like(laplacian)
        S = I - laplacian
        self.f = np.dot(linalg.inv(I - self.alpha * S), y_sorted)
        return self

    def predict(self, X):
        """Predict the labels of unlabeled samples.

        Parameters
        ----------
        X : ndarray, shape = (n_samples, n_features)
            The same set of samples used in fitting.
            Smoothness maximization is not supposed to work on new
            samples.

        Returns
        -------
        C : ndarray, shape = (n_samples,)
            Labels with markers for unlabeled samples replaced by
            predictions.

        """
        check_is_fitted(self, ["X_", "indirect_sort", "f"])
        X = check_array(X)

        # find the mapping to "undo" the partition
        n_samples = self.f.shape[0]
        undo_sort = np.zeros(n_samples, dtype=int)
        undo_sort[self.indirect_sort] = np.arange(n_samples)

        classes = np.sign(self.f[undo_sort])
        classes[classes == -1] = 0
        return classes

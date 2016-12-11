# -*- coding: utf-8 -*-
"""Generic feature selection mixin"""

# Authors: G. Varoquaux, A. Gramfort, L. Buitinck, J. Nothman
# License: BSD 3 clause

from abc import ABCMeta, abstractmethod
from warnings import warn

import numpy as np
from scipy.sparse import issparse, csc_matrix

from ..base import TransformerMixin
from ..utils import check_array, safe_mask, check_X_y
from ..externals import six


class SelectorMixin(six.with_metaclass(ABCMeta, TransformerMixin)):
    """
    Tranformer mixin that performs feature selection given a support mask

    This mixin provides a feature selector implementation with `transform` and
    `inverse_transform` functionality given an implementation of
    `_get_support_mask`.
    """

    def get_support(self, indices=False):
        """
        Get a mask, or integer index, of the features selected

        Parameters
        ----------
        indices : boolean (default False)
            If True, the return value will be an array of integers, rather
            than a boolean mask.

        Returns
        -------
        support : array
            An index that selects the retained features from a feature vector.
            If `indices` is False, this is a boolean array of shape
            [# input features], in which an element is True iff its
            corresponding feature is selected for retention. If `indices` is
            True, this is an integer array of shape [# output features] whose
            values are indices into the input feature vector.
        """
        mask = self._get_support_mask()
        return mask if not indices else np.where(mask)[0]

    @abstractmethod
    def _get_support_mask(self):
        """
        Get the boolean mask indicating which features are selected

        Returns
        -------
        support : boolean array of shape [# input features]
            An element is True iff its corresponding feature is selected for
            retention.
        """

    def transform(self, X):
        """Reduce X to the selected features.

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            The input samples.

        Returns
        -------
        X_r : array of shape [n_samples, n_selected_features]
            The input samples with only the selected features.
        """
        X = check_array(X, accept_sparse='csr')
        mask = self.get_support()
        if not mask.any():
            warn("No features were selected: either the data is"
                 " too noisy or the selection test too strict.",
                 UserWarning)
            return np.empty(0).reshape((X.shape[0], 0))
        if len(mask) != X.shape[1]:
            raise ValueError("X has a different shape than during fitting.")
        return X[:, safe_mask(X, mask)]

    def inverse_transform(self, X):
        """
        Reverse the transformation operation

        Parameters
        ----------
        X : array of shape [n_samples, n_selected_features]
            The input samples.

        Returns
        -------
        X_r : array of shape [n_samples, n_original_features]
            `X` with columns of zeros inserted where features would have
            been removed by `transform`.
        """
        if issparse(X):
            X = X.tocsc()
            # insert additional entries in indptr:
            # e.g. if transform changed indptr from [0 2 6 7] to [0 2 3]
            # col_nonzeros here will be [2 0 1] so indptr becomes [0 2 2 3]
            it = self.inverse_transform(np.diff(X.indptr).reshape(1, -1))
            col_nonzeros = it.ravel()
            indptr = np.concatenate([[0], np.cumsum(col_nonzeros)])
            Xt = csc_matrix((X.data, X.indices, indptr),
                            shape=(X.shape[0], len(indptr) - 1), dtype=X.dtype)
            return Xt

        support = self.get_support()
        X = check_array(X)
        if support.sum() != X.shape[1]:
            raise ValueError("X has a different shape than during fitting.")

        if X.ndim == 1:
            X = X[None, :]
        Xt = np.zeros((X.shape[0], support.size), dtype=X.dtype)
        Xt[:, support] = X
        return Xt


def wrapper_scorer(score_func, X, y=None):
    """ A wrapper function around score functions. This function takes as
        input the feature matrix X and/or the target vector y and sends it
        as parameter to the score function.

    Parameters
    ----------
    score_func : callable
        Function taking array(s) X and/or y, and returning a pair of arrays
        (scores, pvalues) or a single array with scores.
        
    X : array-like, shape = [n_samples, n_features]
        The training input samples/Feature matrix.

    y : array-like or None, shape = [n_samples]
        The target values (class labels in classification, real numbers in
        regression).

    Notes
    -----
    The negative score values returned by a score function are changed to
    zero as both negative and zero score values signify the same thing.
    E.g. Mutual information between two random variables is a non-negative
    value, which measures the dependency between the variables. It is equal
    to zero if and only if two random variables are independent, and higher
    values mean higher dependency.
    """

    if not callable(score_func):
        raise TypeError("The score function should be a callable, %s (%s) "
                        "was passed."
                        % (score_func, type(score_func)))

    if y == None:
        X = check_array(X, ('csr', 'csc'))
    else:
        X, y = check_X_y(X, y, ['csr', 'csc'], multi_output=True)

    if y == None:
        score_func_ret = score_func(X)
    else:
        score_func_ret = score_func(X, y)

    if isinstance(score_func_ret, (list, tuple)):
        scores, pvalues = score_func_ret
        pvalues = np.asarray(pvalues)
    else:
        scores = score_func_ret
        pvalues = None
    scores = np.asarray(scores)
    scores[scores < 0.0] = 0.0

    if pvalues == None:
        return scores
    else:
        return scores, pvalues
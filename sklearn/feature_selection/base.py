# -*- coding: utf-8 -*-
"""Generic feature selection mixin"""

# Authors: G. Varoquaux, A. Gramfort, L. Buitinck, J. Nothman
# License: BSD 3 clause

from abc import ABCMeta, abstractmethod
from warnings import warn

import numpy as np
from scipy.sparse import issparse, csc_matrix

from ..base import TransformerMixin
from ..utils import check_array, check_X_y, safe_mask
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


def featurewise_scorer(score_func, **kwargs):
    """ A wrapper function around score functions.

    Parameters
    ----------
    score_func : callable
        Function taking array(s) X and/or y, and returning a pair of arrays
        (scores, pvalues) or a single array with scores. This function is also
        allowed to take other parameters as input.

    Returns
    -------
    scores : array-like, shape (n_features,)
        Score values returned by the scoring function.
    p_vals : array-like, shape (n_features,)
        The set of p-values returned by the scoring function.

    Notes
    -----
    This wrapper function wraps around scoring functions like `spearmanr`,
    `pearsonr` etc. from the scipy.stats module and makes it usable for
    feature selection algorithms like `SelectKBest`. Also, this wrapper
    function takes the absolute value of the scores returned i.e. a score of
    +1 is same as -1.

    Example
    -------
    >>> from sklearn.feature_selection import featurewise_scorer, SelectKBest
    >>> from scipy.stats import spearmanr
    >>> from sklearn.datasets import make_classification
    >>> X, y = make_classification(random_state=0)
    >>> skb = SelectKBest(featurewise_scorer(spearmanr, nan_policy='propagate'), k=10)
    >>> skb.fit(X, y) #doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    SelectKBest(k=10, score_func=...)
    >>> new_X = skb.transform(X)

    """
    def call_scorer(X, y=None):
        if y is None:
            X = check_array(X, ('csr', 'csc'))
        else:
            X, y = check_X_y(X, y, ('csr', 'csc'), multi_output=True)

        scores = []
        p_vals = []

        for i in six.moves.range(X.shape[1]):
            if y is None:
                score_func_ret = score_func(X[:, i], **kwargs)
            else:
                score_func_ret = score_func(X[:, i], y, **kwargs)

            if isinstance(score_func_ret, (list, tuple)):
                score, p_val = score_func_ret
                p_vals.append(p_val)
            else:
                score = score_func_ret
            scores.append(abs(score))

        scores = np.asarray(scores)
        if len(p_vals) > 0:
            p_vals = np.asarray(p_vals)
            return (scores, p_vals)
        else:
            return scores

    return call_scorer

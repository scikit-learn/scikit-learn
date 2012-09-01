# -*- coding: utf-8 -*-
"""Univariate features selection."""

# Authors: V. Michel, B. Thirion, G. Varoquaux, A. Gramfort, E. Duchesnay.
#          L. Buitinck
# License: BSD 3 clause

from abc import ABCMeta, abstractmethod

import numpy as np
from scipy import stats
from scipy.sparse import issparse

from ..base import BaseEstimator, TransformerMixin
from ..preprocessing import LabelBinarizer
from ..utils import array2d, atleast2d_or_csr, deprecated, \
        check_arrays, safe_asarray, safe_sqr, safe_mask
from ..utils.extmath import safe_sparse_dot

######################################################################
# Scoring functions


# The following function is a rewriting of scipy.stats.f_oneway
# Contrary to the scipy.stats.f_oneway implementation it does not
# copy the data while keeping the inputs unchanged.
def f_oneway(*args):
    """Performs a 1-way ANOVA.

    The one-way ANOVA tests the null hypothesis that 2 or more groups have
    the same population mean.  The test is applied to samples from two or
    more groups, possibly with differing sizes.

    Parameters
    ----------
    sample1, sample2, ... : array_like, sparse matrices
        The sample measurements should be given as arguments.

    Returns
    -------
    F-value : float
        The computed F-value of the test
    p-value : float
        The associated p-value from the F-distribution

    Notes
    -----
    The ANOVA test has important assumptions that must be satisfied in order
    for the associated p-value to be valid.

    1. The samples are independent
    2. Each sample is from a normally distributed population
    3. The population standard deviations of the groups are all equal. This
       property is known as homoscedasticity.

    If these assumptions are not true for a given set of data, it may still be
    possible to use the Kruskal-Wallis H-test (`scipy.stats.kruskal`_) although
    with some loss of power.

    The algorithm is from Heiman[2], pp.394-7.

    See ``scipy.stats.f_oneway`` that should give the same results while
    being less efficient.

    References
    ----------

    .. [1] Lowry, Richard.  "Concepts and Applications of Inferential
           Statistics". Chapter 14.
           http://faculty.vassar.edu/lowry/ch14pt1.html

    .. [2] Heiman, G.W.  Research Methods in Statistics. 2002.

    """
    n_classes = len(args)
    args = [safe_asarray(a) for a in args]
    n_samples_per_class = np.array([a.shape[0] for a in args])
    n_samples = np.sum(n_samples_per_class)
    ss_alldata = reduce(lambda x, y: x + y,
            [safe_sqr(a).sum(axis=0) for a in args])
    sums_args = [a.sum(axis=0) for a in args]
    square_of_sums_alldata = safe_sqr(reduce(lambda x, y: x + y, sums_args))
    square_of_sums_args = [safe_sqr(s) for s in sums_args]
    sstot = ss_alldata - square_of_sums_alldata / float(n_samples)
    ssbn = 0.
    for k, _ in enumerate(args):
        ssbn += square_of_sums_args[k] / n_samples_per_class[k]
    ssbn -= square_of_sums_alldata / float(n_samples)
    sswn = sstot - ssbn
    dfbn = n_classes - 1
    dfwn = n_samples - n_classes
    msb = ssbn / float(dfbn)
    msw = sswn / float(dfwn)
    f = msb / msw
    # flatten matrix to vector in sparse case
    f = np.asarray(f).ravel()
    prob = stats.fprob(dfbn, dfwn, f)
    return f, prob


def f_classif(X, y):
    """Compute the Anova F-value for the provided sample

    Parameters
    ----------
    X : {array-like, sparse matrix} shape = [n_samples, n_features]
        The set of regressors that will tested sequentially
    y : array of shape(n_samples)
        The data matrix

    Returns
    -------
    F : array, shape = [n_features,]
        The set of F values
    pval : array, shape = [n_features,]
        The set of p-values
    """
    X, y = check_arrays(X, y)
    args = [X[safe_mask(X, y == k)] for k in np.unique(y)]
    return f_oneway(*args)


def chi2(X, y):
    """Compute χ² (chi-squared) statistic for each class/feature combination.

    This score can be used to select the n_features features with the
    highest values for the χ² (chi-square) statistic from either boolean or
    multinomially distributed data (e.g., term counts in document
    classification) relative to the classes.

    Recall that the χ² statistic measures dependence between stochastic
    variables, so using this function "weeds out" the features that are the
    most likely to be independent of class and therefore irrelevant for
    classification.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape = [n_samples, n_features_in]
        Sample vectors.

    y : array-like, shape = n_samples
        Target vector (class labels).

    Returns
    -------
    chi2 : array, shape = [n_features,]
        chi2 statistics of each feature
    pval : array, shape = [n_features,]
        p-values of each feature

    Notes
    -----
    Complexity of this algorithm is O(n_classes * n_features).
    """

    # XXX: we might want to do some of the following in logspace instead for
    # numerical stability.
    X = atleast2d_or_csr(X)
    Y = LabelBinarizer().fit_transform(y)
    if Y.shape[1] == 1:
        Y = np.append(1 - Y, Y, axis=1)

    observed = safe_sparse_dot(Y.T, X)          # n_classes * n_features

    feature_count = array2d(X.sum(axis=0))
    class_prob = array2d(Y.mean(axis=0))
    expected = safe_sparse_dot(class_prob.T, feature_count)

    return stats.chisquare(observed, expected)


def f_regression(X, y, center=True):
    """Univariate linear regression tests

    Quick linear model for testing the effect of a single regressor,
    sequentially for many regressors.

    This is done in 3 steps:
    1. the regressor of interest and the data are orthogonalized
    wrt constant regressors
    2. the cross correlation between data and regressors is computed
    3. it is converted to an F score then to a p-value

    Parameters
    ----------
    X : {array-like, sparse matrix}  shape = [n_samples, n_features]
        The set of regressors that will tested sequentially
    y : array of shape(n_samples)
        The data matrix

    center : True, bool,
        If true, X and y will be centered

    Returns
    -------
    F : array, shape=[m,]
        The set of F values
    pval : array, shape=[m,]
        The set of p-values
    """
    if issparse(X) and center:
        raise ValueError("center=True only allowed for dense data")
    X, y = check_arrays(X, y, dtype=np.float)
    y = y.ravel()
    if center:
        y = y - np.mean(y)
        X = X.copy('F')  # faster in fortran
        X -= X.mean(axis=0)

    # compute the correlation
    corr = safe_sparse_dot(y, X)
    corr /= np.asarray(np.sqrt(safe_sqr(X).sum(axis=0))).ravel()
    corr /= np.asarray(np.sqrt(safe_sqr(y).sum())).ravel()

    # convert to p-value
    dof = y.size - 2
    F = corr ** 2 / (1 - corr ** 2) * dof
    pv = stats.f.sf(F, 1, dof)
    return F, pv


######################################################################
# General class for filter univariate selection

class _AbstractUnivariateFilter(BaseEstimator, TransformerMixin):
    __metaclass__ = ABCMeta

    def __init__(self, score_func):
        """ Initialize the univariate feature selection.

        Parameters
        ===========
        score_func: callable
            Function taking two arrays X and y, and returning 2 arrays:
            both scores and pvalues
        """
        if not callable(score_func):
            raise TypeError(
                "The score function should be a callable, '%s' (type %s) "
                "was passed." % (score_func, type(score_func)))
        self.score_func = score_func

    def fit(self, X, y):
        """
        Evaluate the function
        """
        scores = self.score_func(X, y)
        self.scores_ = scores[0]
        self.pvalues_ = scores[1]
        return self

    @property
    @deprecated('``_scores`` is deprecated and will be removed in '
                'version 0.13. Please use ``scores_`` instead.')
    def _scores(self):
        return self.scores_

    @property
    @deprecated('``_pvalues`` is deprecated and will be removed in '
                'version 0.13. Please use ``pvalues_`` instead.')
    def _pvalues(self):
        return self.pvalues_

    def get_support(self, indices=False):
        """
        Return a mask, or list, of the features/indices selected.
        """
        mask = self._get_support_mask()
        return mask if not indices else np.where(mask)[0]

    @abstractmethod
    def _get_support_mask(self):
        """
        Must return a boolean mask indicating which features are selected.
        """

    def transform(self, X):
        """
        Transform a new matrix using the selected features
        """
        X = atleast2d_or_csr(X)
        mask = self._get_support_mask()
        if len(mask) != X.shape[1]:
            raise ValueError("X has a different shape than during fitting.")
        return atleast2d_or_csr(X)[:, safe_mask(X, mask)]

    def inverse_transform(self, X):
        """
        Transform a new matrix using the selected features
        """
        support_ = self.get_support()
        if X.ndim == 1:
            X = X[None, :]
        Xt = np.zeros((X.shape[0], support_.size))
        Xt[:, support_] = X
        return Xt


######################################################################
# Specific filters
######################################################################

class SelectPercentile(_AbstractUnivariateFilter):
    """Filter: Select the best percentile of the p_values

    Parameters
    ===========
    score_func: callable
        Function taking two arrays X and y, and returning 2 arrays:
        both scores and pvalues

    percentile: int, optional
        Percent of features to keep

    """

    def __init__(self, score_func=f_classif, percentile=10):
        self.percentile = percentile
        _AbstractUnivariateFilter.__init__(self, score_func)

    def _get_support_mask(self):
        percentile = self.percentile
        if percentile > 100:
            raise ValueError("percentile should be between 0 and 100"
                             " (%f given)" % (percentile))
        # Cater for Nans
        if percentile == 100:
            return np.ones(len(self.pvalues_), dtype=np.bool)
        elif percentile == 0:
            return np.zeros(len(self.pvalues_), dtype=np.bool)
        alpha = stats.scoreatpercentile(self.pvalues_, percentile)
        return (self.pvalues_ <= alpha)


class SelectKBest(_AbstractUnivariateFilter):
    """Filter: Select the k lowest p-values.

    Parameters
    ----------
    score_func: callable
        Function taking two arrays X and y, and returning 2 arrays:
        both scores and pvalues

    k: int, optional
        Number of top features to select.

    Notes
    -----
    Ties between features with equal p-values will be broken in an unspecified
    way.

    """

    def __init__(self, score_func=f_classif, k=10):
        self.k = k
        _AbstractUnivariateFilter.__init__(self, score_func)

    def _get_support_mask(self):
        k = self.k
        if k > len(self.pvalues_):
            raise ValueError("cannot select %d features among %d"
                             % (k, len(self.pvalues_)))

        # XXX This should be refactored; we're getting an array of indices
        # from argsort, which we transform to a mask, which we probably
        # transform back to indices later.
        mask = np.zeros(self.pvalues_.shape, dtype=bool)
        mask[np.argsort(self.pvalues_)[:k]] = 1
        return mask


class SelectFpr(_AbstractUnivariateFilter):
    """Filter: Select the pvalues below alpha based on a FPR test.

    FPR test stands for False Positive Rate test. It controls the total
    amount of false detections.

    Parameters
    ----------
    score_func: callable
        Function taking two arrays X and y, and returning 2 arrays:
        both scores and pvalues

    alpha: float, optional
        The highest p-value for features to be kept
    """

    def __init__(self, score_func=f_classif, alpha=5e-2):
        self.alpha = alpha
        _AbstractUnivariateFilter.__init__(self, score_func)

    def _get_support_mask(self):
        alpha = self.alpha
        return self.pvalues_ < alpha


class SelectFdr(_AbstractUnivariateFilter):
    """Filter: Select the p-values for an estimated false discovery rate

    This uses the Benjamini-Hochberg procedure. ``alpha`` is the target false
    discovery rate.

    Parameters
    ----------
    score_func: callable
        Function taking two arrays X and y, and returning 2 arrays:
        both scores and pvalues

    alpha: float, optional
        The highest uncorrected p-value for features to keep

    """

    def __init__(self, score_func=f_classif, alpha=5e-2):
        self.alpha = alpha
        _AbstractUnivariateFilter.__init__(self, score_func)

    def _get_support_mask(self):
        alpha = self.alpha
        sv = np.sort(self.pvalues_)
        threshold = sv[sv < alpha * np.arange(len(self.pvalues_))].max()
        return self.pvalues_ <= threshold


class SelectFwe(_AbstractUnivariateFilter):
    """Filter: Select the p-values corresponding to Family-wise error rate

    Parameters
    ----------
    score_func: callable
        Function taking two arrays X and y, and returning 2 arrays:
        both scores and pvalues

    alpha: float, optional
        The highest uncorrected p-value for features to keep

    """

    def __init__(self, score_func=f_classif, alpha=5e-2):
        self.alpha = alpha
        _AbstractUnivariateFilter.__init__(self, score_func)

    def _get_support_mask(self):
        alpha = self.alpha
        return (self.pvalues_ < alpha / len(self.pvalues_))


######################################################################
# Generic filter
######################################################################

class GenericUnivariateSelect(_AbstractUnivariateFilter):
    """Univariate feature selector with configurable strategy

    Parameters
    ----------
    score_func: callable
        Function taking two arrays X and y, and returning 2 arrays:
        both scores and pvalues

    mode: {'percentile', 'k_best', 'fpr', 'fdr', 'fwe'}
        Feature selection mode

    param: float or int depending on the feature selection mode
        Parameter of the corresponding mode
    """

    _selection_modes = {'percentile':   SelectPercentile,
                        'k_best':       SelectKBest,
                        'fpr':          SelectFpr,
                        'fdr':          SelectFdr,
                        'fwe':          SelectFwe,
                        }

    def __init__(self, score_func=f_classif, mode='percentile', param=1e-5):
        if not callable(score_func):
            raise TypeError(
                "The score function should be a callable, '%s' (type %s) "
                "was passed." % (score_func, type(score_func)))
        if mode not in self._selection_modes:
            raise ValueError(
                "The mode passed should be one of %s, '%s', (type %s) "
                "was passed." % (
                        self._selection_modes.keys(),
                        mode, type(mode)))
        self.score_func = score_func
        self.mode = mode
        self.param = param

    def _get_support_mask(self):
        selector = self._selection_modes[self.mode](lambda x: x)
        selector.pvalues_ = self.pvalues_
        selector.scores_ = self.scores_
        # Now perform some acrobatics to set the right named parameter in
        # the selector
        possible_params = selector._get_param_names()
        possible_params.remove('score_func')
        selector.set_params(**{possible_params[0]: self.param})
        return selector._get_support_mask()

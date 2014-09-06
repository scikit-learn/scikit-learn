"""Univariate features selection."""

# Authors: V. Michel, B. Thirion, G. Varoquaux, A. Gramfort, E. Duchesnay.
#          L. Buitinck, A. Joly
# License: BSD 3 clause


import numpy as np
from scipy import special, stats
from scipy.sparse import issparse

from ..base import BaseEstimator
from ..preprocessing import LabelBinarizer
from ..utils import (as_float_array, check_array, check_X_y, safe_sqr,
                     safe_mask)
from ..utils.extmath import norm, safe_sparse_dot
from .base import SelectorMixin


def _clean_nans(scores):
    """
    Fixes Issue #1240: NaNs can't be properly compared, so change them to the
    smallest value of scores's dtype. -inf seems to be unreliable.
    """
    # XXX where should this function be called? fit? scoring functions
    # themselves?
    scores = as_float_array(scores, copy=True)
    scores[np.isnan(scores)] = np.finfo(scores.dtype).min
    return scores


######################################################################
# Scoring functions


# The following function is a rewriting of scipy.stats.f_oneway
# Contrary to the scipy.stats.f_oneway implementation it does not
# copy the data while keeping the inputs unchanged.
def f_oneway(*args):
    """Performs a 1-way ANOVA.

    The one-way ANOVA tests the null hypothesis that 2 or more groups have
    the same population mean. The test is applied to samples from two or
    more groups, possibly with differing sizes.

    Parameters
    ----------
    sample1, sample2, ... : array_like, sparse matrices
        The sample measurements should be given as arguments.

    Returns
    -------
    F-value : float
        The computed F-value of the test.
    p-value : float
        The associated p-value from the F-distribution.

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
    args = [as_float_array(a) for a in args]
    n_samples_per_class = np.array([a.shape[0] for a in args])
    n_samples = np.sum(n_samples_per_class)
    ss_alldata = sum(safe_sqr(a).sum(axis=0) for a in args)
    sums_args = [np.asarray(a.sum(axis=0)) for a in args]
    square_of_sums_alldata = sum(sums_args) ** 2
    square_of_sums_args = [s ** 2 for s in sums_args]
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
        The set of regressors that will tested sequentially.

    y : array of shape(n_samples)
        The data matrix.

    Returns
    -------
    F : array, shape = [n_features,]
        The set of F values.

    pval : array, shape = [n_features,]
        The set of p-values.
    """
    X, y = check_X_y(X, y, ['csr', 'csc', 'coo'])
    args = [X[safe_mask(X, y == k)] for k in np.unique(y)]
    return f_oneway(*args)


def _chisquare(f_obs, f_exp):
    """Fast replacement for scipy.stats.chisquare.

    Version from https://github.com/scipy/scipy/pull/2525 with additional
    optimizations.
    """
    f_obs = np.asarray(f_obs, dtype=np.float64)

    k = len(f_obs)
    # Reuse f_obs for chi-squared statistics
    chisq = f_obs
    chisq -= f_exp
    chisq **= 2
    chisq /= f_exp
    chisq = chisq.sum(axis=0)
    return chisq, special.chdtrc(k - 1, chisq)


def chi2(X, y):
    """Compute chi-squared statistic for each class/feature combination.

    This score can be used to select the n_features features with the
    highest values for the test chi-squared statistic from X, which must
    contain booleans or frequencies (e.g., term counts in document
    classification), relative to the classes.

    Recall that the chi-square test measures dependence between stochastic
    variables, so using this function "weeds out" the features that are the
    most likely to be independent of class and therefore irrelevant for
    classification.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape = (n_samples, n_features_in)
        Sample vectors.

    y : array-like, shape = (n_samples,)
        Target vector (class labels).

    Returns
    -------
    chi2 : array, shape = (n_features,)
        chi2 statistics of each feature.
    pval : array, shape = (n_features,)
        p-values of each feature.

    Notes
    -----
    Complexity of this algorithm is O(n_classes * n_features).
    """

    # XXX: we might want to do some of the following in logspace instead for
    # numerical stability.
    X = check_array(X, accept_sparse='csr')
    if np.any((X.data if issparse(X) else X) < 0):
        raise ValueError("Input X must be non-negative.")

    Y = LabelBinarizer().fit_transform(y)
    if Y.shape[1] == 1:
        Y = np.append(1 - Y, Y, axis=1)

    observed = safe_sparse_dot(Y.T, X)          # n_classes * n_features

    feature_count = check_array(X.sum(axis=0))
    class_prob = check_array(Y.mean(axis=0))
    expected = np.dot(class_prob.T, feature_count)

    return _chisquare(observed, expected)


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
    X : {array-like, sparse matrix}  shape = (n_samples, n_features)
        The set of regressors that will tested sequentially.

    y : array of shape(n_samples).
        The data matrix

    center : True, bool,
        If true, X and y will be centered.

    Returns
    -------
    F : array, shape=(n_features,)
        F values of features.

    pval : array, shape=(n_features,)
        p-values of F-scores.
    """
    if issparse(X) and center:
        raise ValueError("center=True only allowed for dense data")
    X, y = check_X_y(X, y, ['csr', 'csc', 'coo'], dtype=np.float)
    if center:
        y = y - np.mean(y)
        X = X.copy('F')  # faster in fortran
        X -= X.mean(axis=0)

    # compute the correlation
    corr = safe_sparse_dot(y, X)
    # XXX could use corr /= row_norms(X.T) here, but the test doesn't pass
    corr /= np.asarray(np.sqrt(safe_sqr(X).sum(axis=0))).ravel()
    corr /= norm(y)

    # convert to p-value
    degrees_of_freedom = y.size - (2 if center else 1)
    F = corr ** 2 / (1 - corr ** 2) * degrees_of_freedom
    pv = stats.f.sf(F, 1, degrees_of_freedom)
    return F, pv


######################################################################
# Base classes

class _BaseFilter(BaseEstimator, SelectorMixin):
    """Initialize the univariate feature selection.

    Parameters
    ----------
    score_func : callable
        Function taking two arrays X and y, and returning a pair of arrays
        (scores, pvalues).
    """

    def __init__(self, score_func):
        self.score_func = score_func

    def fit(self, X, y):
        """Run score function on (X, y) and get the appropriate features.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values (class labels in classification, real numbers in
            regression).

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y, ['csr', 'csc', 'coo'])

        if not callable(self.score_func):
            raise TypeError("The score function should be a callable, %s (%s) "
                            "was passed."
                            % (self.score_func, type(self.score_func)))

        self._check_params(X, y)

        self.scores_, self.pvalues_ = self.score_func(X, y)
        self.scores_ = np.asarray(self.scores_)
        self.pvalues_ = np.asarray(self.pvalues_)
        return self

    def _check_params(self, X, y):
        pass


######################################################################
# Specific filters
######################################################################
class SelectPercentile(_BaseFilter):
    """Select features according to a percentile of the highest scores.

    Parameters
    ----------
    score_func : callable
        Function taking two arrays X and y, and returning a pair of arrays
        (scores, pvalues).

    percentile : int, optional, default=10
        Percent of features to keep.

    Attributes
    ----------
    scores_ : array-like, shape=(n_features,)
        Scores of features.

    pvalues_ : array-like, shape=(n_features,)
        p-values of feature scores.

    Notes
    -----
    Ties between features with equal scores will be broken in an unspecified
    way.

    """

    def __init__(self, score_func=f_classif, percentile=10):
        super(SelectPercentile, self).__init__(score_func)
        self.percentile = percentile

    def _check_params(self, X, y):
        if not 0 <= self.percentile <= 100:
            raise ValueError("percentile should be >=0, <=100; got %r"
                             % self.percentile)

    def _get_support_mask(self):
        # Cater for NaNs
        if self.percentile == 100:
            return np.ones(len(self.scores_), dtype=np.bool)
        elif self.percentile == 0:
            return np.zeros(len(self.scores_), dtype=np.bool)

        scores = _clean_nans(self.scores_)
        treshold = stats.scoreatpercentile(scores,
                                           100 - self.percentile)
        mask = scores > treshold
        ties = np.where(scores == treshold)[0]
        if len(ties):
            max_feats = len(scores) * self.percentile // 100
            kept_ties = ties[:max_feats - mask.sum()]
            mask[kept_ties] = True
        return mask


class SelectKBest(_BaseFilter):
    """Select features according to the k highest scores.

    Parameters
    ----------
    score_func : callable
        Function taking two arrays X and y, and returning a pair of arrays
        (scores, pvalues).

    k : int or "all", optional, default=10
        Number of top features to select.
        The "all" option bypasses selection, for use in a parameter search.

    Attributes
    ----------
    scores_ : array-like, shape=(n_features,)
        Scores of features.

    pvalues_ : array-like, shape=(n_features,)
        p-values of feature scores.

    Notes
    -----
    Ties between features with equal scores will be broken in an unspecified
    way.

    """

    def __init__(self, score_func=f_classif, k=10):
        super(SelectKBest, self).__init__(score_func)
        self.k = k

    def _check_params(self, X, y):
        if not (self.k == "all" or 0 <= self.k <= X.shape[1]):
            raise ValueError("k should be >=0, <= n_features; got %r."
                             "Use k='all' to return all features."
                             % self.k)

    def _get_support_mask(self):
        if self.k == 'all':
            return np.ones(self.scores_.shape, dtype=bool)
        elif self.k == 0:
            return np.zeros(self.scores_.shape, dtype=bool)
        else:
            scores = _clean_nans(self.scores_)
            mask = np.zeros(scores.shape, dtype=bool)

            # Request a stable sort. Mergesort takes more memory (~40MB per
            # megafeature on x86-64).
            mask[np.argsort(scores, kind="mergesort")[-self.k:]] = 1
            return mask


class SelectFpr(_BaseFilter):
    """Filter: Select the pvalues below alpha based on a FPR test.

    FPR test stands for False Positive Rate test. It controls the total
    amount of false detections.

    Parameters
    ----------
    score_func : callable
        Function taking two arrays X and y, and returning a pair of arrays
        (scores, pvalues).

    alpha : float, optional
        The highest p-value for features to be kept.

    Attributes
    ----------
    scores_ : array-like, shape=(n_features,)
        Scores of features.

    pvalues_ : array-like, shape=(n_features,)
        p-values of feature scores.
    """

    def __init__(self, score_func=f_classif, alpha=5e-2):
        super(SelectFpr, self).__init__(score_func)
        self.alpha = alpha

    def _get_support_mask(self):
        return self.pvalues_ < self.alpha


class SelectFdr(_BaseFilter):
    """Filter: Select the p-values for an estimated false discovery rate

    This uses the Benjamini-Hochberg procedure. ``alpha`` is the target false
    discovery rate.

    Parameters
    ----------
    score_func : callable
        Function taking two arrays X and y, and returning a pair of arrays
        (scores, pvalues).

    alpha : float, optional
        The highest uncorrected p-value for features to keep.


    Attributes
    ----------
    scores_ : array-like, shape=(n_features,)
        Scores of features.

    pvalues_ : array-like, shape=(n_features,)
        p-values of feature scores.
    """

    def __init__(self, score_func=f_classif, alpha=5e-2):
        super(SelectFdr, self).__init__(score_func)
        self.alpha = alpha

    def _get_support_mask(self):
        alpha = self.alpha
        sv = np.sort(self.pvalues_)
        threshold = sv[sv < alpha * np.arange(len(self.pvalues_))].max()
        return self.pvalues_ <= threshold


class SelectFwe(_BaseFilter):
    """Filter: Select the p-values corresponding to Family-wise error rate

    Parameters
    ----------
    score_func : callable
        Function taking two arrays X and y, and returning a pair of arrays
        (scores, pvalues).

    alpha : float, optional
        The highest uncorrected p-value for features to keep.

    Attributes
    ----------
    scores_ : array-like, shape=(n_features,)
        Scores of features.

    pvalues_ : array-like, shape=(n_features,)
        p-values of feature scores.
    """

    def __init__(self, score_func=f_classif, alpha=5e-2):
        super(SelectFwe, self).__init__(score_func)
        self.alpha = alpha

    def _get_support_mask(self):
        return (self.pvalues_ < self.alpha / len(self.pvalues_))


######################################################################
# Generic filter
######################################################################

# TODO this class should fit on either p-values or scores,
# depending on the mode.
class GenericUnivariateSelect(_BaseFilter):
    """Univariate feature selector with configurable strategy.

    Parameters
    ----------
    score_func : callable
        Function taking two arrays X and y, and returning a pair of arrays
        (scores, pvalues).

    mode : {'percentile', 'k_best', 'fpr', 'fdr', 'fwe'}
        Feature selection mode.

    param : float or int depending on the feature selection mode
        Parameter of the corresponding mode.

    Attributes
    ----------
    scores_ : array-like, shape=(n_features,)
        Scores of features.

    pvalues_ : array-like, shape=(n_features,)
        p-values of feature scores.
    """

    _selection_modes = {'percentile':   SelectPercentile,
                        'k_best':       SelectKBest,
                        'fpr':          SelectFpr,
                        'fdr':          SelectFdr,
                        'fwe':          SelectFwe}

    def __init__(self, score_func=f_classif, mode='percentile', param=1e-5):
        super(GenericUnivariateSelect, self).__init__(score_func)
        self.mode = mode
        self.param = param

    def _make_selector(self):
        selector = self._selection_modes[self.mode](score_func=self.score_func)

        # Now perform some acrobatics to set the right named parameter in
        # the selector
        possible_params = selector._get_param_names()
        possible_params.remove('score_func')
        selector.set_params(**{possible_params[0]: self.param})

        return selector

    def _check_params(self, X, y):
        if self.mode not in self._selection_modes:
            raise ValueError("The mode passed should be one of %s, %r,"
                             " (type %s) was passed."
                             % (self._selection_modes.keys(), self.mode,
                                type(self.mode)))

        self._make_selector()._check_params(X, y)

    def _get_support_mask(self):
        selector = self._make_selector()
        selector.pvalues_ = self.pvalues_
        selector.scores_ = self.scores_
        return selector._get_support_mask()

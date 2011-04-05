"""
Univariate features selection.
"""

# Authors: V. Michel, B. Thirion, G. Varoquaux, A. Gramfort, E. Duchesnay
# License: BSD 3 clause

import numpy as np
from scipy import stats

from ..base import BaseEstimator, TransformerMixin

######################################################################
# Scoring functions

# The following function is a rewriting of scipy.stats.f_oneway
# Contrary to the scipy.stats.f_oneway implementation it does not
# copy the data while keeping the inputs unchanged.
def f_oneway(*args):
    """
    Performs a 1-way ANOVA.

    The on-way ANOVA tests the null hypothesis that 2 or more groups have
    the same population mean.  The test is applied to samples from two or
    more groups, possibly with differing sizes.

    Parameters
    ----------
    sample1, sample2, ... : array_like
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
    3. The population standard deviations of the groups are all equal.  This
       property is known as homocedasticity.

    If these assumptions are not true for a given set of data, it may still be
    possible to use the Kruskal-Wallis H-test (`stats.kruskal`_) although with
    some loss of power

    The algorithm is from Heiman[2], pp.394-7.

    See scipy.stats.f_oneway that should give the same results while
    being less efficient

    References
    ----------
    .. [1] Lowry, Richard.  "Concepts and Applications of Inferential
           Statistics". Chapter 14.
           http://faculty.vassar.edu/lowry/ch14pt1.html

    .. [2] Heiman, G.W.  Research Methods in Statistics. 2002.

    """
    n_classes = len(args)
    n_samples_per_class = np.array([len(a) for a in args])
    n_samples = np.sum(n_samples_per_class)
    ss_alldata = reduce(lambda x, y: x+y, [np.sum(a**2, axis=0) for a in args])
    sums_args = [np.sum(a, axis=0) for a in args]
    square_of_sums_alldata = reduce(lambda x, y: x+y, sums_args)**2
    square_of_sums_args = [s**2 for s in sums_args]
    sstot = ss_alldata - square_of_sums_alldata / float(n_samples)
    ssbn = 0
    for k, _ in enumerate(args):
        ssbn += square_of_sums_args[k] / n_samples_per_class[k]
    ssbn -= square_of_sums_alldata / float(n_samples)
    sswn = sstot - ssbn
    dfbn = n_classes - 1
    dfwn = n_samples - n_classes
    msb = ssbn / float(dfbn)
    msw = sswn / float(dfwn)
    f = msb / msw
    prob = stats.fprob(dfbn, dfwn, f)
    return f, prob


def f_classif(X, y):
    """
    Compute the Anova F-value for the provided sample

    Parameters
    ----------
    X : array of shape (n_samples, n_features)
        the set of regressors sthat will tested sequentially
    y : array of shape(n_samples)
        the data matrix

    Returns
    -------
    F : array of shape (m),
        the set of F values
    pval : array of shape(m),
        the set of p-values
    """
    X = np.atleast_2d(X)
    y = np.atleast_1d(y)
    if y.ndim > 1:
        y = y.ravel()
    args = [X[y==k] for k in np.unique(y)]
    return f_oneway(*args)


def f_regression(X, y, center=True):
    """
    Quick linear model for testing the effect of a single regressor,
    sequentially for many regressors
    This is done in 3 steps:
    1. the regressor of interest and the data are orthogonalized
    wrt constant regressors
    2. the cross correlation between data and regressors is computed
    3. it is converted to an F score then to a p-value

    Parameters
    ----------
    X : array of shape (n_samples, n_features)
        the set of regressors sthat will tested sequentially
    y : array of shape(n_samples)
        the data matrix

    center : True, bool,
        If true, X and y are centered

    Returns
    -------
    F : array of shape (m),
        the set of F values
    pval : array of shape(m)
        the set of p-values
    """

    # orthogonalize everything wrt to confounds
    y = y.copy().ravel()
    X = X.copy()
    if center:
        y -= np.mean(y)
        X -= np.mean(X, 0)

    # compute the correlation
    X /= np.sqrt(np.sum(X**2, 0))
    y /= np.sqrt(np.sum(y**2))
    corr = np.dot(y, X)

    # convert to p-value
    dof = y.size - 2
    F = corr**2 / (1 - corr**2) * dof
    pv = stats.f.sf(F, 1, dof)
    return F, pv

######################################################################
# General class for filter univariate selection

class _AbstractUnivariateFilter(BaseEstimator, TransformerMixin):
    """ Abstract class, not meant to be used directly
    """

    def __init__(self, score_func):
        """ Initialize the univariate feature selection.

        Parameters
        ===========
        score_func: callable
            function taking two arrays X and y, and returning 2 arrays:
            both scores and pvalues
        """
        assert callable(score_func), ValueError(
                "The score function should be a callable, '%s' (type %s) "
                "was passed." % (score_func, type(score_func)))
        self.score_func = score_func

    def fit(self, X, y, **params):
        """
        Evaluate the function
        """
        self._set_params(**params)
        _scores = self.score_func(X, y)
        self._scores = _scores[0]
        self._pvalues = _scores[1]
        return self

    def transform(self, X, **params):
        """
        Transform a new matrix using the selected features
        """
        self._set_params(**params)
        return X[:, self.get_support()]

    def inverse_transform(self, X, **params):
        """
        Transform a new matrix using the selected features
        """
        self._set_params(**params)
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
    """
    Filter : Select the best percentile of the p_values
    """

    def __init__(self, score_func, percentile=10):
        """ Initialize the univariate feature selection.

        Parameters
        ===========
        score_func: callable
            function taking two arrays X and y, and returning 2 arrays:
            both scores and pvalues
        percentile: int, optional
            percent of features to keep
        """
        self.percentile = percentile
        _AbstractUnivariateFilter.__init__(self, score_func)

    def get_support(self):
        percentile = self.percentile
        assert percentile<=100, ValueError('percentile should be \
                            between 0 and 100 (%f given)' %(percentile))
        # Cater for Nans
        if percentile == 100:
            return np.ones(len(self._pvalues), dtype=np.bool)
        elif percentile == 0:
            return np.zeros(len(self._pvalues), dtype=np.bool)
        alpha = stats.scoreatpercentile(self._pvalues, percentile)
        return (self._pvalues <= alpha)


class SelectKBest(_AbstractUnivariateFilter):
    """
    Filter : Select the k lowest p-values
    """

    def __init__(self, score_func, k=10):
        """ Initialize the univariate feature selection.

        Parameters
        ===========
        score_func: callable
            function taking two arrays X and y, and returning 2 arrays:
            both scores and pvalues
        percentile: int, optional
            percent of features to keep
        """
        self.k = k
        _AbstractUnivariateFilter.__init__(self, score_func)

    def get_support(self):
        k = self.k
        assert k<=len(self._pvalues), ValueError('cannot select %d features'
                                    ' among %d ' % (k, len(self._pvalues)))
        alpha = np.sort(self._pvalues)[k-1]
        return (self._pvalues <= alpha)


class SelectFpr(_AbstractUnivariateFilter):
    """
    Filter : Select the pvalues below alpha
    """

    def __init__(self, score_func, alpha=5e-2):
        """ Initialize the univariate feature selection.

        Parameters
        ===========
        score_func: callable
            function taking two arrays X and y, and returning 2 arrays:
            both scores and pvalues
        alpha: float, optional
            the highest p-value for features to keep
        """
        self.alpha = alpha
        _AbstractUnivariateFilter.__init__(self, score_func)

    def get_support(self):
        alpha = self.alpha
        return (self._pvalues < alpha)


class SelectFdr(_AbstractUnivariateFilter):
    """
    Filter : Select the p-values corresponding to an estimated false
    discovery rate of alpha. This uses the Benjamini-Hochberg procedure
    """

    def __init__(self, score_func, alpha=5e-2):
        """ Initialize the univariate feature selection.

        Parameters
        ===========
        score_func: callable
            function taking two arrays X and y, and returning 2 arrays:
            both scores and pvalues
        alpha: float, optional
            the highest uncorrected p-value for features to keep
        """
        self.alpha = alpha
        _AbstractUnivariateFilter.__init__(self, score_func)

    def get_support(self):
        alpha = self.alpha
        sv = np.sort(self._pvalues)
        threshold = sv[sv < alpha*np.arange(len(self._pvalues))].max()
        return (self._pvalues < threshold)


class SelectFwe(_AbstractUnivariateFilter):
    """
    Filter : Select the p-values corresponding to a corrected p-value of alpha
    """

    def __init__(self, score_func, alpha=5e-2):
        """ Initialize the univariate feature selection.

        Parameters
        ===========
        score_func: callable
            function taking two arrays X and y, and returning 2 arrays:
            both scores and pvalues
        alpha: float, optional
            the highest uncorrected p-value for features to keep
        """
        self.alpha = alpha
        _AbstractUnivariateFilter.__init__(self, score_func)

    def get_support(self):
        alpha = self.alpha
        return (self._pvalues < alpha/len(self._pvalues))


######################################################################
# Generic filter
######################################################################

class GenericUnivariateSelect(_AbstractUnivariateFilter):
    _selection_modes = {'percentile':   SelectPercentile,
                        'k_best':       SelectKBest,
                        'fpr':          SelectFpr,
                        'fdr':          SelectFdr,
                        'fwe':          SelectFwe,
                        }

    def __init__(self, score_func, mode='percentile', param=1e-5):
        """ Initialize the univariate feature selection.

        Parameters
        ===========
        score_func: callable
            Function taking two arrays X and y, and returning 2 arrays:
            both scores and pvalues
        mode: {%s}
            Feature selection mode
        param: float or int depending on the feature selection mode
            Parameter of the corresponding mode
        """ % self._selection_modes.keys()
        assert callable(score_func), ValueError(
                "The score function should be a callable, '%s' (type %s) "
                "was passed." % (score_func, type(score_func)))
        assert mode in self._selection_modes, ValueError(
                "The mode passed should be one of %s, '%s', (type %s) "
                "was passed." % (
                        self._selection_modes.keys(),
                        mode, type(mode)))
        self.score_func = score_func
        self.mode = mode
        self.param = param

    def get_support(self):
        selector = self._selection_modes[self.mode](lambda x: x)
        selector._pvalues = self._pvalues
        selector._scores = self._scores
        # Now make some acrobaties to set the right named parameter in
        # the selector
        possible_params = selector._get_param_names()
        possible_params.remove('score_func')
        selector._set_params(**{possible_params[0]: self.param})
        return selector.get_support()

"""
Univariate features selection.
"""

# Authors: V. Michel, B. Thirion, G. Varoquaux, A. Gramfort, E. Duchesnay
# License: BSD 3 clause

import numpy as np
from scipy import stats

from ..base import BaseEstimator

######################################################################
# Scoring functions
######################################################################

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
    args = [X[y==k] for k in np.unique(y)]
    return stats.f_oneway(*args)


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
######################################################################


class _AbstractUnivariateFilter(BaseEstimator):
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
                "was passed." % (score_func, type(score_func))
            )
        self.score_func = score_func


    def fit(self, X, y):
        """
        Evaluate the function
        """
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
                "was passed." % (score_func, type(score_func))
            )
        assert mode in self._selection_modes, ValueError(
                "The mode passed should be one of %s, '%s', (type %s) "
                "was passed." % (
                        self._selection_modes.keys(),
                        mode, type(mode)))
        self.score_func = score_func
        self.mode = mode
        self.param = param


    def get_support(self):
        selector = self._selection_modes[self.mode](lambda x:x)
        selector._pvalues = self._pvalues
        selector._scores  = self._scores
        # Now make some acrobaties to set the right named parameter in
        # the selector
        possible_params = selector._get_param_names()
        possible_params.remove('score_func')
        selector._set_params(**{possible_params[0]: self.param})
        return selector.get_support()




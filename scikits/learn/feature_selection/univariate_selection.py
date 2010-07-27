"""
Univariate features selection.
"""

# Authors: V. Michel, B. Thirion, G. Varoquaux, A. Gramfort, E. Duchesnay
# License: BSD 3 clause

import numpy as np
from scipy import stats


######################################################################
# General class for filter univariate selection
######################################################################

class UnivariateFilter(object):
    """
    General class for filter univariate selection
    """

    def __init__(self,ranking,score_func):
        """
        Initialize the univariate feature selection.
        score_func : function taking two arrays X and y, and returning an array.
                     Returns both scores and pvalues.
        ranking : ranking function
        """
        assert callable(score_func), ValueError(
                "The score function should be a callable, '%s' (type %s) "
                "was passed." % (score_func, type(score_func))
            )
        self.score_func = score_func
        self.ranking = ranking

    def fit(self,X,y):
        """
        Evaluate the function
        """
        self._scores = self.score_func(X, y)
        #self._scores = _scores[0]
        #self._pvalues = _scores[1]
        #self._rank = np.argsort(self._pvalues)
        return self

    def transform(self,X,**kwargs):
        """
        Transform a new matrix using the selected features
        """
        self.support = self.ranking.support(self._scores,**kwargs)
        return X[:,self.support]





######################################################################
# Specific rankings
######################################################################

class SelectPercentile(object):
    """
    Filter : Select the best percentile of the p_values
    """
    def __init__(self,percentile):
        self.percentile = percentile

    def support(self,scores,percentile=None):
        if percentile is not None:
                self.percentile = percentile 
        assert self.percentile<=100, ValueError('percentile should be \
                            between 0 and 100 (%f given)' %(self.percentile))
        alpha = stats.scoreatpercentile(scores[1], self.percentile)
        return (scores[1] <= alpha)



class SelectKBest(object):
    """
    Filter : Select the k lowest p-values
    """
    def __init__(self,k):
        self.k = k

    def support(self,scores,k=None):
          if k is not None:
                self.k=k
          assert self.k<=len(scores[1]), ValueError('cannot select %d features'
                                  ' among %d ' % (self.k, len(scores[1])))
          alpha = np.sort(scores[1])[self.k-1]
          return (scores[1] <= alpha)



class SelectFpr(object):
    """
    Filter : Select the pvalues below alpha
    """
    def __init__(self,alpha):
        self.alpha = alpha

    def support(self,alpha = None):
          if alpha is not None:
                self.alpha=alpha
        return (self._pvalues < self.alpha)


class SelectFdr(object):
    """
    Filter : Select the p-values corresponding to an estimated false
    discovery rate of alpha. This uses the Benjamini-Hochberg procedure
    """
    def __init__(self,alpha):
        self.alpha = alpha

    def support(self,alpha = None):
          if alpha is not None:
                self.alpha=alpha
        sv = np.sort(self._pvalues)
        threshold = sv[sv < self.alpha*np.arange(len(self._pvalues))].max()
        return (self._pvalues < threshold)


class SelectFwe(object):
    """
    Filter : Select the p-values corresponding to a corrected p-value of alpha
    """
    def __init__(self,alpha):
        self.alpha = alpha

    def support(self,alpha = None):
          if alpha is not None:
                self.alpha=alpha
        return (self._pvalues < self.alpha/len(self._pvalues))



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
    X = np.asanyarray(X)
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
    y = y.copy()
    X = X.copy()
    if center:
        y -= np.mean(y)
        X -= np.mean(X, 0)

    # compute the correlation
    X /= np.sqrt(np.sum(X**2,0))
    y /= np.sqrt(np.sum(y**2))
    corr = np.dot(y, X)

    # convert to p-value
    dof = y.size-2
    F = corr**2/(1-corr**2)*dof
    pv = stats.f.sf(F, 1, dof)
    return F, pv






if __name__ == "__main__":
    import scikits.learn.datasets.samples_generator as sg
    from scikits.learn.svm import SVR, SVC

    X,y = sg.sparse_uncorrelated(50,100)
    univ_filter = UnivariateFilter(SelectKBest(k=5),f_regression)
    X_r_5 = univ_filter.fit(X, y).transform(X)
    X_r_10 = univ_filter.transform(X,k=10)
    univ_filter.ranking.k = 20
    X_r_20 = univ_filter.fit(X, y).transform(X)
    univ_filter.ranking = SelectPercentile(percentile = 50)
    X_r_50 = univ_filter.fit(X, y).transform(X)

    #clf = SVR(kernel='linear', C=1.)
    #y_ = clf.fit(X_r, y).predict(X_r)
    #print sel

    #### now change k
    #X_r = univariate_filter.transform(X, k=2)
    #y_ = clf.fit(X_r, y).predict(X)
    #print sel



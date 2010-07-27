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

class SelectFpr(UnivariateFilter):
    """
    Filter : Select the pvalues below alpha
    """
    def support(self,alpha):
        return (self._pvalues < alpha)


class SelectFdr(UnivariateFilter):
    """
    Filter : Select the p-values corresponding to an estimated false
    discovery rate of alpha. This uses the Benjamini-Hochberg procedure
    """
    def support(self,alpha):
        sv = np.sort(self._pvalues)
        threshold = sv[sv < alpha*np.arange(len(self._pvalues))].max()
        return (self._pvalues < threshold)


class SelectFwe(UnivariateFilter):
    """
    Filter : Select the p-values corresponding to a corrected p-value of alpha
    """
    def support(self,alpha):
        return (self._pvalues < alpha/len(self._pvalues))



if __name__ == "__main__":
    import scikits.learn.datasets.samples_generator as sg
    from scikits.learn.svm import SVR, SVC
    import scikits.learn.feature_selection.univ_scores as us

    X,y = sg.sparse_uncorrelated(50,100)
    univ_filter = UnivariateFilter(SelectKBest(k=5),us.f_regression)
    X_r_5 = univ_filter.fit(X, y).transform(X)
    X_r_10 = univ_filter.fit(X, y).transform(X,k=10)
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



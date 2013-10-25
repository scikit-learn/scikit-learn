# Author: Andrea Bravi <a.bravi@uottawa.ca>
# License: 3-clause BSD

import numpy as np
from ..base import BaseEstimator
from .base import SelectorMixin
from ..metrics.cluster.supervised import mutual_info_score
from ..utils.validation import array2d


class MinRedundancyMaxRelevance(BaseEstimator, SelectorMixin):
    """
    Select the subset of features with minimal redundancy and maximal
    relevance (mRMR) with the outcome.

    IMPORTANT: This version only supports data in categorical or integer form.

    Attributes
    ----------
    k : int, default=2
        Number of features to select (selected_features)
    mask : list, len=selected_features
           Integer list of the features ordered by maximal relevance and
           minimal redundancy
    score : array, shape=[selected_features]
            mRMR score associated to each entry in mask
    relevance : array, shape=[n_features]
                Relevance of all the features
    redundancy : array, shape=[n_features]
                 Redundancy of all the features
    rule : string, default='diff'
           Rule to combine relevance and redundancy, either
           'diff' - difference between the two
           'prod' - product between the two
    X : array, shape=[n_samples, n_features]
        Input dataset, must be either integer or categorical
    y : array, shape=[n_samples]
        Label vector, must be either integer or categorical

    Methods
    -------
    _compute_mRMR(X, y)
        Computes the minimal relevance maximal redundancy of each feature
        returning mask and score

    References
    ----------
    .. [1] H. Peng, F. Long, and C. Ding, "Feature selection based on mutual
       information: criteria of max-dependency, max-relevance, and
       min-redundancy", IEEE Transactions on Pattern Analysis and Machine
       Intelligence, Vol. 27, No. 8, pp.1226-1238, 2005.
    """
    def __init__(self, k=2, rule='diff'):
        """
        Parameters
        ----------
        k : int, default=2
            Number of features to select
        rule : string, default='diff'
               Rule to combine relevance and redundancy, either
               'diff' - difference between the two
               'prod' - product between the two
        """
        self.k = k
        self.rule = rule
        self._rule_function = self._get_rule_function(rule)

    def fit(self, X, y):
        """
        Parameters
        ----------
        X : array, shape=[n_samples, n_features]
            Input dataset, must be either integer or categorical
        y : array, shape=[n_samples]
            Label vector, must be either integer or categorical
        """
        X = array2d(X)

        self.X = X
        self.y = y
        self.mask, self.score = self._compute_mRMR(X, y)
        return self

    def _get_support_mask(self):
        """
        Returns
        -------
        support : array, dype=bool, shape=[n_features]
                  Boolean mask with True the selected features
        """

        support = np.zeros(self.n_features, dtype=bool)
        support[[self.mask]] = True
        return support

    def _compute_mRMR(self, X, y):
        """
        Parameters
        ----------
        X : array, shape=[n_samples, n_features]
            Input dataset, must be either integer or categorical
        y : array, shape=[n_samples]
            Label vector, must be either integer or categorical

        Returns
        -------
        mask : list, len=selected_features
               Integer list of the features ordered by maximal relevance and
               minimal redundancy
        score : list, len=selected_features
                mRMR score associated to each entry in mask
        """
        M = X.shape[1]  # Number of features
        self.n_features = M

        # Computation of relevance and redundancy
        relevance = np.zeros(M)
        redundancy = np.zeros([M, M])
        for m1 in range(0, M):
            relevance[m1] = mutual_info_score(X[:, m1], y)
            for m2 in range(m1+1, M):
                redundancy[m1, m2] = mutual_info_score(X[:, m1],
                                                       X[:, m2])
                redundancy[m2, m1] = redundancy[m1, m2]

        self.relevance = relevance
        self.redundancy = redundancy

        # Sequential search optimization
        mask = []
        score = []
        search_space = range(0, M)

        score.append(max(relevance))
        ind = int(relevance.argmax(0))  # Optimal feature
        mask.append(ind)
        search_space.pop(ind)

        fun = self._rule_function
        for m in range(0, self.k-1):
            tmp_score = fun(relevance[search_space],
                            np.mean(redundancy[:, search_space].
                                    take(mask, axis=0), 0))
            score.append(max(tmp_score))
            ind = tmp_score.argmax(0)
            mask.append(search_space[ind])
            search_space.pop(ind)

        return mask, score

    def _get_rule_function(self, rule):
        """
        Returns
        -------
        fun : function
              Function used to combine relevance (k) and redundancy (h) arrays
        """
        if rule == 'diff':
            def fun(a, b):
                return a+b
        elif rule == 'prod':
            def fun(a, b):
                return a*b
        else:
            raise ValueError("rule should be either 'diff' or 'prod'")

        return fun

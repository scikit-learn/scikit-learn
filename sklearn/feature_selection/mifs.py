"""
Parallelized Mutual Information based Feature Selection module.
"""

# Author: Daniel Homola <dani.homola@gmail.com>
# License: BSD 3 clause


import numpy as np
from scipy import signal
from ..utils import check_X_y
from ..preprocessing import StandardScaler
import bottleneck as bn
import mi


class MutualInformationFeatureSelector(object):
    """
    MI_FS stands for Mutual Information based Feature Selection.
    This class contains routines for selecting features using both
    continuous and discrete y variables. Three selection algorithms are
    implemented: JMI, JMIM and MRMR.

    This implementation tries to mimic the scikit-learn interface, so use fit,
    transform or fit_transform, to run the feature selection.

    Parameters
    ----------

    method : string, default = 'JMI'
        Which mutual information based feature selection method to use:
        - 'JMI' : Joint Mutual Information [1]
        - 'JMIM' : Joint Mutual Information Maximisation [2]
        - 'MRMR' : Max-Relevance Min-Redundancy [3]

    k : int, default = 5
        Sets the number of samples to use for the kernel density estimation
        with the kNN method. Kraskov et al. recommend a small integer between
        3 and 10.

    n_features : int or string, default = 'auto'
        If int, it sets the number of features that has to be selected from X.
        If 'auto' this is determined automatically based on the amount of
        mutual information the previously selected features share with y.

    categorical : Boolean, default = True
        If True, y is assumed to be a categorical class label. If False, y is
        treated as a continuous. Consequently this parameter determines the
        method of estimation of the MI between the predictors in X and y.

    verbose : int, default=0
        Controls verbosity of output:
        - 0: no output
        - 1: displays selected features
        - 2: displays selected features and mutual information

    Attributes
    ----------

    n_features_ : int
        The number of selected features.

    support_ : array of length X.shape[1]
        The mask array of selected features.

    ranking_ : array of shape n_features
        The feature ranking of the selected features, with the first being
        the first feature selected with largest marginal MI with y, followed by
        the others with decreasing MI.

    mi_ : array of shape n_features
        The JMIM of the selected features. Usually this a monotone decreasing
        array of numbers converging to 0. One can use this to estimate the
        number of features to select. In fact this is what n_features='auto'
        tries to do heuristically.

    Examples
    --------

    import pandas as pd
    import mifs

    # load X and y
    X = pd.read_csv('my_X_table.csv', index_col=0).values
    y = pd.read_csv('my_y_vector.csv', index_col=0).values

    # define MI_FS feature selection method
    feat_selector = mifs.MutualInformationFeatureSelector()

    # find all relevant features
    feat_selector.fit(X, y)

    # check selected features
    feat_selector.support_

    # check ranking of features
    feat_selector.ranking_

    # call transform() on X to filter it down to selected features
    X_filtered = feat_selector.transform(X)

    References
    ----------

    [1] H. Yang and J. Moody, "Data Visualization and Feature Selection: New
        Algorithms for Nongaussian Data"
        NIPS 1999
    [2] Bennasar M., Hicks Y., Setchi R., "Feature selection using Joint Mutual
        Information Maximisation"
        Expert Systems with Applications, Vol. 42, Issue 22, Dec 2015
    [3] H. Peng, Fulmi Long, C. Ding, "Feature selection based on mutual
        information criteria of max-dependency, max-relevance,
        and min-redundancy"
        Pattern Analysis & Machine Intelligence 2005
    """

    def __init__(self, method='JMI', k=5, n_features='auto', categorical=True,
                 verbose=0):
        self.method = method
        self.k = k
        self.n_features = n_features
        self.categorical = categorical
        self.verbose = verbose

    def fit(self, X, y):
        """
        Fits the MI_FS feature selection with the chosen MI_FS method.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values.
        """

        return self._fit(X, y)

    def transform(self, X):
        """
        Reduces the input X to the features selected by chosen MI_FS method.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The training input samples.

        Returns
        -------
        X : array-like, shape = [n_samples, n_features_]
            The input matrix X's columns are reduced to the features which were
            selected by the chosen MI_FS method.
        """

        return self._transform(X)

    def fit_transform(self, X, y):
        """
        Fits MI_FS, then reduces the input X to the selected features.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values.

        Returns
        -------
        X : array-like, shape = [n_samples, n_features_]
            The input matrix X's columns are reduced to the features which were
            selected by the chosen MI_FS method .
        """

        self._fit(X, y)
        return self._transform(X)

    def _fit(self, X, y):
        self.X, y = self._check_params(X, y)
        n, p = X.shape
        self.y = y.reshape((n, 1))

        # list of selected features
        S = []
        # list of all features
        F = range(p)

        if self.n_features != 'auto':
            feature_mi_matrix = np.zeros((self.n_features, p))
        else:
            feature_mi_matrix = np.zeros((n, p))
        feature_mi_matrix[:] = np.nan
        S_mi = []

        # ----------------------------------------------------------------------
        # FIND FIRST FEATURE
        # ----------------------------------------------------------------------

        # check a range of ks (3-10), and choose the one with the max median MI
        k_min = 3
        k_max = 11
        xy_MI = np.zeros((k_max - k_min, p))
        xy_MI[:] = np.nan
        for i, k in enumerate(range(k_min, k_max)):
            xy_MI[i, :] = mi.get_first_mi_vector(self, k)
        xy_MI = bn.nanmedian(xy_MI, axis=0)

        # choose the best, add it to S, remove it from F
        S, F = self._add_remove(S, F, bn.nanargmax(xy_MI))
        S_mi.append(bn.nanmax(xy_MI))

        # notify user
        if self.verbose > 0:
            self._print_results(S, S_mi)

        # ----------------------------------------------------------------------
        # FIND SUBSEQUENT FEATURES
        # ----------------------------------------------------------------------

        while len(S) < self.n_features:
            # loop through the remaining unselected features and calculate MI
            s = len(S) - 1
            feature_mi_matrix[s, F] = mi.get_mi_vector(self, F, s)

            # make decision based on the chosen FS algorithm
            fmm = feature_mi_matrix[:len(S), F]
            if self.method == 'JMI':
                selected = F[bn.nanargmax(bn.nansum(fmm, axis=0))]
            elif self.method == 'JMIM':
                selected = F[bn.nanargmax(bn.nanmin(fmm, axis=0))]
            elif self.method == 'MRMR':
                MRMR = xy_MI[F] - bn.nanmean(fmm, axis=0)
                selected = F[bn.nanargmax(MRMR)]

            # record the JMIM of the newly selected feature and add it to S
            S_mi.append(bn.nanmax(bn.nanmin(fmm, axis=0)))
            S, F = self._add_remove(S, F, selected)

            # notify user
            if self.verbose > 0:
                self._print_results(S, S_mi)

            # if n_features == 'auto', let's check the S_mi to stop
            if self.n_features == 'auto' and len(S) > 10:
                # smooth the 1st derivative of the MI values of previously sel
                MI_dd = signal.savgol_filter(S_mi[1:], 9, 2, 1)
                # does the mean of the last 5 converge to 0?
                if np.abs(np.mean(MI_dd[-5:])) < 1e-3:
                    break

        # ----------------------------------------------------------------------
        # SAVE RESULTS
        # ----------------------------------------------------------------------

        self.n_features_ = len(S)
        self.support_ = np.zeros(p, dtype=np.bool)
        self.support_[S] = 1
        self.ranking_ = S
        self.mi_ = S_mi

        return self

    def _transform(self, X):
        # sanity check
        try:
            self.ranking_
        except AttributeError:
            raise ValueError('You need to call the fit(X, y) method first.')
        X = X[:, self.support_]
        return X

    def _isinteger(self, x):
        return np.all(np.equal(np.mod(x, 1), 0))

    def _check_params(self, X, y):
        # checking input data and scaling it if y is continuous
        X, y = check_X_y(X, y)

        if not self.categorical:
            ss = StandardScaler()
            X = ss.fit_transform(X)
            y = ss.fit_transform(y)

        # sanity checks
        methods = ['JMI', 'JMIM', 'MRMR']
        if self.method not in methods:
            raise ValueError('Please choose one of the following methods:\n' +
                             '\n'.join(methods))

        if not isinstance(self.k, int):
            raise ValueError("k must be an integer.")
        if self.k < 1:
            raise ValueError('k must be larger than 0.')
        if self.categorical and np.any(self.k > np.bincount(y)):
            raise ValueError('k must be smaller than your smallest class.')

        if not isinstance(self.categorical, bool):
            raise ValueError('Categorical must be Boolean.')
        if self.categorical and np.unique(y).shape[0] > 5:
            print 'Are you sure y is categorical? It has more than 5 levels.'
        if not self.categorical and self._isinteger(y):
            print 'Are you sure y is continuous? It seems to be discrete.'
        if self._isinteger(X):
            print ('The values of X seem to be discrete. MI_FS will treat them'
                   'as continuous.')
        return X, y

    def _add_remove(self, S, F, i):
        """
        Helper function: removes ith element from F and adds it to S.
        """

        S.append(i)
        F.remove(i)
        return S, F

    def _print_results(self, S, MIs):
        out = ''
        if self.n_features == 'auto':
            out += 'Auto selected feature #' + str(len(S)) + ' : ' + str(S[-1])
        else:
            out += ('Selected feature #' + str(len(S)) + ' / ' +
                    str(self.n_features) + ' : ' + str(S[-1]))

        if self.verbose > 1:
            out += ', JMIM: ' + str(MIs[-1])
        print out

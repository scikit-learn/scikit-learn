"""
Boruta all-relevant feaure selection algorithm.
"""

# Author: Daniel Homola <dani.homola@gmail.com>
# Original code and method by: Miron B Kursa, https://m2.icm.edu.pl/boruta/
# License: BSD 3 clause


import numpy as np
import scipy as sp
from statsmodels.sandbox.stats.multicomp import multipletests as multicor
from ..utils import check_X_y
from bottleneck import nanrankdata


class BorutaPy(object):
    """
    Improved Python implementation of the Boruta R package.

    The improvements of this implementation include:
    - Faster run times:
        Thanks to scikit-learn's fast implementation of the ensemble methods.
    - Scikit-learn like interface:
        Use BorutaPy just like any other scikit learner: fit, fit_transform and
        transform are all implemented in a similar fashion.
    - Modularity:
        Any ensemble method could be used: random forest, extra trees
        classifier, even gradient boosted trees.
    - Automatic tree number:
        Setting the n_estimator to 'auto' will calculate the number of trees
        in each itartion based on the number of features under investigation.
        This way more trees are used when the training data has many feautres
        and less when most of the features have been rejected.
    - Ranking of features:
        After fitting BorutaPy it provides the user with ranking of features.
        Confirmed ones are 1, Tentatives are 2, and the rejected are ranked
        starting from 3, based on their feautre importance history through
        the iterations.

    We highly recommend using pruned trees with a depth between 3-7.

    For more, see the docs of these functions, and the examples below.

    Original code and method by: Miron B Kursa, https://m2.icm.edu.pl/boruta/

    Boruta is an all relevant feature selection method, while most other are
    minimal optimal; this means it tries to find all features carrying
    information usable for prediction, rather than finding a possibly compact
    subset of features on which some classifier has a minimal error.

    Why bother with all relevant feature selection?
    When you try to understand the phenomenon that made your data, you should
    care about all factors that contribute to it, not just the bluntest signs
    of it in context of your methodology (yes, minimal optimal set of features
    by definition depends on your classifier choice).

    Parameters
    ----------

    estimator : object
        A supervised learning estimator, with a 'fit' method that returns the
        feature_importances_ attribute. Important features must correspond to
        high absolute values in the feature_importances_.

    n_estimators : int or string, default = 1000
        If int sets the number of estimators in the chosen ensemble method.
        If 'auto' this is determined automatically based on the size of the
        dataset. The other parameters of the used estimators need to be set
        with initialisation.

    multi_corr_method : string, default = 'bonferroni'
        Method for correcting for multiple testing during the feature selection
        process. statsmodels' multiple test is used, so one of the following:
        - 'bonferroni' : one-step correction
        - 'sidak' : one-step correction
        - 'holm-sidak' : step down method using Sidak adjustments
        - 'holm' : step-down method using Bonferroni adjustments
        - 'simes-hochberg' : step-up method  (independent)
        - 'hommel' : closed method based on Simes tests (non-negative)
        - 'fdr_bh' : Benjamini/Hochberg  (non-negative)
        - 'fdr_by' : Benjamini/Yekutieli (negative)
        - 'fdr_tsbh' : two stage fdr correction (non-negative)
        - 'fdr_tsbky' : two stage fdr correction (non-negative)

    multi_alpha : float, default = 0.01
        Level at which the corrected p-values will get rejected.

    max_iter : int, default = 100
        The number of maximum iterations to perform.

    verbose : int, default=0
        Controls verbosity of output:
        - 0: no output
        - 1: displays iteration number
        - 2: which features have been selected already

    Attributes
    ----------

    n_features_ : int
        The number of selected features.

    support_ : array of shape [n_features]

        The mask of selected features - only confirmed ones are True.

    support_weak_ : array of shape [n_features]

        The mask of selected tentative features, which haven't gained enough
        support during the max_iter number of iterations..

    ranking_ : array of shape [n_features]

        The feature ranking, such that ``ranking_[i]`` corresponds to the
        ranking position of the i-th feature. Selected (i.e., estimated
        best) features are assigned rank 1 and tentative features are assigned
        rank 2.

    Examples
    --------

    import pandas
    from sklearn.ensemble import RandomForestClassifier
    from boruta_py import boruta_py

    # load X and y
    X = pd.read_csv('my_X_table.csv', index_col=0).values
    y = pd.read_csv('my_y_vector.csv', index_col=0).values

    # define random forest classifier, with utilising all cores and
    # sampling in proportion to y labels
    rf = RandomForestClassifier(n_jobs=-1, class_weight='auto', max_depth=5)

    # define Boruta feature selection method
    feat_selector = boruta_py.BorutaPy(rf, n_estimators='auto', verbose=2)

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

    [1] Kursa M., Rudnicki W., "Feature Selection with the Boruta Package"
        Journal of Statistical Software, Vol. 36, Issue 11, Sep 2010
    """

    def __init__(self, estimator, n_estimators=1000,
                 multi_corr_method='bonferroni', multi_alpha=0.01,
                 max_iter=100, verbose=0):
        self.estimator = estimator
        self.n_estimators = n_estimators
        self.multi_corr_method = multi_corr_method
        self.multi_alpha = multi_alpha
        self.max_iter = max_iter
        self.verbose = verbose

    def fit(self, X, y):
        """
        Fits the Boruta feature selection with the provided estimator.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values.
        """

        return self._fit(X, y)

    def transform(self, X, weak=False):
        """
        Reduces the input X to the features selected by Boruta.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The training input samples.

        weak: boolean, default = False
            If set to true, the tentative features are also used to reduce X.

        Returns
        -------
        X : array-like, shape = [n_samples, n_features_]
            The input matrix X's columns are reduced to the features which were
            selected by Boruta.
        """

        return self._transform(X, weak)

    def fit_transform(self, X, y, weak=False):
        """
        Fits Boruta, then reduces the input X to the selected features.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values.

        weak: boolean, default = False
            If set to true, the tentative features are also used to reduce X.

        Returns
        -------
        X : array-like, shape = [n_samples, n_features_]
            The input matrix X's columns are reduced to the features which were
            selected by Boruta.
        """

        self._fit(X, y)
        return self._transform(X, weak)

    def _fit(self, X, y):
        # check input params
        self._check_params(X, y)

        # setup variables for Boruta
        n_sample, n_feat = X.shape
        iter = 1
        # holds the decision about each feature:
        # 0  - default state = tentative in original code
        # 1  - accepted in original code
        # -1 - rejected in original code
        dec_reg = np.zeros(n_feat, dtype=np.int)
        # counts how many times a given feature was more important than
        # the best of the shadow features
        hit_reg = np.zeros(n_feat, dtype=np.int)
        # these record the history of the iterations
        imp_history = np.zeros(n_feat, dtype=np.float)
        sha_max_history = []

        # set n_estimators
        if self.n_estimators != 'auto':
            self.estimator.set_params(n_estimators=self.n_estimators)

        # main feature selection loop
        while np.any(dec_reg == 0) and iter < self.max_iter:
            # find optimal number of trees
            if self.n_estimators == 'auto':
                not_rejected = np.where(dec_reg >= 0)[0]
                sample_size = np.min([n_sample, not_rejected.shape[0]])
                n_tree = self._get_tree_num(sample_size)
                self.estimator.set_params(n_estimators=n_tree)

            # make sure we start with a new tree in each iteration
            rnd_st = np.random.randint(1, 1e6, 1)[0]
            self.estimator.set_params(random_state=rnd_st)

            # add shadow attributes, shuffle them and train estimator, get imps
            cur_imp = self._add_shadows_get_imps(X, y, dec_reg)

            # record importance history
            imp_sha_max = np.max(cur_imp[1])
            sha_max_history.append(imp_sha_max)
            imp_history = np.vstack((imp_history, cur_imp[0]))

            # register which feature is more imp than the max of shadows
            hit_reg = self._assign_hits(hit_reg, cur_imp, imp_sha_max)

            # based on hit_reg we check if a feature is doing better than
            # expected by chance
            dec_reg = self._do_tests(dec_reg, hit_reg, iter)

            # print out confirmed features
            if self.verbose > 0 and iter < self.max_iter:
                self._print_results(dec_reg, iter, 0)
            if iter < self.max_iter:
                iter += 1

        # we automatically apply R package's rough fix for tentative ones
        confirmed = np.where(dec_reg == 1)[0]
        tentative = np.where(dec_reg == 0)[0]
        # ignore the first row of zeros
        tentative_median = np.median(imp_history[1:, tentative], axis=0)
        # which tentative to keep
        tentative_confirmed = np.where(tentative_median
                                       > np.median(sha_max_history))[0]
        tentative = tentative[tentative_confirmed]

        # basic result variables
        self.n_features_ = confirmed.shape[0]
        self.support_ = np.zeros(n_feat, dtype=np.bool)
        self.support_[confirmed] = 1
        self.support_weak_ = np.zeros(n_feat, dtype=np.bool)
        self.support_weak_[tentative] = 1

        # ranking, confirmed variables are rank 1
        self.ranking_ = np.ones(n_feat, dtype=np.int)
        # tentative variables are rank 2
        self.ranking_[tentative] = 2
        # selected = confirmed and tentative
        selected = np.hstack((confirmed, tentative))
        # all rejected features are sorted by importance history
        not_selected = np.setdiff1d(np.arange(n_feat), selected)
        # large importance values should rank higher = lower ranks -> *(-1)
        imp_history_rejected = imp_history[1:, not_selected] * -1
        # calculate ranks in each iteration, then median of ranks across feats
        iter_ranks = nanrankdata(imp_history_rejected, axis=1)
        rank_medians = np.nanmedian(iter_ranks, axis=0)
        ranks = nanrankdata(rank_medians)
        # set smallest rank to 3 if there are tentative feats
        if tentative.shape[0] > 0:
            ranks = ranks - np.min(ranks) + 3
        else:
            # and 2 otherwise
            ranks = ranks - np.min(ranks) + 2
        self.ranking_[not_selected] = ranks

        # notify user
        if self.verbose > 0:
            self._print_results(dec_reg, iter, 1)
        return self

    def _transform(self, X, weak=False):
        # sanity check
        try:
            self.ranking_
        except AttributeError:
            raise ValueError('You need to call the fit(X, y) method first.')

        if weak:
            X = X[:, self.support_ + self.support_weak_]
        else:
            X = X[:, self.support_]
        return X

    def _check_params(self, X, y):
        X, y = check_X_y(X, y)
        multi_corr_methods = ['bonferroni', 'sidak', 'holm-sidak', 'holm',
                              'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by',
                              'fdr_tsbh', 'fdr_tsbky']
        if self.multi_corr_method not in multi_corr_methods:
            raise ValueError('For multiple testing correction method, please '
                             'choose one of the following:\n' +
                             '\n'.join(multi_corr_methods))
        if self.multi_alpha <= 0 or self.multi_alpha > 1:
            raise ValueError('Multi_alpha should be between 0 and 1.')

    def _print_results(self, dec_reg, iter, flag):
        n_iter = str(iter) + ' / ' + str(self.max_iter)
        n_confirmed = np.where(dec_reg == 1)[0].shape[0]
        n_rejected = np.where(dec_reg == -1)[0].shape[0]
        cols = ['Iteration: ', 'Confirmed: ', 'Tentative: ', 'Rejected: ']

        # still in feature selection
        if flag == 0:
            n_tentative = np.where(dec_reg == 0)[0].shape[0]
            content = map(str, [n_iter, n_confirmed, n_tentative, n_rejected])
            if self.verbose == 1:
                output = cols[0] + n_iter
            elif self.verbose > 1:
                output = '\n'.join([x[0] + '\t' + x[1]
                                    for x in zip(cols, content)])

        # Boruta finished running and tentatives have been filtered
        else:
            n_tentative = np.sum(self.support_weak_)
            content = map(str, [n_iter, n_confirmed, n_tentative, n_rejected])
            result = '\n'.join([x[0] + '\t' + x[1]
                                for x in zip(cols, content)])
            output = "\n\nBorutaPy finished running.\n\n" + result
        print output

    def _get_tree_num(self, n_feat):
        depth = self.estimator.get_params()['max_depth']
        if depth is None:
            depth = 10
        # how many times a feature should be considered on average
        f_repr = 100
        # 2 because the training matrix is extended with n shadow features
        multi = ((n_feat * 2) / float(np.sqrt(n_feat * 2) * depth))
        n_estimators = int(multi * f_repr)
        return int(n_estimators)

    def _get_imp(self, X, y):
        try:
            self.estimator.fit(X, y)
        except Exception as e:
            raise ValueError('Please check your X and y variable. The provided'
                             'estimator cannot be fitted to your data.\n' + e)
        try:
            imp = self.estimator.feature_importances_
        except Exception:
            raise ValueError('Only methods with feature_importance_ attribute '
                             'are currently supported in BorutaPy.')
        return imp

    def _get_shuffle(self, seq):
        np.random.shuffle(seq)
        return seq

    def _add_shadows_get_imps(self, X, y, dec_reg):
        # find features that are tentative still
        x_cur_ind = np.where(dec_reg >= 0)[0]
        x_cur = np.copy(X[:, x_cur_ind])
        x_cur_w = x_cur.shape[1]

        # deep copy the matrix for the shadow matrix
        x_sha = np.copy(x_cur)
        # make sure there's at least 5 columns in the shadow matrix for
        while (x_sha.shape[1] < 5):
            x_sha = np.hstack((x_sha, x_sha))
        # shuffle xSha
        x_sha = np.apply_along_axis(self._get_shuffle, 0, x_sha)
        # get importance of the merged matrix
        imp = self._get_imp(np.hstack((x_cur, x_sha)), y)

        # separate importances of real and shadow features
        imp_sha = imp[x_cur_w:]
        imp_real = np.zeros(X.shape[1])
        imp_real[:] = np.nan
        imp_real[x_cur_ind] = imp[:x_cur_w]
        return (imp_real, imp_sha)

    def _assign_hits(self, hit_reg, cur_imp, imp_sha_max):
        # register hits for feautres that did better than the best of shadows
        hits = np.where(cur_imp[0] > imp_sha_max)[0]
        hit_reg[hits] += 1
        return hit_reg

    def _do_tests(self, dec_reg, hit_reg, iter):
        # get uncorrected p values based on hit_reg
        to_accept_ps = sp.stats.binom.sf(hit_reg - 1, iter, .5).flatten()
        to_reject_ps = sp.stats.binom.cdf(hit_reg, iter, .5).flatten()

        # correct p values for multiple testing
        to_accept = np.where(multicor(to_accept_ps, alpha=self.multi_alpha,
                                      method=self.multi_corr_method)[0])[0]
        to_reject = np.where(multicor(to_reject_ps, alpha=self.multi_alpha,
                                      method=self.multi_corr_method)[0])[0]

        # finding those to_accept and to_reject that are 0, and setting them
        dec_reg[to_accept[np.where(dec_reg[to_accept] == 0)]] = 1
        dec_reg[to_reject[np.where(dec_reg[to_reject] == 0)]] = -1
        return dec_reg

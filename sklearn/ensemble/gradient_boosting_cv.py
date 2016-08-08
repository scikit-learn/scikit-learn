"""
This module contains Gradient Boosting Classifier and Regressor with
Cross Validation.
"""

# Authors: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Vighnesh Birodkar <vighneshbirodkar@nyu.edu>
#          Raghav RV <rvraghav93@gmail.com>
# License: BSD 3 clause

import warnings
import numpy as np

from ..base import BaseEstimator
from ..metrics.scorer import check_scoring
from ..utils.validation import check_X_y

from .gradient_boosting import GradientBoostingClassifier
from .gradient_boosting import GradientBoostingRegressor
from ..tree._tree import DTYPE

from ..model_selection import check_cv
from ..model_selection import ParameterGrid
from ..model_selection import ParameterSampler
from ..model_selection._search import _CVScoreTuple

from ..externals.joblib import Parallel, delayed

__all__ = ["GradientBoostingClassifierCV", "GradientBoostingRegressorCV"]


class GradientBoostingClassifierCV(BaseEstimator):
    """Gradient Boosting Classifier with Cross Validation

    This class implements GradientBoostingClassifier with an additional support
    to select the best hyper-parameter values from the list of given values.

    For each parameter setting, the best number of stages (``n_estimators``) is
    found by iterating until the last ``n_iter_no_change`` stages show no
    improvement on the previous scores when evaluated on the validation set.

    This approach will be significantly faster compared to
    :class:`sklearn.model_selection.GridSearchCV` due to early stopping.
    All parameters which are not parsed by this class are passed directly to
    the underlying Gradient Boosting estimator.

    Read more in the :ref:`User Guide <gradient_boosting>`.

    Parameters
    ----------
    n_stop_rounds : int, optional, default=10
        If the score on the test set rounded off to ``score_precision`` decimal
        places does not change for ``n_stop_rounds`` iterations, the gradient
        boosting is halted.

        Set this value to -1 to disable early stopping.

    score_precision : int, optional, default=2
        The number of decimal places to round the score off by before
        comparing it with the scores of previous iterations.

    max_iter : int, optinal, default=10000
        The maximum number of estimators that will be added to the Gradient
        Boosting class. If the model converges, the number of estimators will
        be lower than this.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
          - None, to use the default 3-fold cross validation,
          - integer, to specify the number of folds in a ``(Stratified)KFold``,
          - An object to be used as a cross-validation generator.
          - An iterable yielding train, validation splits.

        For integer/None inputs, ``StratifiedKFold`` is used for classification
        tasks, when ``y`` is binary or multiclass.

        See the :mod:`sklearn.model_selection` module for the list of
        cross-validation strategies that can be used here.

        Also refer :ref:`cross-validation documentation <cross_validation>`

    scoring : string, callable or None, default=None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.
        If ``None``, the ``score`` method of the estimator is used.

    refit : boolean, default=True
        If set, the data will be fit on a new Gradient Boosting model with
        the parameters that gave the best mean score. If not set, it is not
        possbile to make predections later with this class.

    n_jobs : int, default=1
        Number of jobs to run in parallel.

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    loss : string, list of strings, optional (default='deviance')
        Loss function to be optimized. 'deviance' refers to
        deviance (= logistic regression) for classification
        with probabilistic outputs. For loss 'exponential' gradient
        boosting recovers the AdaBoost algorithm. If a list is given, all
        loss functions are evaluated to choose the best one.

    learning_rate : float, list of floats, optional (default=0.1)
        Learning rate shrinks the contribution of each tree by
        ``learning_rate``.
        There is a trade-off between ``learning_rate`` and ``n_estimators``.
        If a list is given, all learning rates will be evaluated.

    max_depth : integer, list of integers, optional (default=3)
        Maximum depth of the individual regression estimators. The maximum
        depth limits the number of nodes in the tree. Tune this parameter
        for best performance; the best value depends on the interaction
        of the input variables.
        Ignored if ``max_leaf_nodes`` is not None.
        If a list is given, each item in the list will be evaluated.

    min_samples_split : int, float, list of int or float, optional (default=2)
        The minimum number of samples required to split an internal node:
        - If int, then consider ``min_samples_split`` as the minimum number.
        - If float, then ``min_samples_split`` is a percentage and
          ``ceil(min_samples_split * n_samples)`` are the minimum
          number of samples for each split.
        If a list is given, each item in the list will be evaluated.

    min_samples_leaf : int, float, list of int or floats, optional (default=1)
        The minimum number of samples required to be at a leaf node:
        - If int, then consider ``min_samples_leaf`` as the minimum number.
        - If float, then ``min_samples_leaf`` is a percentage and
          ``ceil(min_samples_leaf * n_samples)`` are the minimum
          number of samples for each node.
        If a list is given, each item in the list will be evaluated.

    min_weight_fraction_leaf : float, list of floats, optional (default=0.)
        The minimum weighted fraction of the input samples required to be at a
        leaf node. If a list is given, each item in the list will be evaluated.

    subsample : float, list of floats, optional (default=1.0)
        The fraction of samples to be used for fitting the individual base
        learners. If smaller than 1.0 this results in Stochastic Gradient
        Boosting. ``subsample`` interacts with the parameter ``n_estimators``.
        Choosing ``subsample < 1.0`` leads to a reduction of variance
        and an increase in bias. If a list is given, each item in the list
        will be evaluated.

    max_features : int, float, string or None, list, optional (default=None)
        The number of features to consider when looking for the best split:
        - If int, then consider ``max_features`` features at each split.
        - If float, then ``max_features`` is a percentage and
          ``int(max_features * n_features)`` features are considered at each
          split.
        - If ``"auto"``, then ``max_features=sqrt(n_features)``.
        - If ``"sqrt"``, then ``max_features=sqrt(n_features)``.
        - If ``"log2"``, then ``max_features=log2(n_features)``.
        - If ``None``, then ``max_features=n_features``.
        Choosing ``max_features < n_features`` leads to a reduction of variance
        and an increase in bias.
        Note: the search for a split does not stop until at least one
        valid partition of the node samples is found, even if it requires to
        effectively inspect more than ``max_features`` features. If a list is
        given, each item in the list will be evaluated.

    max_leaf_nodes : int, list of int, or None, optional (default=None)
        Grow trees with ``max_leaf_nodes`` in best-first fashion.
        Best nodes are defined as relative reduction in impurity.
        If None then unlimited number of leaf nodes.
        If not None then ``max_depth`` will be ignored.
        If a list is given, each item in the list will be evaluated.

    init : BaseEstimator, None, list, optional (default=None)
        An estimator object that is used to compute the initial
        predictions. ``init`` has to provide ``fit`` and ``predict``.
        If None it uses ``loss.init_estimator``.
        If a list is given, each item in the list will be evaluated.

    random_state : int, RandomState instance or None,
                   list of int or RandomState, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by ``np.random``.

    presort : bool or 'auto', list of string or bool optional (default='auto')
        Whether to presort the data to speed up the finding of best splits in
        fitting. Auto mode by default will use presorting on dense data and
        default to normal sorting on sparse data. Setting presort to true on
        sparse data will raise an error. If a list is given, each item in the
        list will be evaluated.

    Attributes
    ----------
    grid_scores_ : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.
        Each named tuple has the attributes:

        - ``parameters``, a dict of parameter settings
        - ``mean_validation_score``, the mean score over the
          cross-validation folds
        - ``cv_validation_scores``, the list of scores for each fold

    best_score_tuple_ : named tuple
        The attributes of the model which gave the best score.

        - ``parameters``, a dict of parameter settings, which includes the
          ``n_estimators`` parameter
        - ``mean_validation_score``, the mean score over the
          cross-validation folds
        - ``cv_validation_scores``, the list of scores for each fold

    best_params_ : dict
        The parameters which gave the best score, along with the
        ``n_estimators`` parameter.

    best_estimator_ : GradientBoostingClassifier
        If ``refit`` was ``True``, the model with the best score fit to the
        entire dataset.

    """

    _estimator_class = GradientBoostingClassifier

    def __init__(self, n_stop_rounds=10, score_precision=2,
                 max_iterations=10000, cv=3, scoring=None, refit=True,
                 n_jobs=1, pre_dispatch='2*n_jobs', verbose=0, loss='deviance',
                 learning_rate=0.1, subsample=1.0, min_samples_split=2,
                 min_samples_leaf=1, min_weight_fraction_leaf=0., max_depth=3,
                 init=None, random_state=None, max_features=None,
                 max_leaf_nodes=None, presort='auto'):

        self.n_stop_rounds = n_stop_rounds
        self.score_precision = score_precision
        self.max_iterations = max_iterations
        self.cv = cv
        self.scoring = scoring
        self.refit = refit
        self.n_jobs = n_jobs
        self.pre_dispatch = pre_dispatch
        self.verbose = verbose

        self.loss = loss
        self.learning_rate = learning_rate
        self.subsample = subsample
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_fraction_leaf = min_weight_fraction_leaf
        self.max_depth = max_depth
        self.init = init
        self.random_state = random_state
        self.max_features = max_features
        self.max_leaf_nodes = max_leaf_nodes
        self.presort = presort

    def fit(self, X, y):
        """Run fit with all sets of parameters till convergence.

        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.
        y : array-like, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)
            For classification, labels must correspond to classes.
        """

        params = self._get_params()
        X, y = check_X_y(X, y, dtype=DTYPE)

        cv = check_cv(self.cv, y)
        param_iter = ParameterGrid(params)

        parallel = Parallel(n_jobs=self.n_jobs, verbose=self.verbose,
                            backend='threading')

        out = parallel(delayed(_fit_single_param)
                       (self._estimator_class, X, y, train, validation, params,
                        self.n_stop_rounds, self.max_iterations, self.scoring,
                        self.score_precision, self.random_state)
                       for train, validation in cv.split(X)
                       for params in param_iter)

        n_splits = int(len(out) / len(param_iter))

        grid_scores = []
        i = 0
        best_mean = -np.inf
        for params in param_iter:

            score_list = []
            n_est_list = []
            for idx in range(n_splits):
                score, n_estimators = out[i]
                n_est_list.append(n_estimators)
                score_list.append(score)
                i += 1

            scores = np.array(score_list)
            mean = np.mean(scores)
            params['n_estimators'] = int(np.mean(n_est_list))
            score_tuple = _CVScoreTuple(params, mean, scores)
            grid_scores.append(_CVScoreTuple(params, mean, scores))

            if mean > best_mean:
                best_tuple = score_tuple
                best_params = params

        self.grid_scores_ = grid_scores
        self.best_score_tuple_ = best_tuple
        self.best_params_ = best_params

        if self.refit:
            self.best_params_['random_state'] = self.random_state
            gb = self._estimator_class(**self.best_params_)
            gb.fit(X, y)
            self.best_estimator_ = gb

        return self

    def predict(self, X):
        """Call ``predict`` on the best estimator.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.
        Returns
        -------
        y: array of shape = ["n_samples]
            The predicted values.
        """
        return self.best_estimator_.predict(X)

    def predict_proba(self, X):
        """Call ``predict_proba`` on the best estimator.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Raises
        ------
        AttributeError
            If the ``loss`` does not support probabilities.

        Returns
        -------
        p : array of shape = [n_samples]
            The class probabilities of the input samples. The order of the
            classes corresponds to that in the attribute ``classes_``.
        """
        return self.best_estimator_.predict_proba(X)

    def predict_log_proba(self, X):
        """Call ``predict_log_proba`` on the best estimator.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.
        Raises
        ------

        AttributeError
            If the ``loss`` does not support probabilities.

        Returns
        -------
        p : array of shape = [n_samples]
            The class probabilities of the input samples. The order of the
            classes corresponds to that in the attribute ``classes_``.
        """
        return self.best_estimator_.predict_log_proba(X)

    def _get_params(self):
        "Return parametrs to iterate over for underlying gradient boosting"

        params = {
            'loss': self.loss,
            'learning_rate': self.learning_rate,
            'subsample': self.subsample,
            'min_samples_split': self.min_samples_split,
            'min_samples_leaf': self.min_samples_leaf,
            'min_weight_fraction_leaf': self.min_weight_fraction_leaf,
            'max_depth': self.max_depth,
            'init': self.init,
            'max_features': self.max_features,
            'max_leaf_nodes': self.max_leaf_nodes,
            'presort': self.presort,
        }

        # Pre processing to ensure every parameter is a list
        # which ``ParameterGrid`` can iterate over
        for key, value in params.items():
            params[key] = np.atleast_1d(value)

        return params


class GradientBoostingRegressorCV(BaseEstimator):
    """Gradient Boosting Regressor with Cross Validation

    This class implements GradientBoostingRegressor with an additional support
    to select the best hyper-parameter values from the list of given values.

    For each parameter setting, the best number of stages (``n_estimators``) is
    found by iterating until the last ``n_iter_no_change`` stages show no
    improvement on the previous scores when evaluated on the validation set.

    This approach will be significantly faster compared to
    :class:`sklearn.model_selection.GridSearchCV` due to early stopping.
    All parameters which are not parsed by this class are passed directly to
    the underlying Gradient Boosting estimator.

    Read more in the :ref:`User Guide <gradient_boosting>`.

    Parameters
    ----------
    n_stop_rounds : int, optional, default=10
        If the score on the test set rounded off to ``score_precision`` decimal
        places does not change for ``n_stop_rounds`` iterations, the gradient
        boosting is halted.

        Set this value to -1 to disable early stopping.

    score_precision : int, optional, default=2
        The number of decimal places to round the score off by before
        comparing it with the scores of previous iterations.

    max_iter : int, optinal, default=10000
        The maximum number of estimators that will be added to the Gradient
        Boosting class. If the model converges, the number of estimators will
        be lower than this.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
          - ``None``, to use the default 3-fold cross validation,
          - integer, to specify the number of folds in a ``(Stratified)KFold``,
          - An object to be used as a cross-validation generator.
          - An iterable yielding train, validation splits.

        For integer/``None`` inputs, ``StratifiedKFold`` is used for
        classification tasks, when ``y`` is binary or multiclass.

        See the :mod:`sklearn.model_selection` module for the list of
        cross-validation strategies that can be used here.

        Also refer :ref:`cross-validation documentation <cross_validation>`

    scoring : string, callable or None, default=None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.
        If ``None``, the ``score`` method of the estimator is used.

    refit : boolean, default=True
        If set, the data will be fit on a new Gradient Boosting model with
        the parameters that gave the best mean score. If not set, it is not
        possbile to make predections later with this class.

    n_jobs : int, default=1
        Number of jobs to run in parallel.

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    loss : {'ls', 'lad', 'huber', 'quantile'}, optional (default='ls')
        loss function to be optimized. ``'ls'`` refers to least squares
        regression. ``'lad'`` (least absolute deviation) is a highly robust
        loss function solely based on order information of the input
        variables. ``'huber'`` is a combination of the two. ``'quantile'``
        allows quantile regression (use ``alpha`` to specify the quantile).

    learning_rate : float, list of floats, optional (default=0.1)
        learning rate shrinks the contribution of each tree by
        ``learning_rate``. There is a trade-off between ``learning_rate`` and
        ``n_estimators``. If a list is given, all learning rates will be
        evaluated.

    max_depth : integer, list of integers, optional (default=3)
        maximum depth of the individual regression estimators. The maximum
        depth limits the number of nodes in the tree. Tune this parameter
        for best performance; the best value depends on the interaction
        of the input variables.
        Ignored if ``max_leaf_nodes`` is not ``None``.
        If a list is given, each item in the list will be evaluated.

    min_samples_split : int, float, list of int or float, optional (default=2)
        The minimum number of samples required to split an internal node:
        - If int, then consider ``min_samples_split`` as the minimum number.
        - If float, then ``min_samples_split`` is a percentage and
          ``ceil(min_samples_split * n_samples)`` are the minimum
          number of samples for each split.
        If a list is given, each item in the list will be evaluated.

    min_samples_leaf : int, float, list of int or floats, optional (default=1)
        The minimum number of samples required to be at a leaf node:
        - If int, then consider ``min_samples_leaf`` as the minimum number.
        - If float, then ``min_samples_leaf`` is a percentage and
          ``ceil(min_samples_leaf * n_samples)`` are the minimum
          number of samples for each node.
        If a list is given, each item in the list will be evaluated.

    min_weight_fraction_leaf : float, list of floats, optional (default=0.)
        The minimum weighted fraction of the input samples required to be at a
        leaf node. If a list is given, each item in the list will be evaluated.

    subsample : float, list of floats, optional (default=1.0)
        The fraction of samples to be used for fitting the individual base
        learners. If smaller than 1.0 this results in Stochastic Gradient
        Boosting. ``subsample`` interacts with the parameter ``n_estimators``.
        Choosing ``subsample < 1.0`` leads to a reduction of variance
        and an increase in bias. If a list is given, each item in the list
        will be evaluated.

    max_features : int, float, string or None, list, optional (default=None)
        The number of features to consider when looking for the best split:
        - If int, then consider ``max_features`` features at each split.
        - If float, then ``max_features`` is a percentage and
          ``int(max_features * n_features)`` features are considered at each
          split.
        - If "auto", then ``max_features=sqrt(n_features)``.
        - If "sqrt", then ``max_features=sqrt(n_features)``.
        - If "log2", then ``max_features=log2(n_features)``.
        - If None, then ``max_features=n_features``.
        Choosing ``max_features < n_features`` leads to a reduction of variance
        and an increase in bias.
        Note: the search for a split does not stop until at least one
        valid partition of the node samples is found, even if it requires to
        effectively inspect more than ``max_features`` features. If a list is
        given, each item in the list will be evaluated.

    max_leaf_nodes : int, list of int, or None, optional (default=None)
        Grow trees with ``max_leaf_nodes`` in best-first fashion.
        Best nodes are defined as relative reduction in impurity.
        If None then unlimited number of leaf nodes.
        If not None then ``max_depth`` will be ignored.
        If a list is given, each item in the list will be evaluated.

    init : BaseEstimator, None, list, optional (default=None)
        An estimator object that is used to compute the initial
        predictions. ``init`` has to provide ``fit`` and ``predict``.
        If None it uses ``loss.init_estimator``.
        If a list is given, each item in the list will be evaluated.

    random_state : int, RandomState instance or None,
                   list of int or RandomState, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by ``np.random``.

    presort : bool or 'auto', list of string or bool optional (default='auto')
        Whether to presort the data to speed up the finding of best splits in
        fitting. Auto mode by default will use presorting on dense data and
        default to normal sorting on sparse data. Setting presort to true on
        sparse data will raise an error. If a list is given, each item in the
        list will be evaluated.

    Attributes
    ----------
    grid_scores_ : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.
        Each named tuple has the attributes:

            * ``parameters``, a dict of parameter settings
            * ``mean_validation_score``, the mean score over the
              cross-validation folds
            * ``cv_validation_scores``, the list of scores for each fold

    best_score_tuple_ : named tuple
        The attributes of the model which gave the best score.

        * ``parameters``, a dict of parameter settings, which includes the
          ``n_estimators`` parameter
        * ``mean_validation_score``, the mean score over the
          cross-validation folds
        * ``cv_validation_scores``, the list of scores for each fold

    best_params_ : dict
        The parameters which gave the best score, along with the
        ``n_estimators`` parameter.

    best_estimator_ : GradientBoostingClassifier
        If ``refit`` was ``True``, the model with the best score fit to the
        entire dataset.

    """

    _estimator_class = GradientBoostingRegressor

    def __init__(self, n_stop_rounds=10, score_precision=2,
                 max_iterations=10000, cv=3, scoring=None, refit=True,
                 n_jobs=1, pre_dispatch='2*n_jobs', verbose=0, loss='ls',
                 learning_rate=0.1, subsample=1.0, min_samples_split=2,
                 min_samples_leaf=1, min_weight_fraction_leaf=0., max_depth=3,
                 init=None, random_state=None, max_features=None,
                 max_leaf_nodes=None, presort='auto'):

        self.n_stop_rounds = n_stop_rounds
        self.score_precision = score_precision
        self.max_iterations = max_iterations
        self.cv = cv
        self.scoring = scoring
        self.refit = refit
        self.n_jobs = n_jobs
        self.pre_dispatch = pre_dispatch
        self.verbose = verbose

        self.loss = loss
        self.learning_rate = learning_rate
        self.subsample = subsample
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_fraction_leaf = min_weight_fraction_leaf
        self.max_depth = max_depth
        self.init = init
        self.random_state = random_state
        self.max_features = max_features
        self.max_leaf_nodes = max_leaf_nodes
        self.presort = presort

    def fit(self, X, y):
        """Run fit with all sets of parameters till convergence.

        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.
        y : array-like, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)
            For classification, labels must correspond to classes.
        """

        params = self._get_params()
        X, y = check_X_y(X, y, dtype=DTYPE)

        cv = check_cv(self.cv, y)
        param_iter = ParameterGrid(params)

        parallel = Parallel(n_jobs=self.n_jobs, verbose=self.verbose,
                            backend='threading')

        out = parallel(delayed(_fit_single_param)
                       (self._estimator_class, X, y, train, validation, params,
                        self.n_stop_rounds, self.max_iterations, self.scoring,
                        self.score_precision, self.random_state)
                       for train, validation in cv.split(X)
                       for params in param_iter)

        n_splits = int(len(out) / len(param_iter))

        grid_scores = []
        i = 0
        best_mean = -np.inf
        for params in param_iter:

            score_list = []
            n_est_list = []
            for idx in range(n_splits):
                score, n_estimators = out[i]
                n_est_list.append(n_estimators)
                score_list.append(score)
                i += 1

            scores = np.array(score_list)
            mean = np.mean(scores)
            params['n_estimators'] = int(np.mean(n_est_list))
            score_tuple = _CVScoreTuple(params, mean, scores)
            grid_scores.append(_CVScoreTuple(params, mean, scores))

            if mean > best_mean:
                best_tuple = score_tuple
                best_params = params

        self.grid_scores_ = grid_scores
        self.best_score_tuple_ = best_tuple
        self.best_params_ = best_params

        if self.refit:
            self.best_params_['random_state'] = self.random_state
            gb = self._estimator_class(**self.best_params_)
            gb.fit(X, y)
            self.best_estimator_ = gb

        return self

    def predict(self, X):
        """Call ``predict`` on the best estimator.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.
        Returns
        -------
        y: array of shape = ["n_samples]
            The predicted values.
        """
        return self.best_estimator_.predict(X)

    def _get_params(self):
        "Return parametrs to iterate over for underlying gradient boosting"

        params = {
            'loss': self.loss,
            'learning_rate': self.learning_rate,
            'subsample': self.subsample,
            'min_samples_split': self.min_samples_split,
            'min_samples_leaf': self.min_samples_leaf,
            'min_weight_fraction_leaf': self.min_weight_fraction_leaf,
            'max_depth': self.max_depth,
            'init': self.init,
            'max_features': self.max_features,
            'max_leaf_nodes': self.max_leaf_nodes,
            'presort': self.presort,
        }

        # Pre processing to ensure every parameter is a list
        # which ``ParameterGrid`` can iterate over
        for key, value in params.items():
            params[key] = np.atleast_1d(value)

        return params


def _fit_single_param(estimator, X, y, train, validation, params, stop_rounds,
                      max_iter, scoring, score_precision, random_state):
    """ Fit a single estimator till stopping criteria is reached.

    Parameters
    ---------
    estimator : GradientBoostingRegressor or GradientBoostingClassifier
        The class of the model to instantiate.

    X : array-like, shape = [n_samples, n_features]
        Training vectors, where n_samples is the number of samples
        and n_features is the number of features.

    y : array-like, shape = [n_samples]
        Target values (integers in classification, real numbers in
        regression)
        For classification, labels must correspond to classes.

    train : array-like
        The indices of the training samples.

    validation : array-like
        The indices of the validation samples.

    params : dict
        The parameters to instantiate the Gradient Boosting model with

    stop_rounds : int
        The number of iterations for which the score should remain constant
        over to stop training.

        Set to -1 to disable early stopping.

    max_iter : int
        The maximum number of estimators to fit during gradient boosting.

    scoring : string, callable or None, default=None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.
        If ``None``, the ``score`` method of the estimator is used.

    score_precision : int
        The number of decimal places the score during each iterations is
        rounded off to before comparing with the previous ones.

    random_state : int, RandomState instance or None,
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by ``np.random``.


    Returns
    -------
    score: float
        The score when the model converged.

    n_iterations:
        The number of iterations it took to converge. This is the same as
        parameter ``n_estimators`` for the gradient boosting model
    """

    params['random_state'] = random_state
    params['warm_start'] = True
    gb = estimator(**params)

    scorer = check_scoring(estimator, scoring=scoring)
    scores = []
    rounded_scores = []

    X_train = X[train]
    y_train = y[train]
    X_validation = X[validation]
    y_validation = y[validation]

    for i in range(1, max_iter + 1):
        gb.n_estimators = i
        gb.fit(X_train, y_train)

        this_score = scorer(gb, X_validation, y_validation)
        this_score_rounded = np.round(this_score, score_precision)

        scores.append(this_score)
        rounded_scores.append(this_score_rounded)

        # Check if we need to stop early
        if (stop_rounds > 0 and i >= stop_rounds and
                np.all(this_score_rounded <= rounded_scores[-stop_rounds:])):
            break

    if i == max_iter:
        warnings.warn(str(gb) + ' failed to converge')

    return scores[-1], i

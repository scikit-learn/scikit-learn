"""
The :mod:`sklearn.model_selection.validate` module includes
"""
#TODO Complete docstring

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>,
#         Gael Varoquaux <gael.varoquaux@normalesup.org>,
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Alexander Fabisch <afabisch@informatik.uni-bremen.de>
# License: BSD 3 clause

from __future__ import print_function
from __future__ import division

import numbers
import time
import warnings

import numpy as np

from sklearn.base import is_classifier, clone
from sklearn.utils import check_arrays, check_random_state, safe_mask
from sklearn.utils.validation import _num_samples
from sklearn.externals.joblib import Parallel, delayed, logger
from sklearn.metrics.scorer import check_scoring
from sklearn.utils import check_arrays
from sklearn.utils.fixes import astype
from .partition import _check_cv

def cross_val_score(estimator, X, y=None, scoring=None, cv=None, n_jobs=1,
                    verbose=0, fit_params=None, score_func=None,
                    pre_dispatch='2*n_jobs'):
    """Evaluate a score by cross-validation

    Parameters
    ----------
    estimator : estimator object implementing 'fit'
        The object to use to fit the data.

    X : array-like
        The data to fit. Can be, for example a list, or an array at least 2d.

    y : array-like, optional, default: None
        The target variable to try to predict in the case of
        supervised learning.

    scoring : string, callable or None, optional, default: None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    cv : cross-validation generator, optional, default: None
        A cross-validation generator. If None, a 3-fold cross
        validation is used or 3-fold stratified cross-validation
        when y is supplied and estimator is a classifier.

    n_jobs : integer, optional
        The number of CPUs to use to do the computation. -1 means
        'all CPUs'.

    verbose : integer, optional
        The verbosity level.

    fit_params : dict, optional
        Parameters to pass to the fit method of the estimator.

    pre_dispatch : int, or string, optional
        Controls the number of jobs that get dispatched during parallel
        execution. Reducing this number can be useful to avoid an
        explosion of memory consumption when more jobs get dispatched
        than CPUs can process. This parameter can be:

            - None, in which case all the jobs are immediately
              created and spawned. Use this for lightweight and
              fast-running jobs, to avoid delays due to on-demand
              spawning of the jobs

            - An int, giving the exact number of total jobs that are
              spawned

            - A string, giving an expression as a function of n_jobs,
              as in '2*n_jobs'

    Returns
    -------
    scores : array of float, shape=(len(list(cv)),)
        Array of scores of the estimator for each run of the cross validation.
    """
    X, y = check_arrays(X, y, sparse_format='csr', allow_lists=True,
                        allow_nans=True, allow_nd=True)
    if y is not None:
        y = np.asarray(y)

    cv = _check_cv(cv, X, y, classifier=is_classifier(estimator))
    scorer = check_scoring(estimator, score_func=score_func, scoring=scoring)
    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    parallel = Parallel(n_jobs=n_jobs, verbose=verbose,
                        pre_dispatch=pre_dispatch)
    scores = parallel(delayed(_fit_and_score)(clone(estimator), X, y, scorer,
                                              train, test, verbose, None,
                                              fit_params)
                      for train, test in cv)
    return np.array(scores)[:, 0]


def _fit_and_score(estimator, X, y, scorer, train, test, verbose, parameters,
                   fit_params, return_train_score=False,
                   return_parameters=False):
    """Fit estimator and compute scores for a given dataset split.

    Parameters
    ----------
    estimator : estimator object implementing 'fit'
        The object to use to fit the data.

    X : array-like of shape at least 2D
        The data to fit.

    y : array-like, optional, default: None
        The target variable to try to predict in the case of
        supervised learning.

    scoring : callable
        A scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    train : array-like, shape = (n_train_samples,)
        Indices of training samples.

    test : array-like, shape = (n_test_samples,)
        Indices of test samples.

    verbose : integer
        The verbosity level.

    parameters : dict or None
        Parameters to be set on the estimator.

    fit_params : dict or None
        Parameters that will be passed to ``estimator.fit``.

    return_train_score : boolean, optional, default: False
        Compute and return score on training set.

    return_parameters : boolean, optional, default: False
        Return parameters that has been used for the estimator.

    Returns
    -------
    train_score : float, optional
        Score on training set, returned only if `return_train_score` is `True`.

    test_score : float
        Score on test set.

    n_test_samples : int
        Number of test samples.

    scoring_time : float
        Time spent for fitting and scoring in seconds.

    parameters : dict or None, optional
        The parameters that have been evaluated.
    """
    if verbose > 1:
        if parameters is None:
            msg = "no parameters to be set"
        else:
            msg = '%s' % (', '.join('%s=%s' % (k, v)
                          for k, v in parameters.items()))
        print("[CV] %s %s" % (msg, (64 - len(msg)) * '.'))

    # Adjust lenght of sample weights
    n_samples = _num_samples(X)
    fit_params = fit_params if fit_params is not None else {}
    fit_params = dict([(k, np.asarray(v)[train]
                       if hasattr(v, '__len__') and len(v) == n_samples else v)
                       for k, v in fit_params.items()])

    if parameters is not None:
        estimator.set_params(**parameters)

    start_time = time.time()

    X_train, y_train = _safe_split(estimator, X, y, train)
    X_test, y_test = _safe_split(estimator, X, y, test, train)
    if y_train is None:
        estimator.fit(X_train, **fit_params)
    else:
        estimator.fit(X_train, y_train, **fit_params)
    test_score = _score(estimator, X_test, y_test, scorer)
    if return_train_score:
        train_score = _score(estimator, X_train, y_train, scorer)

    scoring_time = time.time() - start_time

    if verbose > 2:
        msg += ", score=%f" % test_score
    if verbose > 1:
        end_msg = "%s -%s" % (msg, logger.short_format_time(scoring_time))
        print("[CV] %s %s" % ((64 - len(end_msg)) * '.', end_msg))

    ret = [train_score] if return_train_score else []
    ret.extend([test_score, _num_samples(X_test), scoring_time])
    if return_parameters:
        ret.append(parameters)
    return ret


def _safe_split(estimator, X, y, indices, train_indices=None):
    """Create subset of dataset and properly handle kernels."""
    if hasattr(estimator, 'kernel') and callable(estimator.kernel):
        # cannot compute the kernel values with custom function
        raise ValueError("Cannot use a custom kernel function. "
                         "Precompute the kernel matrix instead.")

    if not hasattr(X, "shape"):
        if getattr(estimator, "_pairwise", False):
            raise ValueError("Precomputed kernels or affinity matrices have "
                             "to be passed as arrays or sparse matrices.")
        X_subset = [X[idx] for idx in indices]
    else:
        if getattr(estimator, "_pairwise", False):
            # X is a precomputed square kernel matrix
            if X.shape[0] != X.shape[1]:
                raise ValueError("X should be a square kernel matrix")
            if train_indices is None:
                X_subset = X[np.ix_(indices, indices)]
            else:
                X_subset = X[np.ix_(indices, train_indices)]
        else:
            X_subset = X[safe_mask(X, indices)]

    if y is not None:
        y_subset = y[safe_mask(y, indices)]
    else:
        y_subset = None

    return X_subset, y_subset


def _score(estimator, X_test, y_test, scorer):
    """Compute the score of an estimator on a given test set."""
    if y_test is None:
        score = scorer(estimator, X_test)
    else:
        score = scorer(estimator, X_test, y_test)
    if not isinstance(score, numbers.Number):
        raise ValueError("scoring must return a number, got %s (%s) instead."
                         % (str(score), type(score)))
    return score


def _permutation_test_score(estimator, X, y, cv, scorer):
    """Auxiliary function for permutation_test_score"""
    avg_score = []
    for train, test in cv:
        estimator.fit(X[train], y[train])
        avg_score.append(scorer(estimator, X[test], y[test]))
    return np.mean(avg_score)


def _shuffle(y, labels, random_state):
    """Return a shuffled copy of y eventually shuffle among same labels."""
    if labels is None:
        ind = random_state.permutation(len(y))
    else:
        ind = np.arange(len(labels))
        for label in np.unique(labels):
            this_mask = (labels == label)
            ind[this_mask] = random_state.permutation(ind[this_mask])
    return y[ind]


def permutation_test_score(estimator, X, y, score_func=None, cv=None,
                           n_permutations=100, n_jobs=1, labels=None,
                           random_state=0, verbose=0, scoring=None):
    """Evaluate the significance of a cross-validated score with permutations

    Parameters
    ----------
    estimator : estimator object implementing 'fit'
        The object to use to fit the data.

    X : array-like of shape at least 2D
        The data to fit.

    y : array-like
        The target variable to try to predict in the case of
        supervised learning.

    scoring : string, callable or None, optional, default: None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    cv : integer or cross-validation generator, optional
        If an integer is passed, it is the number of fold (default 3).
        Specific cross-validation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects.

    n_permutations : integer, optional
        Number of times to permute ``y``.

    n_jobs : integer, optional
        The number of CPUs to use to do the computation. -1 means
        'all CPUs'.

    labels : array-like of shape [n_samples] (optional)
        Labels constrain the permutation among groups of samples with
        a same label.

    random_state : RandomState or an int seed (0 by default)
        A random number generator instance to define the state of the
        random permutations generator.

    verbose : integer, optional
        The verbosity level.

    Returns
    -------
    score : float
        The true score without permuting targets.

    permutation_scores : array, shape = [n_permutations]
        The scores obtained for each permutations.

    pvalue : float
        The returned value equals p-value if `score_func` returns bigger
        numbers for better scores (e.g., accuracy_score). If `score_func` is
        rather a loss function (i.e. when lower is better such as with
        `mean_squared_error`) then this is actually the complement of the
        p-value:  1 - p-value.

    Notes
    -----
    This function implements Test 1 in:

        Ojala and Garriga. Permutation Tests for Studying Classifier
        Performance.  The Journal of Machine Learning Research (2010)
        vol. 11

    """
    X, y = check_arrays(X, y, sparse_format='csr', allow_nans=True)
    cv = _check_cv(cv, X, y, classifier=is_classifier(estimator))
    scorer = check_scoring(estimator, scoring=scoring, score_func=score_func)
    random_state = check_random_state(random_state)

    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    score = _permutation_test_score(clone(estimator), X, y, cv, scorer)
    permutation_scores = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(_permutation_test_score)(
            clone(estimator), X, _shuffle(y, labels, random_state), cv,
            scorer)
        for _ in range(n_permutations))
    permutation_scores = np.array(permutation_scores)
    pvalue = (np.sum(permutation_scores >= score) + 1.0) / (n_permutations + 1)
    return score, permutation_scores, pvalue


permutation_test_score.__test__ = False  # to avoid a pb with nosetests


def learning_curve(estimator, X, y, train_sizes=np.linspace(0.1, 1.0, 5),
                   cv=None, scoring=None, exploit_incremental_learning=False,
                   n_jobs=1, pre_dispatch="all", verbose=0):
    """Learning curve.

    Determines cross-validated training and test scores for different training
    set sizes.

    A cross-validation generator splits the whole dataset k times in training
    and test data. Subsets of the training set with varying sizes will be used
    to train the estimator and a score for each training subset size and the
    test set will be computed. Afterwards, the scores will be averaged over
    all k runs for each training subset size.

    Parameters
    ----------
    estimator : object type that implements the "fit" and "predict" methods
        An object of that type which is cloned for each validation.

    X : array-like, shape (n_samples, n_features)
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.

    y : array-like, shape (n_samples) or (n_samples, n_features), optional
        Target relative to X for classification or regression;
        None for unsupervised learning.

    train_sizes : array-like, shape (n_ticks,), dtype float or int
        Relative or absolute numbers of training examples that will be used to
        generate the learning curve. If the dtype is float, it is regarded as a
        fraction of the maximum size of the training set (that is determined
        by the selected validation method), i.e. it has to be within (0, 1].
        Otherwise it is interpreted as absolute sizes of the training sets.
        Note that for classification the number of samples usually have to
        be big enough to contain at least one sample from each class.
        (default: np.linspace(0.1, 1.0, 5))

    cv : integer, cross-validation generator, optional
        If an integer is passed, it is the number of folds (defaults to 3).
        Specific cross-validation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects

    scoring : string, callable or None, optional, default: None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    exploit_incremental_learning : boolean, optional, default: False
        If the estimator supports incremental learning, this will be
        used to speed up fitting for different training set sizes.

    n_jobs : integer, optional
        Number of jobs to run in parallel (default 1).

    pre_dispatch : integer or string, optional
        Number of predispatched jobs for parallel execution (default is
        all). The option can reduce the allocated memory. The string can
        be an expression like '2*n_jobs'.

    verbose : integer, optional
        Controls the verbosity: the higher, the more messages.

    Returns
    -------
    train_sizes_abs : array, shape = (n_unique_ticks,), dtype int
        Numbers of training examples that has been used to generate the
        learning curve. Note that the number of ticks might be less
        than n_ticks because duplicate entries will be removed.

    train_scores : array, shape (n_ticks, n_cv_folds)
        Scores on training sets.

    test_scores : array, shape (n_ticks, n_cv_folds)
        Scores on test set.

    Notes
    -----
    See :ref:`examples/plot_learning_curve.py <example_plot_learning_curve.py>`
    """
    if exploit_incremental_learning and not hasattr(estimator, "partial_fit"):
        raise ValueError("An estimator must support the partial_fit interface "
                         "to exploit incremental learning")

    X, y = check_arrays(X, y, sparse_format='csr', allow_lists=True)
    # Make a list since we will be iterating multiple times over the folds
    cv = list(_check_cv(cv, X, y, classifier=is_classifier(estimator)))
    scorer = check_scoring(estimator, scoring=scoring)

    # HACK as long as boolean indices are allowed in cv generators
    if cv[0][0].dtype == bool:
        new_cv = []
        for i in range(len(cv)):
            new_cv.append((np.nonzero(cv[i][0])[0], np.nonzero(cv[i][1])[0]))
        cv = new_cv

    n_max_training_samples = len(cv[0][0])
    # Because the lengths of folds can be significantly different, it is
    # not guaranteed that we use all of the available training data when we
    # use the first 'n_max_training_samples' samples.
    train_sizes_abs = _translate_train_sizes(train_sizes,
                                             n_max_training_samples)
    n_unique_ticks = train_sizes_abs.shape[0]
    if verbose > 0:
        print("[learning_curve] Training set sizes: " + str(train_sizes_abs))

    parallel = Parallel(n_jobs=n_jobs, pre_dispatch=pre_dispatch,
                        verbose=verbose)
    if exploit_incremental_learning:
        classes = np.unique(y) if is_classifier(estimator) else None
        out = parallel(delayed(_incremental_fit_estimator)(
            clone(estimator), X, y, classes, train, test, train_sizes_abs,
            scorer, verbose) for train, test in cv)
    else:
        out = parallel(delayed(_fit_and_score)(
            clone(estimator), X, y, scorer, train[:n_train_samples], test,
            verbose, parameters=None, fit_params=None, return_train_score=True)
            for train, test in cv for n_train_samples in train_sizes_abs)
        out = np.array(out)[:, :2]
        n_cv_folds = out.shape[0] // n_unique_ticks
        out = out.reshape(n_cv_folds, n_unique_ticks, 2)

    out = np.asarray(out).transpose((2, 1, 0))

    return train_sizes_abs, out[0], out[1]


def _translate_train_sizes(train_sizes, n_max_training_samples):
    """Determine absolute sizes of training subsets and validate 'train_sizes'.

    Examples:
        _translate_train_sizes([0.5, 1.0], 10) -> [5, 10]
        _translate_train_sizes([5, 10], 10) -> [5, 10]

    Parameters
    ----------
    train_sizes : array-like, shape (n_ticks,), dtype float or int
        Numbers of training examples that will be used to generate the
        learning curve. If the dtype is float, it is regarded as a
        fraction of 'n_max_training_samples', i.e. it has to be within (0, 1].

    n_max_training_samples : int
        Maximum number of training samples (upper bound of 'train_sizes').

    Returns
    -------
    train_sizes_abs : array, shape (n_unique_ticks,), dtype int
        Numbers of training examples that will be used to generate the
        learning curve. Note that the number of ticks might be less
        than n_ticks because duplicate entries will be removed.
    """
    train_sizes_abs = np.asarray(train_sizes)
    n_ticks = train_sizes_abs.shape[0]
    n_min_required_samples = np.min(train_sizes_abs)
    n_max_required_samples = np.max(train_sizes_abs)
    if np.issubdtype(train_sizes_abs.dtype, np.float):
        if n_min_required_samples <= 0.0 or n_max_required_samples > 1.0:
            raise ValueError("train_sizes has been interpreted as fractions "
                             "of the maximum number of training samples and "
                             "must be within (0, 1], but is within [%f, %f]."
                             % (n_min_required_samples,
                                n_max_required_samples))
        train_sizes_abs = astype(train_sizes_abs * n_max_training_samples,
                                 dtype=np.int, copy=False)
        train_sizes_abs = np.clip(train_sizes_abs, 1,
                                  n_max_training_samples)
    else:
        if (n_min_required_samples <= 0 or
                n_max_required_samples > n_max_training_samples):
            raise ValueError("train_sizes has been interpreted as absolute "
                             "numbers of training samples and must be within "
                             "(0, %d], but is within [%d, %d]."
                             % (n_max_training_samples,
                                n_min_required_samples,
                                n_max_required_samples))

    train_sizes_abs = np.unique(train_sizes_abs)
    if n_ticks > train_sizes_abs.shape[0]:
        warnings.warn("Removed duplicate entries from 'train_sizes'. Number "
                      "of ticks will be less than than the size of "
                      "'train_sizes' %d instead of %d)."
                      % (train_sizes_abs.shape[0], n_ticks), RuntimeWarning)

    return train_sizes_abs


def _incremental_fit_estimator(estimator, X, y, classes, train, test,
                               train_sizes, scorer, verbose):
    """Train estimator on training subsets incrementally and compute scores."""
    train_scores, test_scores = [], []
    partitions = zip(train_sizes, np.split(train, train_sizes)[:-1])
    for n_train_samples, partial_train in partitions:
        train_subset = train[:n_train_samples]
        X_train, y_train = _safe_split(estimator, X, y, train_subset)
        X_partial_train, y_partial_train = _safe_split(estimator, X, y,
                                                       partial_train)
        X_test, y_test = _safe_split(estimator, X, y, test, train_subset)
        if y_partial_train is None:
            estimator.partial_fit(X_partial_train, classes=classes)
        else:
            estimator.partial_fit(X_partial_train, y_partial_train,
                                  classes=classes)
        train_scores.append(_score(estimator, X_train, y_train, scorer))
        test_scores.append(_score(estimator, X_test, y_test, scorer))
    return np.array((train_scores, test_scores)).T


def validation_curve(estimator, X, y, param_name, param_range, cv=None,
                     scoring=None, n_jobs=1, pre_dispatch="all", verbose=0):
    """Validation curve.

    Determine training and test scores for varying parameter values.

    Compute scores for an estimator with different values of a specified
    parameter. This is similar to grid search with one parameter. However, this
    will also compute training scores and is merely a utility for plotting the
    results.

    Parameters
    ----------
    estimator : object type that implements the "fit" and "predict" methods
        An object of that type which is cloned for each validation.

    X : array-like, shape (n_samples, n_features)
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.

    y : array-like, shape (n_samples) or (n_samples, n_features), optional
        Target relative to X for classification or regression;
        None for unsupervised learning.

    param_name : string
        Name of the parameter that will be varied.

    param_range : array-like, shape (n_values,)
        The values of the parameter that will be evaluated.

    cv : integer, cross-validation generator, optional
        If an integer is passed, it is the number of folds (defaults to 3).
        Specific cross-validation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects

    scoring : string, callable or None, optional, default: None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    n_jobs : integer, optional
        Number of jobs to run in parallel (default 1).

    pre_dispatch : integer or string, optional
        Number of predispatched jobs for parallel execution (default is
        all). The option can reduce the allocated memory. The string can
        be an expression like '2*n_jobs'.

    verbose : integer, optional
        Controls the verbosity: the higher, the more messages.

    Returns
    -------
    train_scores : array, shape (n_ticks, n_cv_folds)
        Scores on training sets.

    test_scores : array, shape (n_ticks, n_cv_folds)
        Scores on test set.

    Notes
    -----
    See
    :ref:`examples/plot_validation_curve.py <example_plot_validation_curve.py>`
    """
    X, y = check_arrays(X, y, sparse_format='csr', allow_lists=True)
    cv = _check_cv(cv, X, y, classifier=is_classifier(estimator))
    scorer = check_scoring(estimator, scoring=scoring)

    parallel = Parallel(n_jobs=n_jobs, pre_dispatch=pre_dispatch,
                        verbose=verbose)
    out = parallel(delayed(_fit_and_score)(
        estimator, X, y, scorer, train, test, verbose,
        parameters={param_name: v}, fit_params=None, return_train_score=True)
        for train, test in cv for v in param_range)

    out = np.asarray(out)[:, :2]
    n_params = len(param_range)
    n_cv_folds = out.shape[0] // n_params
    out = out.reshape(n_cv_folds, n_params, 2).transpose((2, 1, 0))

    return out[0], out[1]

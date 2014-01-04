# Author: Alexander Fabisch <afabisch@informatik.uni-bremen.de>
#
# License: BSD 3 clause

import numpy as np
import warnings
from .base import is_classifier, clone
from .cross_validation import _check_cv
from .utils import check_arrays
from .externals.joblib import Parallel, delayed
from .metrics.scorer import _deprecate_loss_and_score_funcs
from .grid_search import _check_scorable, _split_and_score

def learning_curve(estimator, X, y, n_samples_range=np.linspace(0.1, 1.0, 10),
                   cv=None, scoring=None, exploit_incremental_learning=False,
                   n_jobs=1, pre_dispatch="all", verbose=0):
    """Learning curve

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

    n_samples_range : array-like, shape = [n_ticks,], dtype float or int
        Numbers of training examples that will be used to generate the
        learning curve. If the dtype is float, it is regarded as a
        fraction of the maximum size of the training set (that is determined
        by the selected validation method), i.e. it has to be within (0, 1].
        Note that for classification the number of samples usually have to
        be big enough to contain at least one sample from each class.
        (default: np.linspace(0.1, 1.0, 10))

    cv : integer, cross-validation generator or None, optional, default: None
        If an integer is passed, it is the number of folds (default 3).
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
    n_samples_range : array, shape = [n_unique_ticks,], dtype int
        Numbers of training examples that has been used to generate the
        learning curve. Note that the number of ticks might be less
        than n_ticks because duplicate entries will be removed.

    train_scores : array, shape = [n_ticks,]
        Scores on training sets.

    test_scores : array, shape = [n_ticks,]
        Scores on test set.
    """

    if exploit_incremental_learning and not hasattr(estimator, "partial_fit"):
        raise ValueError("An estimator must support the partial_fit interface "
                         "to exploit incremental learning")

    X, y = check_arrays(X, y, sparse_format='csr', allow_lists=True)
    # Make a list since we will be iterating multiple times over the folds
    cv = list(_check_cv(cv, X, y, classifier=is_classifier(estimator)))

    # HACK as long as boolean indices are allowed in cv generators
    if cv[0][0].dtype == bool:
        new_cv = []
        for i in range(len(cv)):
            new_cv.append((np.nonzero(cv[i][0])[0], np.nonzero(cv[i][1])[0]))
        cv = new_cv

    n_max_training_samples = len(cv[0][0])
    n_samples_range, n_unique_ticks = _translate_n_samples_range(
            n_samples_range, n_max_training_samples)
    # Because the lengths of folds can be significantly different, it is
    # not guaranteed that we use all of the available training data when we
    # use the first 'n_max_training_samples' samples.
    if verbose > 0:
        print("[learning_curve] Training set sizes: " + str(n_samples_range))

    _check_scorable(estimator, scoring=scoring)
    scorer = _deprecate_loss_and_score_funcs(scoring=scoring)

    parallel = Parallel(n_jobs=n_jobs, pre_dispatch=pre_dispatch,
                        verbose=verbose)
    if exploit_incremental_learning:
        if is_classifier(estimator):
            classes = np.unique(y)
        else:
            classes = None
        out = parallel(delayed(_incremental_fit_estimator)(
                           estimator, X, y, classes, train, test,
                           n_samples_range, scorer, verbose)
                       for train, test in cv)
    else:
        out = parallel(delayed(_fit_estimator)(
                           estimator, X, y, train, test, n_train_samples,
                           scorer, verbose)
                       for train, test in cv
                       for n_train_samples in n_samples_range)
        out = np.array(out)
        n_cv_folds = out.shape[0]/n_unique_ticks
        out = out.reshape(n_cv_folds, n_unique_ticks, 2)

    avg_over_cv = np.asarray(out).mean(axis=0).reshape(n_unique_ticks, 2)

    return n_samples_range, avg_over_cv[:, 0], avg_over_cv[:, 1]


def _translate_n_samples_range(n_samples_range, n_max_training_samples):
    """Determine range of number of training samples"""
    n_samples_range = np.asarray(n_samples_range)
    n_ticks = n_samples_range.shape[0]
    n_min_required_samples = np.min(n_samples_range)
    n_max_required_samples = np.max(n_samples_range)
    if np.issubdtype(n_samples_range.dtype, np.float):
        if n_min_required_samples <= 0.0 or n_max_required_samples > 1.0:
            raise ValueError("n_samples_range must be within (0, 1], "
                             "but is within [%f, %f]."
                             % (n_min_required_samples,
                                n_max_required_samples))
        n_samples_range = (n_samples_range * n_max_training_samples
                ).astype(np.int)
        n_samples_range = np.clip(n_samples_range, 1, n_max_training_samples)
    else:
        if (n_min_required_samples <= 0 or
            n_max_required_samples > n_max_training_samples):
            raise ValueError("n_samples_range must be within (0, %d], "
                             "but is within [%d, %d]."
                             % (n_max_training_samples,
                                n_min_required_samples,
                                n_max_required_samples))

    n_samples_range = np.unique(n_samples_range)
    n_unique_ticks = n_samples_range.shape[0]
    if n_ticks > n_unique_ticks:
        warnings.warn("Number of ticks will be less than than the size of "
                      "'n_samples_range' (%d instead of %d)."
                      % (n_unique_ticks, n_ticks), RuntimeWarning)

    return n_samples_range, n_unique_ticks


def _fit_estimator(base_estimator, X, y, train, test,
                   n_train_samples, scorer, verbose):
    estimator = clone(base_estimator)
    test_score, _, train_score, _ = _split_and_score(
            estimator, X, y, train=train[:n_train_samples],
            test=test, scorer=scorer, return_train_score=True)
    return train_score, test_score


def _incremental_fit_estimator(base_estimator, X, y, classes, train, test,
                               n_samples_range, scorer, verbose):
    estimator = clone(base_estimator)
    train_scores, test_scores = [], []
    for n_train_samples, partial_train in zip(n_samples_range, np.split(train,
            n_samples_range)[:-1]):
        test_score, _, train_score, _ = _split_and_score(
                estimator, X, y, train=train[:n_train_samples],
                partial_train=partial_train, test=test, scorer=scorer,
                return_train_score=True, classes=classes)
        train_scores.append(train_score)
        test_scores.append(test_score)
    return np.array((train_scores, test_scores)).T

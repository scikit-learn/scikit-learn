import numpy as np
from .base import is_classifier, clone
from .cross_validation import _check_cv
from .utils import check_arrays
from .externals.joblib import Parallel, delayed
from .metrics.scorer import _deprecate_loss_and_score_funcs
from .grid_search import _check_scorable, _split_and_score

def learning_curve(estimator, X, y,
                   n_samples_range=np.linspace(0.1, 1.0, 10), cv=None, scoring=None,
                   n_jobs=1, verbose=0):
    """ TODO document me

    Parameters
    ----------
    estimator : object type that implements the "fit" and "predict" methods
        An object of that type is instantiated for each validation.

    X : array-like, shape = [n_samples, n_features]
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.

    y : array-like, shape = [n_samples] or [n_samples, n_output], optional
        Target relative to X for classification or regression;
        None for unsupervised learning.

    n_samples_range : array-like, shape = [n_ticks,], dtype float or int
        Numbers of training examples that will be used to generate the
        learning curve. If the dtype is float, it is regarded as a
        fraction of n_samples, i.e. it has to be within ]0, 1].
        (default: np.linspace(0.1, 1.0, 10))

    cv : integer, cross-validation generator or None, optional, default: None
        If an integer is passed, it is the number of folds (default 3).
        Specific cross-validation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects

    scoring : string, callable or None, optional, default: None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    n_jobs : integer, optional
        Number of jobs to run in parallel (default 1).

    verbose : integer, optional
        Controls the verbosity: the higher, the more messages.

    Returns
    -------
    n_samples_range : array, shape = [n_ticks,], dtype int
        Numbers of training examples that has been used to generate the
        learning curve.

    train_scores : array, shape = [n_ticks,]
        Scores on training sets.

    test_scores : array, shape = [n_ticks,]
        Scores on test set.
    """
    # TODO tests, doc
    # TODO allow y to be None for unsupervised learning
    # TODO there is an overlap with grid search -> refactoring
    # TODO exploit incremental learning
    # TODO use verbose argument

    X, y = check_arrays(X, y, sparse_format='csr', allow_lists=True)
    # Make a list since we will be iterating multiple times over the folds
    cv = list(_check_cv(cv, X, y, classifier=is_classifier(estimator)))

    # Determine range of number of training samples
    n_max_training_samples = cv[0][0].shape[0]
    n_samples_range = np.asarray(n_samples_range)
    n_min_required_samples = np.min(n_samples_range)
    n_max_required_samples = np.max(n_samples_range)
    if np.issubdtype(n_samples_range.dtype, np.float):
        if n_min_required_samples <= 0.0 or n_max_required_samples > 1.0:
            raise ValueError("n_samples_range must be within ]0, 1], "
                             "but is within [%f, %f]."
                             % (n_min_required_samples,
                                n_max_required_samples))
        n_samples_range = (n_samples_range *
                           n_max_training_samples).astype(np.int)
    else:
        if (n_min_required_samples <= 0 or
            n_max_required_samples > n_max_training_samples):
            raise ValueError("n_samples_range must be within ]0, %d], "
                             "but is within [%d, %d]."
                             % (n_max_training_samples,
                                n_min_required_samples,
                                n_max_required_samples))

    _check_scorable(estimator, scoring=scoring)
    scorer = _deprecate_loss_and_score_funcs(scoring=scoring)

    out = Parallel(
        # TODO use pre_dispatch parameter? what is it good for?
        n_jobs=n_jobs, verbose=verbose)(
            delayed(_fit_estimator)(
                estimator, X, y, train, test, n_train_samples,
                scorer, verbose)
            for train, test in cv for n_train_samples in n_samples_range)

    out = np.asarray(out)
    train_scores = np.zeros(n_samples_range.shape, dtype=np.float)
    test_scores = np.zeros(n_samples_range.shape, dtype=np.float)
    for i, n_train_samples in enumerate(n_samples_range):
        res_indices = np.where(out[:, 0] == n_train_samples)
        train_scores[i], test_scores[i] = out[res_indices[0], 1:].mean(axis=0)

    return n_samples_range, train_scores, test_scores

def _fit_estimator(base_estimator, X, y, train, test, n_train_samples,
                   scorer, verbose):
    test_score, _, train_score, _ = _split_and_score(
            base_estimator, X, y, parameters={}, train=train[:n_train_samples],
            test=test, scorer=scorer, return_train_score=True)
    return n_train_samples, train_score, test_score

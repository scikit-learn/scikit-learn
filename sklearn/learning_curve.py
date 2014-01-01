import numpy as np
from .base import is_classifier, clone
from .cross_validation import _check_cv
from .utils import check_arrays
from .externals.joblib import Parallel, delayed
from .metrics.scorer import _deprecate_loss_and_score_funcs

def learning_curve(estimator, X, y,
                   n_samples_range=np.linspace(0.1, 1.0, 10), cv=None, scoring=None,
                   n_jobs=1, verbose=False, random_state=None):
    """ TODO document me
    Parameters
    ----------
    n_samples_range : array-like with dtype float or int,
        If the dtype is float, it is regarded as a fraction of n_samples, i.e. it has to be within ]0, 1].
    """
    # TODO tests, doc
    # TODO allow y to be None for unsupervised learning
    # TODO there is an overlap with grid search -> refactoring
    # TODO exploit incremental learning

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

    # TODO copied from BaseGridSearch -> move to utils? .base? where?
    if (not hasattr(estimator, 'fit') or
            not (hasattr(estimator, 'predict')
                  or hasattr(estimator, 'score'))):
        raise TypeError("estimator should a be an estimator implementing"
                        " 'fit' and 'predict' or 'score' methods,"
                        " %s (type %s) was passed" %
                        (estimator, type(estimator)))
    if scoring is None:
        if not hasattr(estimator, 'score'):
            raise TypeError(
                "If no scoring is specified, the estimator passed "
                "should have a 'score' method. The estimator %s "
                "does not." % estimator)
        scorer = _deprecate_loss_and_score_funcs(scoring=scoring)

    out = Parallel(
        # TODO use pre_dispatch parameter? what is it good for?
        n_jobs=n_jobs, verbose=verbose)(
            delayed(_fit_estimator)(
                estimator, X, y, train[:n_train_samples], test, scorer,
                verbose)
            for train, test in cv for n_train_samples in n_samples_range)

    out = np.asarray(out)
    train_scores = np.zeros(n_samples_range.shape, dtype=np.float)
    test_scores = np.zeros(n_samples_range.shape, dtype=np.float)
    for i, n_train_samples in enumerate(n_samples_range):
        res_indices = np.where(out[:, 0] == n_train_samples)
        train_scores[i], test_scores[i] = out[res_indices[0], 1:].mean(axis=0)

    return n_samples_range, train_scores, test_scores

def _fit_estimator(base_estimator, X, y, train, test, scorer, verbose):
    # TODO similar to fit_grid_point from grid search, refactor
    estimator = clone(base_estimator)
    estimator.fit(X[train], y[train])
    if scorer is None:
        train_score = estimator.score(X[train], y[train])
        test_score = estimator.score(X[test], y[test])
    else:
        train_score = scorer(estimator, X[train], y[train])
        test_score = scorer(estimator, X[test], y[test])
    return train.shape[0], train_score, test_score

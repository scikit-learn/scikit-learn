import numpy as np
from .base import clone
from .cross_validation import KFold
from .externals.joblib import Parallel, delayed
from .metrics.scorer import _deprecate_loss_and_score_funcs

def learning_curve(estimator, X, y, n_samples_range=None, step_size=1,
                   n_cv_folds=10, loss_func=None, scoring=None,
                   n_jobs=1, verbose=False, random_state=None):
    # TODO tests, doc
    # TODO allow y to be None for unsupervised learning
    # TODO test different n_cv_folds / dataset sizes / etc. (there could be bugs)
    # TODO there is a huge overlap with grid search -> refactoring
    # TODO exploit incremental learning?! (might be a bit complicated with CV)

    n_samples = X.shape[0]
    max_fold_size = n_samples / n_cv_folds

    if n_samples_range is None:
        if step_size is None or step_size < 1:
            raise ValueError("Define either a range of training set sizes or "
                             "a proper step size.")
        n_samples_range = np.arange(n_cv_folds-1, n_samples-max_fold_size+1,
                                    step_size)

    n_max_samples = np.max(n_samples_range)
    n_required_samples = n_max_samples + max_fold_size
    if n_samples < n_required_samples:
        raise ValueError(
                "For %d-fold cross-validation with %d training examples, "
                "%d samples are required (got %d)."
                % (n_cv_folds, n_max_samples, n_required_samples, n_samples))

    # TODO copied from BaseGridSearch -> move to utils? .base? where?
    if (not hasattr(estimator, 'fit') or
            not (hasattr(estimator, 'predict')
                  or hasattr(estimator, 'score'))):
        raise TypeError("estimator should a be an estimator implementing"
                        " 'fit' and 'predict' or 'score' methods,"
                        " %s (type %s) was passed" %
                        (estimator, type(estimator)))
    if scoring is None and loss_func is None:
        if not hasattr(estimator, 'score'):
            raise TypeError(
                "If no scoring is specified, the estimator passed "
                "should have a 'score' method. The estimator %s "
                "does not." % estimator)
        scorer = _deprecate_loss_and_score_funcs(loss_func=loss_func,
                                                 scoring=scoring)

    scores = []
    for n_train_samples in n_samples_range:
        # TODO maybe take random indices instead of the first slice_length?
        fold_size = (n_train_samples+1) / n_cv_folds
        slice_length = n_train_samples + fold_size
        cv = KFold(n=slice_length, n_folds=n_cv_folds,
                   random_state=random_state)

        out = Parallel(
            # TODO set pre_dispatch parameter? what is it good for?
            n_jobs=n_jobs, verbose=verbose)(
                delayed(_fit_estimator)(estimator, X, y, train, test, scorer,
                                       verbose)
                for train, test in cv)
        scores.append(np.mean(out, axis=0))
    scores = np.array(scores)

    return n_samples_range, scores[:, 0], scores[:, 1]

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
    return train_score, test_score

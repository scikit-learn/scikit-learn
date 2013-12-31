import numpy as np
from sklearn import clone # TODO relative import
from sklearn.cross_validation import KFold # TODO relative import
from sklearn.externals.joblib import Parallel, delayed # TODO relative import
#from .metrics import _deprecate_loss_and_score_funcs # TODO relative import

def learning_curve(estimator, X, y, n_samples_range,
                   loss_func=None, scoring='accuracy', n_cv_folds=10, n_jobs=1,
                   verbose=False, random_state=None):
    n_max_samples = np.max(n_samples_range)
    n_samples = X.shape[0]
    n_required_samples = int(n_max_samples * (1 + 1.0 / n_cv_folds))
    if n_samples < n_required_samples:
        # TODO test case
        raise ValueError(
                "For %d-fold cross-validation with %d training examples, "
                "%d samples are required (got %d)."
                % (n_cv_folds, n_max_samples, n_required_samples, n_samples))

    # TODO check if scoring is possible (see BaseGridSearch)

    #scorer = _deprecate_loss_and_score_funcs(loss_func=loss_func, scoring=scoring) # TODO enable with relative import
    from sklearn.metrics import accuracy_score
    scorer = accuracy_score

    scores = []
    for n_train_samples in n_samples_range:
        # TODO maybe take random indices instead of the first slice_length?
        slice_length = int(n_train_samples * (1 + 1.0 / n_cv_folds))
        cv = KFold(n=slice_length, n_folds=n_cv_folds,
                   random_state=random_state)

        out = Parallel(
            n_jobs=n_jobs, verbose=verbose)( # TODO set pre_dispatch parameter?
                delayed(fit_estimator)(estimator, X, y, train, test, scorer, verbose)
                for train, test in cv)
        scores.append(np.mean(out, axis=0))
    scores = np.array(scores)
    print scores
    return scores[:, 0], scores[:, 1]

def fit_estimator(base_estimator, X, y, train, test, scorer, verbose):
    estimator = clone(base_estimator)
    estimator.fit(X[train], y[train])
    y_test_pred = estimator.predict(X[test])
    y_train_pred = estimator.predict(X[train])
    train_score = scorer(y[train], y_train_pred)
    test_score = scorer(y[test], y_test_pred)
    return train_score, test_score

if __name__ == "__main__":
    #from sklearn.linear_model import RidgeClassifier
    #estimator = RidgeClassifier(alpha=100)
    from sklearn.svm import SVC
    estimator = SVC(gamma=0.001)

    from sklearn.datasets import load_digits
    digits = load_digits()
    X, y = digits.data, digits.target

    n_samples_range = np.arange(10, 1611, 100)
    train_scores, test_scores = learning_curve(estimator, X, y, n_samples_range, n_jobs=2, verbose=False)
    import pylab
    pylab.plot(n_samples_range, train_scores)
    pylab.plot(n_samples_range, test_scores)
    pylab.show()

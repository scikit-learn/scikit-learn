import numpy as np
import pytest
from sklearn.datasets import make_classification, make_regression
from sklearn.utils.estimator_checks import check_estimator

from sklearn.experimental import HistGradientBoostingClassifier
from sklearn.experimental import HistGradientBoostingRegressor


X_classification, y_classification = make_classification(random_state=0)
X_regression, y_regression = make_regression(random_state=0)


@pytest.mark.parametrize('GradientBoosting, X, y', [
    (HistGradientBoostingClassifier, X_classification, y_classification),
    (HistGradientBoostingRegressor, X_regression, y_regression)
])
def test_init_parameters_validation(GradientBoosting, X, y):

    with pytest.raises(
            ValueError,
            match="Loss blah is not supported for"):
        GradientBoosting(loss='blah').fit(X, y)

    for learning_rate in (-1, 0):
        with pytest.raises(
                ValueError,
                match="learning_rate={} must be strictly positive".format(
                    learning_rate)):
            GradientBoosting(learning_rate=learning_rate).fit(X, y)

    with pytest.raises(
            ValueError,
            match="max_iter=0 must not be smaller than 1"):
        GradientBoosting(max_iter=0).fit(X, y)

    with pytest.raises(
            ValueError,
            match="max_leaf_nodes=0 should not be smaller than 1"):
        GradientBoosting(max_leaf_nodes=0).fit(X, y)

    with pytest.raises(
            ValueError,
            match="max_depth=0 should not be smaller than 1"):
        GradientBoosting(max_depth=0).fit(X, y)

    with pytest.raises(
            ValueError,
            match="min_samples_leaf=0 should not be smaller than 1"):
        GradientBoosting(min_samples_leaf=0).fit(X, y)

    with pytest.raises(
            ValueError,
            match="l2_regularization=-1 must be positive"):
        GradientBoosting(l2_regularization=-1).fit(X, y)

    for max_bins in (1, 257):
        with pytest.raises(
                ValueError,
                match="max_bins={} should be no smaller than 2 and "
                      "no larger".format(max_bins)):
            GradientBoosting(max_bins=max_bins).fit(X, y)

    with pytest.raises(
            ValueError,
            match="n_iter_no_change=-1 must be positive"):
        GradientBoosting(n_iter_no_change=-1).fit(X, y)

    for validation_fraction in (-1, 0):
        with pytest.raises(
            ValueError,
            match="validation_fraction={} must be strictly positive".format(
                validation_fraction)):
            GradientBoosting(validation_fraction=validation_fraction).fit(X, y)

    with pytest.raises(
            ValueError,
            match="tol=-1 must not be smaller than 0"):
        GradientBoosting(tol=-1).fit(X, y)


def test_invalid_classification_loss():
    binary_clf = HistGradientBoostingClassifier(loss="binary_crossentropy")
    with pytest.raises(
            ValueError,
            match="loss='binary_crossentropy' is not defined for multiclass"
                  " classification with n_classes=3, use"
                  " loss='categorical_crossentropy' instead"):
        binary_clf.fit(np.zeros(shape=(3, 2)), np.arange(3))


@pytest.mark.parametrize(
    'scoring, validation_fraction, n_iter_no_change, tol', [
        ('neg_mean_squared_error', .1, 5, 1e-7),  # use scorer
        ('neg_mean_squared_error', None, 5, 1e-1),  # use scorer on train data
        (None, .1, 5, 1e-7),  # same with default scorer
        (None, None, 5, 1e-1),
        ('loss', .1, 5, 1e-7),  # use loss
        ('loss', None, 5, 1e-1),  # use loss on training data
        (None, None, None, None),  # no early stopping
        ])
def test_early_stopping_regression(scoring, validation_fraction,
                                   n_iter_no_change, tol):

    max_iter = 200

    X, y = make_regression(random_state=0)

    gb = HistGradientBoostingRegressor(verbose=1,  # just for coverage
                                       scoring=scoring,
                                       tol=tol,
                                       validation_fraction=validation_fraction,
                                       max_iter=max_iter,
                                       n_iter_no_change=n_iter_no_change,
                                       random_state=0)
    gb.fit(X, y)

    if n_iter_no_change is not None:
        assert n_iter_no_change <= gb.n_iter_ < max_iter
    else:
        assert gb.n_iter_ == max_iter


@pytest.mark.parametrize('data', (
    make_classification(random_state=0),
    make_classification(n_classes=3, n_clusters_per_class=1, random_state=0)
))
@pytest.mark.parametrize(
    'scoring, validation_fraction, n_iter_no_change, tol', [
        ('accuracy', .1, 5, 1e-7),  # use scorer
        ('accuracy', None, 5, 1e-1),  # use scorer on training data
        (None, .1, 5, 1e-7),  # same with default scorerscor
        (None, None, 5, 1e-1),
        ('loss', .1, 5, 1e-7),  # use loss
        ('loss', None, 5, 1e-1),  # use loss on training data
        (None, None, None, None),  # no early stopping
        ])
def test_early_stopping_classification(data, scoring, validation_fraction,
                                       n_iter_no_change, tol):

    max_iter = 50

    X, y = data

    gb = HistGradientBoostingClassifier(
        verbose=1,  # just for coverage
        scoring=scoring,
        tol=tol,
        validation_fraction=validation_fraction,
        max_iter=max_iter,
        n_iter_no_change=n_iter_no_change,
        random_state=0)
    gb.fit(X, y)

    if n_iter_no_change is not None:
        assert n_iter_no_change <= gb.n_iter_ < max_iter
    else:
        assert gb.n_iter_ == max_iter


def test_should_stop():

    def should_stop(scores, n_iter_no_change, tol):
        gbdt = HistGradientBoostingClassifier(
            n_iter_no_change=n_iter_no_change,
            tol=tol)
        return gbdt._should_stop(scores)

    # not enough iterations
    assert not should_stop([], n_iter_no_change=1, tol=0.001)

    assert not should_stop([1, 1, 1], n_iter_no_change=5, tol=0.001)
    assert not should_stop([1] * 5, n_iter_no_change=5, tol=0.001)

    # still making significant progress up to tol
    assert not should_stop([1, 2, 3, 4, 5, 6], n_iter_no_change=5, tol=0.001)
    assert not should_stop([1, 2, 3, 4, 5, 6], n_iter_no_change=5, tol=0.)
    assert not should_stop([1, 2, 3, 4, 5, 6], n_iter_no_change=5, tol=0.999)
    assert not should_stop([1, 2, 3, 4, 5, 6], n_iter_no_change=5,
                           tol=5 - 1e-5)

    # no significant progress according to tol
    assert should_stop([1] * 6, n_iter_no_change=5, tol=0.)
    assert should_stop([1] * 6, n_iter_no_change=5, tol=0.001)
    assert should_stop([1, 2, 3, 4, 5, 6], n_iter_no_change=5, tol=5)


@pytest.mark.parametrize('Estimator', (
    HistGradientBoostingRegressor(),
    HistGradientBoostingClassifier(),
    ))
def test_estimator_checks(Estimator):
    # Run the check_estimator() test suite on GBRegressor and GBClassifier.
    # Just here for convenience, must be removed before merging since these
    # tests are run in test_common anyways
    check_estimator(Estimator)

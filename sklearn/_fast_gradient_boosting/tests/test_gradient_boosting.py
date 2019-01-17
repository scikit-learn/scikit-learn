import pytest
from sklearn.utils.testing import assert_raises_regex
from sklearn.datasets import make_classification, make_regression
from sklearn.utils.estimator_checks import check_estimator

from sklearn.ensemble import FastGradientBoostingClassifier
from sklearn.ensemble import FastGradientBoostingRegressor


X_classification, y_classification = make_classification(random_state=0)
X_regression, y_regression = make_regression(random_state=0)


@pytest.mark.parametrize('GradientBoosting, X, y', [
    # (GBMClassifier, X_classification, y_classification),  TODO: unskip
    (FastGradientBoostingRegressor, X_regression, y_regression)
])
def test_init_parameters_validation(GradientBoosting, X, y):

    assert_raises_regex(
        ValueError,
        "Loss blah is not supported for",
        GradientBoosting(loss='blah').fit, X, y
    )

    for learning_rate in (-1, 0):
        assert_raises_regex(
            ValueError,
            f"learning_rate={learning_rate} must be strictly positive",
            GradientBoosting(learning_rate=learning_rate).fit, X, y
        )

    assert_raises_regex(
        ValueError,
        f"max_iter=0 must not be smaller than 1",
        GradientBoosting(max_iter=0).fit, X, y
    )

    assert_raises_regex(
        ValueError,
        f"max_leaf_nodes=0 should not be smaller than 1",
        GradientBoosting(max_leaf_nodes=0).fit, X, y
    )

    assert_raises_regex(
        ValueError,
        f"max_depth=0 should not be smaller than 1",
        GradientBoosting(max_depth=0).fit, X, y
    )

    assert_raises_regex(
        ValueError,
        f"min_samples_leaf=0 should not be smaller than 1",
        GradientBoosting(min_samples_leaf=0).fit, X, y
    )

    assert_raises_regex(
        ValueError,
        f"l2_regularization=-1 must be positive",
        GradientBoosting(l2_regularization=-1).fit, X, y
    )

    for max_bins in (1, 257):
        assert_raises_regex(
            ValueError,
            f"max_bins={max_bins} should be no smaller than 2 and no larger",
            GradientBoosting(max_bins=max_bins).fit, X, y
        )

    assert_raises_regex(
        ValueError,
        f"n_iter_no_change=-1 must be positive",
        GradientBoosting(n_iter_no_change=-1).fit, X, y
    )

    for validation_split in (-1, 0):
        assert_raises_regex(
            ValueError,
            f"validation_split={validation_split} must be strictly positive",
            GradientBoosting(validation_split=validation_split).fit, X, y
        )

    assert_raises_regex(
        ValueError,
        f"tol=-1 must not be smaller than 0",
        GradientBoosting(tol=-1).fit, X, y
    )


@pytest.mark.parametrize('scoring, validation_split, n_iter_no_change, tol', [
    ('neg_mean_squared_error', .1, 5, 1e-7),  # use scorer
    ('neg_mean_squared_error', None, 5, 1e-1),  # use scorer on training data
    (None, .1, 5, 1e-7),  # use loss
    (None, None, 5, 1e-1),  # use loss on training data
    (None, None, None, None),  # no early stopping
])
def test_early_stopping_regression(scoring, validation_split,
                                   n_iter_no_change, tol):

    max_iter = 500

    X, y = make_regression(random_state=0)

    gb = FastGradientBoostingRegressor(verbose=1,  # just for coverage
                      scoring=scoring,
                      tol=tol,
                      validation_split=validation_split,
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
@pytest.mark.parametrize('scoring, validation_split, n_iter_no_change, tol', [
    ('accuracy', .1, 5, 1e-7),  # use scorer
    ('accuracy', None, 5, 1e-1),  # use scorer on training data
    (None, .1, 5, 1e-7),  # use loss
    (None, None, 5, 1e-1),  # use loss on training data
    (None, None, None, None),  # no early stopping
])
def test_early_stopping_classification(data, scoring, validation_split,
                                       n_iter_no_change, tol):

    max_iter = 500

    X, y = data

    gb = FastGradientBoostingClassifier(verbose=1,  # just for coverage
                       scoring=scoring,
                       tol=tol,
                       validation_split=validation_split,
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
        gbdt = FastGradientBoostingClassifier(n_iter_no_change=n_iter_no_change,
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
    FastGradientBoostingRegressor(),
    FastGradientBoostingClassifier(scoring=None, validation_split=None, min_samples_leaf=5),
    ))
def test_estimator_checks(Estimator):
    # Run the check_estimator() test suite on GBRegressor and GBClassifier.

    # Notes:
    # - Can't do early stopping with classifier because often
    #   validation_split=.1 leads to test_size=2 < n_classes and
    #   train_test_split raises an error.
    # - Also, need to set a low min_samples_leaf for
    #   check_classifiers_classes() to pass: with only 30 samples on the
    #   dataset, the root is never split with min_samples_leaf=20 and only the
    #   majority class is predicted.
    check_estimator(Estimator)

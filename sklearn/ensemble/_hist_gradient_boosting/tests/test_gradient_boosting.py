import numpy as np
import pytest
from sklearn.base import clone
from sklearn.datasets import make_classification, make_regression
from sklearn.utils.testing import assert_equal, assert_not_equal

# To use this experimental feature, we need to explicitly ask for it:
from sklearn.experimental import enable_hist_gradient_boosting  # noqa
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.ensemble._hist_gradient_boosting.binning import _BinMapper


X_classification, y_classification = make_classification(random_state=0)
X_regression, y_regression = make_regression(random_state=0)

X_classification_large, y_classification_large = make_classification(
    n_samples=20000, random_state=0)
X_regression_large, y_regression_large = make_regression(
    n_samples=20000, random_state=0)


def _assert_predictor_equal(gb_1, gb_2, X):
    """Assert that two HistGBM instances are identical."""
    # Check identical nodes for each tree
    for (pred_ith_1, pred_ith_2) in zip(gb_1._predictors, gb_2._predictors):
        for (predictor_1, predictor_2) in zip(pred_ith_1, pred_ith_2):
            np.testing.assert_array_equal(
                predictor_1.nodes,
                predictor_2.nodes
            )

    # Check identical predictions
    np.testing.assert_allclose(gb_1.predict(X), gb_2.predict(X))


@pytest.mark.parametrize('GradientBoosting, X, y', [
    (HistGradientBoostingClassifier, X_classification, y_classification),
    (HistGradientBoostingRegressor, X_regression, y_regression)
])
@pytest.mark.parametrize(
    'params, err_msg',
    [({'loss': 'blah'}, 'Loss blah is not supported for'),
     ({'learning_rate': 0}, 'learning_rate=0 must be strictly positive'),
     ({'learning_rate': -1}, 'learning_rate=-1 must be strictly positive'),
     ({'max_iter': 0}, 'max_iter=0 must not be smaller than 1'),
     ({'max_leaf_nodes': 0}, 'max_leaf_nodes=0 should not be smaller than 2'),
     ({'max_leaf_nodes': 1}, 'max_leaf_nodes=1 should not be smaller than 2'),
     ({'max_depth': 0}, 'max_depth=0 should not be smaller than 2'),
     ({'max_depth': 1}, 'max_depth=1 should not be smaller than 2'),
     ({'min_samples_leaf': 0}, 'min_samples_leaf=0 should not be smaller'),
     ({'l2_regularization': -1}, 'l2_regularization=-1 must be positive'),
     ({'max_bins': 1}, 'max_bins=1 should be no smaller than 2 and no larger'),
     ({'max_bins': 257}, 'max_bins=257 should be no smaller than 2 and no'),
     ({'n_iter_no_change': -1}, 'n_iter_no_change=-1 must be positive'),
     ({'validation_fraction': -1}, 'validation_fraction=-1 must be strictly'),
     ({'validation_fraction': 0}, 'validation_fraction=0 must be strictly'),
     ({'tol': -1}, 'tol=-1 must not be smaller than 0')]
)
def test_init_parameters_validation(GradientBoosting, X, y, params, err_msg):

    with pytest.raises(ValueError, match=err_msg):
        GradientBoosting(**params).fit(X, y)


def test_invalid_classification_loss():
    binary_clf = HistGradientBoostingClassifier(loss="binary_crossentropy")
    err_msg = ("loss='binary_crossentropy' is not defined for multiclass "
               "classification with n_classes=3, use "
               "loss='categorical_crossentropy' instead")
    with pytest.raises(ValueError, match=err_msg):
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

    gb = HistGradientBoostingRegressor(
        verbose=1,  # just for coverage
        min_samples_leaf=5,  # easier to overfit fast
        scoring=scoring,
        tol=tol,
        validation_fraction=validation_fraction,
        max_iter=max_iter,
        n_iter_no_change=n_iter_no_change,
        random_state=0
    )
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
        min_samples_leaf=5,  # easier to overfit fast
        scoring=scoring,
        tol=tol,
        validation_fraction=validation_fraction,
        max_iter=max_iter,
        n_iter_no_change=n_iter_no_change,
        random_state=0
    )
    gb.fit(X, y)

    if n_iter_no_change is not None:
        assert n_iter_no_change <= gb.n_iter_ < max_iter
    else:
        assert gb.n_iter_ == max_iter


@pytest.mark.parametrize(
    'scores, n_iter_no_change, tol, stopping',
    [
        ([], 1, 0.001, False),  # not enough iterations
        ([1, 1, 1], 5, 0.001, False),  # not enough iterations
        ([1, 1, 1, 1, 1], 5, 0.001, False),  # not enough iterations
        ([1, 2, 3, 4, 5, 6], 5, 0.001, False),  # significant improvement
        ([1, 2, 3, 4, 5, 6], 5, 0., False),  # significant improvement
        ([1, 2, 3, 4, 5, 6], 5, 0.999, False),  # significant improvement
        ([1, 2, 3, 4, 5, 6], 5, 5 - 1e-5, False),  # significant improvement
        ([1] * 6, 5, 0., True),  # no significant improvement
        ([1] * 6, 5, 0.001, True),  # no significant improvement
        ([1] * 6, 5, 5, True),  # no significant improvement
    ]
)
def test_should_stop(scores, n_iter_no_change, tol, stopping):

    gbdt = HistGradientBoostingClassifier(
        n_iter_no_change=n_iter_no_change, tol=tol
    )
    assert gbdt._should_stop(scores) == stopping


def test_binning_train_validation_are_separated():
    # Make sure training and validation data are binned separately.
    # See issue 13926

    rng = np.random.RandomState(0)
    validation_fraction = .2
    gb = HistGradientBoostingClassifier(
        n_iter_no_change=5,
        validation_fraction=validation_fraction,
        random_state=rng
    )
    gb.fit(X_classification, y_classification)
    mapper_training_data = gb.bin_mapper_

    # Note that since the data is small there is no subsampling and the
    # random_state doesn't matter
    mapper_whole_data = _BinMapper(random_state=0)
    mapper_whole_data.fit(X_classification)

    n_samples = X_classification.shape[0]
    assert np.all(mapper_training_data.actual_n_bins_ ==
                  int((1 - validation_fraction) * n_samples))
    assert np.all(mapper_training_data.actual_n_bins_ !=
                  mapper_whole_data.actual_n_bins_)


@pytest.mark.parametrize('GradientBoosting, X, y', [
    (HistGradientBoostingClassifier, X_classification, y_classification),
    (HistGradientBoostingRegressor, X_regression, y_regression)
])
def test_max_iter_with_warm_start_validation(GradientBoosting, X, y):
    # Check that a ValueError is raised when the maximum number of iterations
    # is smaller than the number of iterations from the previous fit when warm
    # start is True.

    estimator = GradientBoosting(max_iter=50, warm_start=True)
    estimator.fit(X, y)
    estimator.set_params(max_iter=25)
    err_msg = ('max_iter=25 must be larger than or equal to n_iter_=50 '
               'when warm_start==True')
    with pytest.raises(ValueError, match=err_msg):
        estimator.fit(X, y)


@pytest.mark.parametrize('GradientBoosting, X, y', [
    (HistGradientBoostingClassifier, X_classification, y_classification),
    (HistGradientBoostingRegressor, X_regression, y_regression)
])
def test_warm_start_yields_identical_results(GradientBoosting, X, y):
    # Make sure that fitting 50 iterations and then 25 with warm start is
    # equivalent to fitting 75 iterations.

    rng = 42
    gb_warm_start = GradientBoosting(
        n_iter_no_change=100, max_iter=50, random_state=rng, warm_start=True
    )
    gb_warm_start.fit(X, y).set_params(max_iter=75).fit(X, y)

    gb_no_warm_start = GradientBoosting(
        n_iter_no_change=100, max_iter=75, random_state=rng, warm_start=False
    )
    gb_no_warm_start.fit(X, y)

    # Check that both predictors are equal
    _assert_predictor_equal(gb_warm_start, gb_no_warm_start, X)


@pytest.mark.parametrize('GradientBoosting, X, y', [
    (HistGradientBoostingClassifier, X_classification, y_classification),
    (HistGradientBoostingRegressor, X_regression, y_regression)
])
def test_warm_start_max_depth(GradientBoosting, X, y):
    # Test if possible to fit trees of different depth in ensemble.
    gb = GradientBoosting(max_iter=100, min_samples_leaf=1,
                          warm_start=True, max_depth=2)
    gb.fit(X, y)
    gb.set_params(max_iter=110, max_depth=3)
    gb.fit(X, y)

    # First 100 trees have max_depth == 2
    for i in range(100):
        assert_equal(gb._predictors[i][0].get_max_depth(), 2)
    # Last 10 trees have max_depth == 3
    for i in range(1, 11):
        assert_equal(gb._predictors[-i][0].get_max_depth(), 3)


@pytest.mark.parametrize('GradientBoosting, X, y', [
    (HistGradientBoostingClassifier, X_classification, y_classification),
    (HistGradientBoostingRegressor, X_regression, y_regression)
])
def test_warm_start_early_stopping(GradientBoosting, X, y):
    # Make sure that early stopping occurs after a small number of iterations
    # when fitting a second time with warm starting.

    n_iter_no_change = 5
    gb = GradientBoosting(
        n_iter_no_change=n_iter_no_change, max_iter=10000,
        random_state=42, warm_start=True, tol=1e-3
    )
    gb.fit(X, y)
    n_iter_first_fit = gb.n_iter_
    gb.fit(X, y)
    n_iter_second_fit = gb.n_iter_
    assert n_iter_second_fit - n_iter_first_fit < n_iter_no_change


@pytest.mark.parametrize('GradientBoosting, X, y', [
    (HistGradientBoostingClassifier, X_classification, y_classification),
    (HistGradientBoostingRegressor, X_regression, y_regression)
])
def test_warm_start_equal_n_estimators(GradientBoosting, X, y):
    # Test if warm start with equal n_estimators does nothing
    gb_1 = GradientBoosting(max_depth=2)
    gb_1.fit(X, y)

    gb_2 = clone(gb_1)
    gb_2.set_params(max_iter=gb_1.max_iter, warm_start=True)
    gb_2.fit(X, y)

    # Check that both predictors are equal
    _assert_predictor_equal(gb_1, gb_2, X)


@pytest.mark.parametrize('GradientBoosting, X, y', [
    (HistGradientBoostingClassifier, X_classification, y_classification),
    (HistGradientBoostingRegressor, X_regression, y_regression)
])
def test_warm_start_clear(GradientBoosting, X, y):
    # Test if fit clears state.
    gb_1 = GradientBoosting(max_depth=2)
    gb_1.fit(X, y)

    gb_2 = GradientBoosting(max_depth=2, warm_start=True)
    gb_2.fit(X, y)  # inits state
    gb_2.set_params(warm_start=False)
    gb_2.fit(X, y)  # clears old state and equals est

    # Check that both predictors are equal
    _assert_predictor_equal(gb_1, gb_2, X)


@pytest.mark.parametrize('GradientBoosting, X, y', [
    (HistGradientBoostingClassifier, X_classification, y_classification),
    (HistGradientBoostingRegressor, X_regression, y_regression)
])
def test_identical_train_val_split_int(GradientBoosting, X, y):
    # Test if identical splits are generated when random_state is an int.
    gb_1 = GradientBoosting(n_iter_no_change=5, random_state=42)
    gb_1.fit(X, y)
    train_val_seed_1 = gb_1._train_val_split_seed

    gb_2 = GradientBoosting(n_iter_no_change=5, random_state=42,
                            warm_start=True)
    gb_2.fit(X, y)  # inits state
    train_val_seed_2 = gb_2._train_val_split_seed
    gb_2.fit(X, y)  # clears old state and equals est
    train_val_seed_3 = gb_2._train_val_split_seed

    # Check that all seeds are equal
    assert_equal(train_val_seed_1, train_val_seed_2)
    assert_equal(train_val_seed_1, train_val_seed_3)


@pytest.mark.parametrize('GradientBoosting, X, y', [
    (HistGradientBoostingClassifier, X_classification, y_classification),
    (HistGradientBoostingRegressor, X_regression, y_regression)
])
def test_identical_train_val_split_random_state(GradientBoosting, X, y):
    # Test if identical splits are generated when random_state is RandomState
    # instance.
    rng = np.random.RandomState(42)
    gb_1 = GradientBoosting(n_iter_no_change=5, random_state=rng)
    gb_1.fit(X, y)
    train_val_seed_1 = gb_1._train_val_split_seed

    rng = np.random.RandomState(42)
    gb_2 = GradientBoosting(n_iter_no_change=5, random_state=rng,
                            warm_start=True)
    gb_2.fit(X, y)  # inits state
    train_val_seed_2 = gb_2._train_val_split_seed
    gb_2.fit(X, y)  # clears old state and equals est
    train_val_seed_3 = gb_2._train_val_split_seed

    # Check that both seeds are equal
    assert_equal(train_val_seed_1, train_val_seed_2)
    assert_equal(train_val_seed_1, train_val_seed_3)


@pytest.mark.parametrize('GradientBoosting, X, y', [
    (HistGradientBoostingClassifier, X_classification, y_classification),
    (HistGradientBoostingRegressor, X_regression, y_regression)
])
def test_different_train_val_splits(GradientBoosting, X, y):
    # Test if two fits with random_state=None have different splits
    gb_1 = GradientBoosting(n_iter_no_change=5)
    gb_1.fit(X, y)
    train_val_seed_1 = gb_1._train_val_split_seed

    gb_2 = GradientBoosting(n_iter_no_change=5)
    gb_2.fit(X, y)
    train_val_seed_2 = gb_2._train_val_split_seed

    gb_3 = GradientBoosting(n_iter_no_change=5, warm_start=True)
    gb_3.fit(X, y)
    train_val_seed_3 = gb_3._train_val_split_seed

    # Check that all seeds are different
    assert_not_equal(train_val_seed_1, train_val_seed_2)
    assert_not_equal(train_val_seed_1, train_val_seed_3)
    assert_not_equal(train_val_seed_2, train_val_seed_3)


@pytest.mark.parametrize('GradientBoosting, X, y', [
    (HistGradientBoostingClassifier, X_classification, y_classification),
    (HistGradientBoostingRegressor, X_regression, y_regression)
])
def test_identical_small_trainset(GradientBoosting, X, y):
    # Test if two fits with random_state=None have different small trainsets
    gb_1 = GradientBoosting(max_iter=2, max_depth=2, n_iter_no_change=5)
    gb_1.fit(X, y)
    small_trainset_seed_1 = gb_1._small_trainset_seed

    gb_2 = GradientBoosting(max_iter=2, max_depth=2, n_iter_no_change=5)
    gb_2.fit(X, y)
    small_trainset_seed_2 = gb_2._small_trainset_seed

    gb_3 = GradientBoosting(max_iter=2, max_depth=2, n_iter_no_change=5,
                            warm_start=True)
    gb_3.fit(X, y)
    small_trainset_seed_3 = gb_3._small_trainset_seed

    # Check that all seeds are different
    assert_not_equal(small_trainset_seed_1, small_trainset_seed_2)
    assert_not_equal(small_trainset_seed_1, small_trainset_seed_3)
    assert_not_equal(small_trainset_seed_2, small_trainset_seed_3)


@pytest.mark.parametrize('GradientBoosting, X, y', [
    (HistGradientBoostingClassifier, X_classification_large,
     y_classification_large),
    (HistGradientBoostingRegressor, X_regression_large,
     y_regression_large)
])
def test_different_small_trainsets(GradientBoosting, X, y):
    # Test if two fits with random_state=None have different splits
    gb_1 = GradientBoosting(max_iter=2, max_depth=2, n_iter_no_change=5)
    gb_1.fit(X, y)
    small_trainset_seed_1 = gb_1._small_trainset_seed

    gb_2 = GradientBoosting(max_iter=2, max_depth=2, n_iter_no_change=5)
    gb_2.fit(X, y)
    small_trainset_seed_2 = gb_2._small_trainset_seed

    gb_3 = GradientBoosting(max_iter=2, max_depth=2, n_iter_no_change=5,
                            warm_start=True)
    gb_3.fit(X, y)
    small_trainset_seed_3 = gb_3._small_trainset_seed

    # Check that all seeds are different
    assert_not_equal(small_trainset_seed_1, small_trainset_seed_2)
    assert_not_equal(small_trainset_seed_1, small_trainset_seed_3)
    assert_not_equal(small_trainset_seed_2, small_trainset_seed_3)

"""
Testing for the bagging ensemble module (sklearn.ensemble.bagging).
"""

# Author: Gilles Louppe
# License: BSD 3 clause
import re
import warnings
from itertools import product

import numpy as np
import joblib
import pytest

from sklearn.base import BaseEstimator
from sklearn.calibration import CalibratedClassifierCV
from sklearn.utils._testing import assert_array_equal
from sklearn.utils._testing import assert_array_almost_equal
from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.model_selection import GridSearchCV, ParameterGrid
from sklearn.ensemble import BaggingClassifier, BaggingRegressor
from sklearn.linear_model import Perceptron, LogisticRegression
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.svm import SVC, SVR
from sklearn.random_projection import SparseRandomProjection
from sklearn.pipeline import make_pipeline
from sklearn.feature_selection import SelectKBest
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_diabetes, load_iris, make_hastie_10_2
from sklearn.utils import check_random_state
from sklearn.utils._mocking import CheckingClassifier
from sklearn.utils.metadata_routing import RequestType
from sklearn.preprocessing import FunctionTransformer, scale
from itertools import cycle

from scipy.sparse import csc_matrix, csr_matrix

rng = check_random_state(0)

# also load the iris dataset
# and randomly permute it
iris = load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the diabetes dataset
# and randomly permute it
diabetes = load_diabetes()
perm = rng.permutation(diabetes.target.size)
diabetes.data = diabetes.data[perm]
diabetes.target = diabetes.target[perm]


def _weighted(estimator):
    return estimator.set_fit_request(sample_weight=True)


def test_classification():
    # Check classification for various parameter settings.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(
        iris.data, iris.target, random_state=rng
    )
    grid = ParameterGrid(
        {
            "max_samples": [0.5, 1.0],
            "max_features": [1, 4],
            "bootstrap": [True, False],
            "bootstrap_features": [True, False],
        }
    )
    estimators = [
        None,
        DummyClassifier(),
        Perceptron(max_iter=20),
        DecisionTreeClassifier(max_depth=2),
        KNeighborsClassifier(),
        SVC(),
    ]
    # Try different parameter settings with different base classifiers without
    # doing the full cartesian product to keep the test durations low.
    for params, estimator in zip(grid, cycle(estimators)):
        BaggingClassifier(
            estimator=estimator,
            random_state=rng,
            n_estimators=2,
            **params,
        ).fit(X_train, y_train).predict(X_test)


@pytest.mark.parametrize(
    "sparse_format, params, method",
    product(
        [csc_matrix, csr_matrix],
        [
            {
                "max_samples": 0.5,
                "max_features": 2,
                "bootstrap": True,
                "bootstrap_features": True,
            },
            {
                "max_samples": 1.0,
                "max_features": 4,
                "bootstrap": True,
                "bootstrap_features": True,
            },
            {"max_features": 2, "bootstrap": False, "bootstrap_features": True},
            {"max_samples": 0.5, "bootstrap": True, "bootstrap_features": False},
        ],
        ["predict", "predict_proba", "predict_log_proba", "decision_function"],
    ),
)
@pytest.mark.parametrize("use_sample_weight", [False, True])
def test_sparse_classification(sparse_format, params, method, use_sample_weight):
    # Check classification for various parameter settings on sparse input.

    class CustomSVC(SVC):
        """SVC variant that records the nature of the training set"""

        def fit(self, X, y, sample_weight=None):
            super().fit(X, y, sample_weight=sample_weight)
            self.data_type_ = type(X)
            return self

    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(
        scale(iris.data), iris.target, random_state=rng
    )

    X_train_sparse = sparse_format(X_train)
    X_test_sparse = sparse_format(X_test)
    base_estimator = CustomSVC(kernel="linear", decision_function_shape="ovr")
    if use_sample_weight:
        base_estimator = _weighted(base_estimator)
    # Trained on sparse format
    sparse_classifier = BaggingClassifier(
        estimator=base_estimator,
        random_state=1,
        **params,
    )
    if use_sample_weight:
        sample_weight = rng.rand(len(y_train))
        sparse_classifier.fit(X_train_sparse, y_train, sample_weight=sample_weight)
    else:
        sparse_classifier.fit(X_train_sparse, y_train)
    sparse_results = getattr(sparse_classifier, method)(X_test_sparse)

    # Trained on dense format
    dense_classifier = BaggingClassifier(
        estimator=base_estimator,
        random_state=1,
        **params,
    )
    if use_sample_weight:
        dense_classifier.fit(X_train, y_train, sample_weight=sample_weight)
    else:
        dense_classifier.fit(X_train, y_train)
    dense_results = getattr(dense_classifier, method)(X_test)
    assert_array_almost_equal(sparse_results, dense_results)

    sparse_type = type(X_train_sparse)
    types = [i.data_type_ for i in sparse_classifier.estimators_]

    assert all([t == sparse_type for t in types])


def test_regression():
    # Check regression for various parameter settings.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(
        diabetes.data[:50], diabetes.target[:50], random_state=rng
    )
    grid = ParameterGrid(
        {
            "max_samples": [0.5, 1.0],
            "max_features": [0.5, 1.0],
            "bootstrap": [True, False],
            "bootstrap_features": [True, False],
        }
    )

    for estimator in [
        None,
        DummyRegressor(),
        DecisionTreeRegressor(),
        KNeighborsRegressor(),
        SVR(),
    ]:
        for params in grid:
            BaggingRegressor(estimator=estimator, random_state=rng, **params).fit(
                X_train, y_train
            ).predict(X_test)


def test_sparse_regression():
    # Check regression for various parameter settings on sparse input.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(
        diabetes.data[:50], diabetes.target[:50], random_state=rng
    )

    class CustomSVR(SVR):
        """SVC variant that records the nature of the training set"""

        def fit(self, X, y):
            super().fit(X, y)
            self.data_type_ = type(X)
            return self

    parameter_sets = [
        {
            "max_samples": 0.5,
            "max_features": 2,
            "bootstrap": True,
            "bootstrap_features": True,
        },
        {
            "max_samples": 1.0,
            "max_features": 4,
            "bootstrap": True,
            "bootstrap_features": True,
        },
        {"max_features": 2, "bootstrap": False, "bootstrap_features": True},
        {"max_samples": 0.5, "bootstrap": True, "bootstrap_features": False},
    ]

    for sparse_format in [csc_matrix, csr_matrix]:
        X_train_sparse = sparse_format(X_train)
        X_test_sparse = sparse_format(X_test)
        for params in parameter_sets:
            # Trained on sparse format
            sparse_classifier = BaggingRegressor(
                estimator=CustomSVR(), random_state=1, **params
            ).fit(X_train_sparse, y_train)
            sparse_results = sparse_classifier.predict(X_test_sparse)

            # Trained on dense format
            dense_results = (
                BaggingRegressor(estimator=CustomSVR(), random_state=1, **params)
                .fit(X_train, y_train)
                .predict(X_test)
            )

            sparse_type = type(X_train_sparse)
            types = [i.data_type_ for i in sparse_classifier.estimators_]

            assert_array_almost_equal(sparse_results, dense_results)
            assert all([t == sparse_type for t in types])
            assert_array_almost_equal(sparse_results, dense_results)


class DummySizeEstimator(BaseEstimator):
    def fit(self, X, y):
        self.training_size_ = X.shape[0]
        self.training_hash_ = joblib.hash(X)

    def predict(self, X):
        return np.ones(X.shape[0])


def test_bootstrap_samples():
    # Test that bootstrapping samples generate non-perfect base estimators.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(
        diabetes.data, diabetes.target, random_state=rng
    )

    estimator = DecisionTreeRegressor().fit(X_train, y_train)

    # without bootstrap, all trees are perfect on the training set
    ensemble = BaggingRegressor(
        estimator=DecisionTreeRegressor(),
        max_samples=1.0,
        bootstrap=False,
        random_state=rng,
    ).fit(X_train, y_train)

    assert estimator.score(X_train, y_train) == ensemble.score(X_train, y_train)

    # with bootstrap, trees are no longer perfect on the training set
    ensemble = BaggingRegressor(
        estimator=DecisionTreeRegressor(),
        max_samples=1.0,
        bootstrap=True,
        random_state=rng,
    ).fit(X_train, y_train)

    assert estimator.score(X_train, y_train) > ensemble.score(X_train, y_train)

    # check that each sampling correspond to a complete bootstrap resample.
    # the size of each bootstrap should be the same as the input data but
    # the data should be different (checked using the hash of the data).
    ensemble = BaggingRegressor(estimator=DummySizeEstimator(), bootstrap=True).fit(
        X_train, y_train
    )
    training_hash = []
    for estimator in ensemble.estimators_:
        assert estimator.training_size_ == X_train.shape[0]
        training_hash.append(estimator.training_hash_)
    assert len(set(training_hash)) == len(training_hash)


def test_bootstrap_features():
    # Test that bootstrapping features may generate duplicate features.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(
        diabetes.data, diabetes.target, random_state=rng
    )

    ensemble = BaggingRegressor(
        estimator=DecisionTreeRegressor(),
        max_features=1.0,
        bootstrap_features=False,
        random_state=rng,
    ).fit(X_train, y_train)

    for features in ensemble.estimators_features_:
        assert diabetes.data.shape[1] == np.unique(features).shape[0]

    ensemble = BaggingRegressor(
        estimator=DecisionTreeRegressor(),
        max_features=1.0,
        bootstrap_features=True,
        random_state=rng,
    ).fit(X_train, y_train)

    for features in ensemble.estimators_features_:
        assert diabetes.data.shape[1] > np.unique(features).shape[0]


def test_probability():
    # Predict probabilities.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(
        iris.data, iris.target, random_state=rng
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        # Normal case
        ensemble = BaggingClassifier(
            estimator=DecisionTreeClassifier(), random_state=rng
        ).fit(X_train, y_train)

        assert_array_almost_equal(
            np.sum(ensemble.predict_proba(X_test), axis=1), np.ones(len(X_test))
        )

        assert_array_almost_equal(
            ensemble.predict_proba(X_test), np.exp(ensemble.predict_log_proba(X_test))
        )

        # Degenerate case, where some classes are missing
        ensemble = BaggingClassifier(
            estimator=LogisticRegression(), random_state=rng, max_samples=5
        ).fit(X_train, y_train)

        assert_array_almost_equal(
            np.sum(ensemble.predict_proba(X_test), axis=1), np.ones(len(X_test))
        )

        assert_array_almost_equal(
            ensemble.predict_proba(X_test), np.exp(ensemble.predict_log_proba(X_test))
        )


def test_oob_score_classification():
    # Check that oob prediction is a good estimation of the generalization
    # error.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(
        iris.data, iris.target, random_state=rng
    )

    for estimator in [DecisionTreeClassifier(), SVC()]:
        clf = BaggingClassifier(
            estimator=estimator,
            n_estimators=100,
            bootstrap=True,
            oob_score=True,
            random_state=rng,
        ).fit(X_train, y_train)

        test_score = clf.score(X_test, y_test)

        assert abs(test_score - clf.oob_score_) < 0.1

        # Test with few estimators
        warn_msg = (
            "Some inputs do not have OOB scores. This probably means too few "
            "estimators were used to compute any reliable oob estimates."
        )
        with pytest.warns(UserWarning, match=warn_msg):
            clf = BaggingClassifier(
                estimator=estimator,
                n_estimators=1,
                bootstrap=True,
                oob_score=True,
                random_state=rng,
            )
            clf.fit(X_train, y_train)


def test_oob_score_regression():
    # Check that oob prediction is a good estimation of the generalization
    # error.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(
        diabetes.data, diabetes.target, random_state=rng
    )

    clf = BaggingRegressor(
        estimator=DecisionTreeRegressor(),
        n_estimators=50,
        bootstrap=True,
        oob_score=True,
        random_state=rng,
    ).fit(X_train, y_train)

    test_score = clf.score(X_test, y_test)

    assert abs(test_score - clf.oob_score_) < 0.1

    # Test with few estimators
    warn_msg = (
        "Some inputs do not have OOB scores. This probably means too few "
        "estimators were used to compute any reliable oob estimates."
    )
    with pytest.warns(UserWarning, match=warn_msg):
        regr = BaggingRegressor(
            estimator=DecisionTreeRegressor(),
            n_estimators=1,
            bootstrap=True,
            oob_score=True,
            random_state=rng,
        )
        regr.fit(X_train, y_train)


def test_single_estimator():
    # Check singleton ensembles.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(
        diabetes.data, diabetes.target, random_state=rng
    )

    clf1 = BaggingRegressor(
        estimator=KNeighborsRegressor(),
        n_estimators=1,
        bootstrap=False,
        bootstrap_features=False,
        random_state=rng,
    ).fit(X_train, y_train)

    clf2 = KNeighborsRegressor().fit(X_train, y_train)

    assert_array_almost_equal(clf1.predict(X_test), clf2.predict(X_test))


def test_error():
    # Test support of decision_function
    X, y = iris.data, iris.target
    base = DecisionTreeClassifier()
    assert not hasattr(BaggingClassifier(base).fit(X, y), "decision_function")


def test_parallel_classification():
    # Check parallel classification.
    X_train, X_test, y_train, y_test = train_test_split(
        iris.data, iris.target, random_state=0
    )

    ensemble = BaggingClassifier(
        DecisionTreeClassifier(), n_jobs=3, random_state=0
    ).fit(X_train, y_train)

    # predict_proba
    y1 = ensemble.predict_proba(X_test)
    ensemble.set_params(n_jobs=1)
    y2 = ensemble.predict_proba(X_test)
    assert_array_almost_equal(y1, y2)

    ensemble = BaggingClassifier(
        DecisionTreeClassifier(), n_jobs=1, random_state=0
    ).fit(X_train, y_train)

    y3 = ensemble.predict_proba(X_test)
    assert_array_almost_equal(y1, y3)

    # decision_function
    ensemble = BaggingClassifier(
        SVC(decision_function_shape="ovr"), n_jobs=3, random_state=0
    ).fit(X_train, y_train)

    decisions1 = ensemble.decision_function(X_test)
    ensemble.set_params(n_jobs=1)
    decisions2 = ensemble.decision_function(X_test)
    assert_array_almost_equal(decisions1, decisions2)

    ensemble = BaggingClassifier(
        SVC(decision_function_shape="ovr"), n_jobs=1, random_state=0
    ).fit(X_train, y_train)

    decisions3 = ensemble.decision_function(X_test)
    assert_array_almost_equal(decisions1, decisions3)


def test_parallel_regression():
    # Check parallel regression.
    rng = check_random_state(0)

    X_train, X_test, y_train, y_test = train_test_split(
        diabetes.data, diabetes.target, random_state=rng
    )

    ensemble = BaggingRegressor(DecisionTreeRegressor(), n_jobs=3, random_state=0).fit(
        X_train, y_train
    )

    ensemble.set_params(n_jobs=1)
    y1 = ensemble.predict(X_test)
    ensemble.set_params(n_jobs=2)
    y2 = ensemble.predict(X_test)
    assert_array_almost_equal(y1, y2)

    ensemble = BaggingRegressor(DecisionTreeRegressor(), n_jobs=1, random_state=0).fit(
        X_train, y_train
    )

    y3 = ensemble.predict(X_test)
    assert_array_almost_equal(y1, y3)


def test_gridsearch():
    # Check that bagging ensembles can be grid-searched.
    # Transform iris into a binary classification task
    X, y = iris.data, iris.target
    y[y == 2] = 1

    # Grid search with scoring based on decision_function
    parameters = {"n_estimators": (1, 2), "estimator__C": (1, 2)}

    GridSearchCV(BaggingClassifier(SVC()), parameters, scoring="roc_auc").fit(X, y)


def test_estimator():
    # Check estimator and its default values.
    rng = check_random_state(0)

    # Classification
    X_train, X_test, y_train, y_test = train_test_split(
        iris.data, iris.target, random_state=rng
    )

    ensemble = BaggingClassifier(None, n_jobs=3, random_state=0).fit(X_train, y_train)

    assert isinstance(ensemble.estimator_, DecisionTreeClassifier)

    ensemble = BaggingClassifier(
        DecisionTreeClassifier(), n_jobs=3, random_state=0
    ).fit(X_train, y_train)

    assert isinstance(ensemble.estimator_, DecisionTreeClassifier)

    ensemble = BaggingClassifier(Perceptron(), n_jobs=3, random_state=0).fit(
        X_train, y_train
    )

    assert isinstance(ensemble.estimator_, Perceptron)

    # Regression
    X_train, X_test, y_train, y_test = train_test_split(
        diabetes.data, diabetes.target, random_state=rng
    )

    ensemble = BaggingRegressor(None, n_jobs=3, random_state=0).fit(X_train, y_train)

    assert isinstance(ensemble.estimator_, DecisionTreeRegressor)

    ensemble = BaggingRegressor(DecisionTreeRegressor(), n_jobs=3, random_state=0).fit(
        X_train, y_train
    )

    assert isinstance(ensemble.estimator_, DecisionTreeRegressor)

    ensemble = BaggingRegressor(SVR(), n_jobs=3, random_state=0).fit(X_train, y_train)
    assert isinstance(ensemble.estimator_, SVR)


def test_bagging_with_pipeline():
    estimator = BaggingClassifier(
        make_pipeline(SelectKBest(k=1), DecisionTreeClassifier()), max_features=2
    )
    estimator.fit(iris.data, iris.target)
    assert isinstance(estimator[0].steps[-1][1].random_state, int)


class DummyZeroEstimator(BaseEstimator):
    def fit(self, X, y):
        self.classes_ = np.unique(y)
        return self

    def predict(self, X):
        return self.classes_[np.zeros(X.shape[0], dtype=int)]


@pytest.mark.parametrize("estimator_cls", [BaggingClassifier, BaggingRegressor])
def test_bagging_sample_weight_unsupported_but_passed(estimator_cls):
    estimator = estimator_cls(DummyZeroEstimator())
    rng = check_random_state(0)

    estimator.fit(iris.data, iris.target).predict(iris.data)
    msg = (
        "fit got unexpected argument(s) {'sample_weight'}, which are not requested "
        "metadata in any object."
    )
    with pytest.raises(TypeError, match=re.escape(msg)):
        estimator.fit(
            iris.data,
            iris.target,
            sample_weight=rng.randint(10, size=(iris.data.shape[0])),
        )


@pytest.mark.parametrize("estimator_cls", [BaggingClassifier, BaggingRegressor])
def test_sample_weight_none_and_base_estimator_supports_it_directly(estimator_cls):
    # Using iris even for regression, which is fine for the purpose of this test
    X, y = iris.data, iris.target
    # sample bagging is achieved through the use of sample weights because the
    # base estimator supports it, thus expect sample weights
    base_estimator = CheckingClassifier(expected_sample_weight=True)
    estimator = estimator_cls(base_estimator)

    # there are no warnings and no errors
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        estimator.fit(X, y)


@pytest.mark.parametrize("estimator_cls", [BaggingClassifier, BaggingRegressor])
def test_sample_weight_none_and_base_estimator_supports_it_indirectly(estimator_cls):
    # Even though "sample_weight" is not in the fit signature, it is still
    # passed. This test exists because bagging used to determine if
    # sample_weight should be passed by inspecting the fit signature.
    from sklearn.utils.metadata_routing import RequestType

    class MyEstimator(BaseEstimator):
        # This attribute indicates that the estimator indeed supports
        # sample_weight
        __metadata_request__fit = {"sample_weight": RequestType.ERROR_IF_PASSED}

        def fit(self, X, y, **kwargs):
            assert kwargs.get("sample_weight") is not None
            return self

        def predict(self, X):
            return np.ones(len(X))

    # Using iris even for regression, which is fine for the purpose of this test
    X, y = iris.data, iris.target
    # sample bagging is achieved through the use of sample weights because the
    # base estimator supports it
    base_estimator = MyEstimator()
    estimator = estimator_cls(base_estimator)

    # there are no warnings and no errors
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        estimator.fit(X, y)


@pytest.mark.parametrize("estimator_cls", [BaggingClassifier, BaggingRegressor])
def test_sample_weight_none_and_base_estimator_does_not_support_it(estimator_cls):
    # sample bagging is achieved through indexing because the base estimator
    # does not support it

    # Using iris even for regression, which is fine for the purpose of this test
    X, y = iris.data, iris.target
    n_samples = len(y)

    class MyEstimator(BaseEstimator):
        def fit(self, X, y, **fit_params):
            # check that the indexing method was used to achieve subsampling,
            # i.e. no sample_weight are passed but instead fewer samples
            assert "sample_weight" not in fit_params
            assert len(y) == 0.5 * n_samples
            return self

        def predict(self, X):
            return np.ones(len(X))

    base_estimator = MyEstimator()
    estimator = estimator_cls(base_estimator, max_samples=0.5)

    # there are no warnings and no errors
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        estimator.fit(X, y)


@pytest.mark.parametrize("estimator_cls", [BaggingClassifier, BaggingRegressor])
def test_sample_weight_not_none_and_base_estimator_supports_it_directly(estimator_cls):
    # sample bagging is achieved through the use of sample weights because the
    # base estimator supports it (but has to request it)

    # Using iris even for regression, which is fine for the purpose of this test
    X, y = iris.data, iris.target
    sample_weight = np.ones_like(y)

    base_estimator = _weighted(CheckingClassifier(expected_sample_weight=False))
    estimator = estimator_cls(base_estimator)

    # there are no warnings and no errors
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        estimator.fit(X, y, sample_weight=sample_weight)


@pytest.mark.parametrize("estimator_cls", [BaggingClassifier, BaggingRegressor])
def test_sample_weight_not_none_and_base_estimator_supports_it_indirectly(
    estimator_cls,
):
    # sample bagging is achieved through the use of sample weights because the
    # base estimator supports it indirectly (but has to request it)

    class MyEstimator(BaseEstimator):
        __metadata_request__fit = {"sample_weight": RequestType.ERROR_IF_PASSED}

        def fit(self, X, y, **kwargs):
            sw = kwargs.get("sample_weight")
            # to test that subsampling via sample_weight is used, not only check
            # that sample_weight was passed but also that it's not just 1s (the
            # original sample_weight).
            assert sw is not None
            assert not np.allclose(sw, 1.0)
            return self

        def predict(self, X):
            return np.ones(len(X))

    # Using iris even for regression, which is fine for the purpose of this test
    X, y = iris.data, iris.target
    sample_weight = np.ones_like(y)
    # sample bagging is achieved through the use of sample weights because the
    # base estimator supports it (but has to request it)
    base_estimator = _weighted(MyEstimator())
    estimator = estimator_cls(base_estimator)

    # there are no warnings and no errors
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        estimator.fit(X, y, sample_weight=sample_weight)


@pytest.mark.parametrize("estimator_cls", [BaggingClassifier, BaggingRegressor])
@pytest.mark.parametrize("use_sw", [True, False])
@pytest.mark.parametrize("base_estimator_supports_sw", [True, False])
def test_sample_weight_supported_and_used_and_base_estimator_is_metaestimator(
    estimator_cls, use_sw, base_estimator_supports_sw
):
    # Check bagging when the base estimator is itself a meta estimator. The
    # logic is that as long as that meta estimator supports SW, it is assumed
    # that it works correctly and passes it on to its own base estimator if
    # necessary. We need to write custom estimator classes to check that if
    # supported, the sample weight method is used for bagging, not the indexing
    # method, and vice versa.

    # Using iris even for regression, which is fine for the purpose of this test
    X, y = iris.data, iris.target
    num_samples = len(y)
    sample_weight = np.ones_like(y)

    from sklearn.base import MetaEstimatorMixin
    from sklearn.utils.metadata_routing import (
        MetadataRouter,
        MethodMapping,
        process_routing,
    )

    class MetaEstimatorWithoutSW(BaseEstimator, MetaEstimatorMixin):
        def __init__(self, estimator):
            self.estimator = estimator

        def check(self, X, y, **fit_params):
            assert "sample_weight" not in fit_params
            assert len(y) == 0.5 * num_samples  # check that indexing was used

        def fit(self, X, y, **fit_params):
            self.check(X, y, **fit_params)
            routed_params = process_routing(
                obj=self,
                method="fit",
                other_params=fit_params,
            )
            self.estimator.fit(X, y, **routed_params.estimator.fit)
            return self

        def predict(self, X):
            return np.ones(len(X))

        def get_metadata_routing(self):
            router = (
                MetadataRouter(owner=self.__class__.__name__)
                .add(
                    estimator=self.estimator,
                    method_mapping=MethodMapping().add(callee="fit", caller="fit"),
                )
                .warn_on(child="estimator", method="fit", params=["sample_weight"])
            )
            return router

    class MetaEstimatorWithSW(MetaEstimatorWithoutSW):
        def check(self, X, y, **fit_params):
            sw = fit_params.get("sample_weight")
            assert sw is not None
            # ensure that indexing was not used
            assert len(sw) == num_samples
            # ensure that the SW method was used for bagging by checking that
            # the passed SW are not just the original one, which was all 1s
            assert not np.allclose(sw, 1.0)

        def fit(self, X, y, sample_weight=None, **fit_params):
            self.check(X, y, sample_weight=sample_weight, **fit_params)
            return super().fit(X, y, sample_weight=sample_weight, **fit_params)

    class BaseEstimatorWithoutSW(BaseEstimator):
        def fit(self, X, y):
            assert len(y) == 0.5 * num_samples  # check that indexing was used
            return self

        def predict(self, X):
            return np.ones(len(X))

    class BaseEstimatorWithSW(BaseEstimatorWithoutSW):
        def fit(self, X, y, sample_weight=None):
            assert sample_weight is not None
            # ensure that indexing was not used
            assert len(sample_weight) == num_samples
            # ensure that the SW method was used for bagging by checking that
            # the passed SW are not just the original one, which was all 1s
            assert not np.allclose(sample_weight, 1.0)
            return self

    if base_estimator_supports_sw:
        base_base_estimator = BaseEstimatorWithSW().set_fit_request(sample_weight=True)
        # TODO: should line below warn when not explicitly requesting SW for fit?
        base_estimator = MetaEstimatorWithSW(base_base_estimator)
    else:
        base_base_estimator = BaseEstimatorWithoutSW()
        base_estimator = MetaEstimatorWithoutSW(base_base_estimator)

    estimator = estimator_cls(base_estimator, max_samples=0.5)

    if use_sw and not base_estimator_supports_sw:
        # this case is covered by metadata routing, not specific to bagging
        with pytest.raises(TypeError):
            estimator.fit(X, y, sample_weight=sample_weight)
        return

    # there are no warnings and no errors
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        if use_sw:
            estimator.fit(X, y, sample_weight=sample_weight)
        else:
            estimator.fit(X, y)


def test_warm_start(random_state=42):
    # Test if fitting incrementally with warm start gives a forest of the
    # right size and the same results as a normal fit.
    X, y = make_hastie_10_2(n_samples=20, random_state=1)

    clf_ws = None
    for n_estimators in [5, 10]:
        if clf_ws is None:
            clf_ws = BaggingClassifier(
                n_estimators=n_estimators, random_state=random_state, warm_start=True
            )
        else:
            clf_ws.set_params(n_estimators=n_estimators)
        clf_ws.fit(X, y)
        assert len(clf_ws) == n_estimators

    clf_no_ws = BaggingClassifier(
        n_estimators=10, random_state=random_state, warm_start=False
    )
    clf_no_ws.fit(X, y)

    assert set([tree.random_state for tree in clf_ws]) == set(
        [tree.random_state for tree in clf_no_ws]
    )


def test_warm_start_smaller_n_estimators():
    # Test if warm start'ed second fit with smaller n_estimators raises error.
    X, y = make_hastie_10_2(n_samples=20, random_state=1)
    clf = BaggingClassifier(n_estimators=5, warm_start=True)
    clf.fit(X, y)
    clf.set_params(n_estimators=4)
    with pytest.raises(ValueError):
        clf.fit(X, y)


def test_warm_start_equal_n_estimators():
    # Test that nothing happens when fitting without increasing n_estimators
    X, y = make_hastie_10_2(n_samples=20, random_state=1)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=43)

    clf = BaggingClassifier(n_estimators=5, warm_start=True, random_state=83)
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    # modify X to nonsense values, this should not change anything
    X_train += 1.0

    warn_msg = "Warm-start fitting without increasing n_estimators does not"
    with pytest.warns(UserWarning, match=warn_msg):
        clf.fit(X_train, y_train)
    assert_array_equal(y_pred, clf.predict(X_test))


def test_warm_start_equivalence():
    # warm started classifier with 5+5 estimators should be equivalent to
    # one classifier with 10 estimators
    X, y = make_hastie_10_2(n_samples=20, random_state=1)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=43)

    clf_ws = BaggingClassifier(n_estimators=5, warm_start=True, random_state=3141)
    clf_ws.fit(X_train, y_train)
    clf_ws.set_params(n_estimators=10)
    clf_ws.fit(X_train, y_train)
    y1 = clf_ws.predict(X_test)

    clf = BaggingClassifier(n_estimators=10, warm_start=False, random_state=3141)
    clf.fit(X_train, y_train)
    y2 = clf.predict(X_test)

    assert_array_almost_equal(y1, y2)


def test_warm_start_with_oob_score_fails():
    # Check using oob_score and warm_start simultaneously fails
    X, y = make_hastie_10_2(n_samples=20, random_state=1)
    clf = BaggingClassifier(n_estimators=5, warm_start=True, oob_score=True)
    with pytest.raises(ValueError):
        clf.fit(X, y)


def test_oob_score_removed_on_warm_start():
    X, y = make_hastie_10_2(n_samples=100, random_state=1)

    clf = BaggingClassifier(n_estimators=5, oob_score=True)
    clf.fit(X, y)

    clf.set_params(warm_start=True, oob_score=False, n_estimators=10)
    clf.fit(X, y)

    with pytest.raises(AttributeError):
        getattr(clf, "oob_score_")


def test_oob_score_consistency():
    # Make sure OOB scores are identical when random_state, estimator, and
    # training data are fixed and fitting is done twice
    X, y = make_hastie_10_2(n_samples=200, random_state=1)
    bagging = BaggingClassifier(
        KNeighborsClassifier(),
        max_samples=0.5,
        max_features=0.5,
        oob_score=True,
        random_state=1,
    )
    assert bagging.fit(X, y).oob_score_ == bagging.fit(X, y).oob_score_


def test_estimators_samples():
    # Check that format of estimators_samples_ is correct and that results
    # generated at fit time can be identically reproduced at a later time
    # using data saved in object attributes.
    X, y = make_hastie_10_2(n_samples=200, random_state=1)
    bagging = BaggingClassifier(
        LogisticRegression(),
        max_samples=0.5,
        max_features=0.5,
        random_state=1,
        bootstrap=False,
    )
    bagging.fit(X, y)

    # Get relevant attributes
    estimators_samples = bagging.estimators_samples_
    estimators_features = bagging.estimators_features_
    estimators = bagging.estimators_

    # Test for correct formatting
    assert len(estimators_samples) == len(estimators)
    assert len(estimators_samples[0]) == len(X) // 2
    assert estimators_samples[0].dtype.kind == "i"

    # Re-fit single estimator to test for consistent sampling
    estimator_index = 0
    estimator_samples = estimators_samples[estimator_index]
    estimator_features = estimators_features[estimator_index]
    estimator = estimators[estimator_index]

    X_train = (X[estimator_samples])[:, estimator_features]
    y_train = y[estimator_samples]

    orig_coefs = estimator.coef_
    estimator.fit(X_train, y_train)
    new_coefs = estimator.coef_

    assert_array_almost_equal(orig_coefs, new_coefs)


def test_estimators_samples_deterministic():
    # This test is a regression test to check that with a random step
    # (e.g. SparseRandomProjection) and a given random state, the results
    # generated at fit time can be identically reproduced at a later time using
    # data saved in object attributes. Check issue #9524 for full discussion.

    iris = load_iris()
    X, y = iris.data, iris.target

    base_pipeline = make_pipeline(
        SparseRandomProjection(n_components=2), LogisticRegression()
    )
    clf = BaggingClassifier(estimator=base_pipeline, max_samples=0.5, random_state=0)
    clf.fit(X, y)
    pipeline_estimator_coef = clf.estimators_[0].steps[-1][1].coef_.copy()

    estimator = clf.estimators_[0]
    estimator_sample = clf.estimators_samples_[0]
    estimator_feature = clf.estimators_features_[0]

    X_train = (X[estimator_sample])[:, estimator_feature]
    y_train = y[estimator_sample]

    estimator.fit(X_train, y_train)
    assert_array_equal(estimator.steps[-1][1].coef_, pipeline_estimator_coef)


def test_max_samples_consistency():
    # Make sure validated max_samples and original max_samples are identical
    # when valid integer max_samples supplied by user
    max_samples = 100
    X, y = make_hastie_10_2(n_samples=2 * max_samples, random_state=1)
    bagging = BaggingClassifier(
        KNeighborsClassifier(),
        max_samples=max_samples,
        max_features=0.5,
        random_state=1,
    )
    bagging.fit(X, y)
    assert bagging._max_samples == max_samples


def test_set_oob_score_label_encoding():
    # Make sure the oob_score doesn't change when the labels change
    # See: https://github.com/scikit-learn/scikit-learn/issues/8933
    random_state = 5
    X = [[-1], [0], [1]] * 5
    Y1 = ["A", "B", "C"] * 5
    Y2 = [-1, 0, 1] * 5
    Y3 = [0, 1, 2] * 5
    x1 = (
        BaggingClassifier(oob_score=True, random_state=random_state)
        .fit(X, Y1)
        .oob_score_
    )
    x2 = (
        BaggingClassifier(oob_score=True, random_state=random_state)
        .fit(X, Y2)
        .oob_score_
    )
    x3 = (
        BaggingClassifier(oob_score=True, random_state=random_state)
        .fit(X, Y3)
        .oob_score_
    )
    assert [x1, x2] == [x3, x3]


def replace(X):
    X = X.astype("float", copy=True)
    X[~np.isfinite(X)] = 0
    return X


def test_bagging_regressor_with_missing_inputs():
    # Check that BaggingRegressor can accept X with missing/infinite data
    X = np.array(
        [
            [1, 3, 5],
            [2, None, 6],
            [2, np.nan, 6],
            [2, np.inf, 6],
            [2, np.NINF, 6],
        ]
    )
    y_values = [
        np.array([2, 3, 3, 3, 3]),
        np.array(
            [
                [2, 1, 9],
                [3, 6, 8],
                [3, 6, 8],
                [3, 6, 8],
                [3, 6, 8],
            ]
        ),
    ]
    for y in y_values:
        regressor = DecisionTreeRegressor()
        pipeline = make_pipeline(FunctionTransformer(replace), regressor)
        pipeline.fit(X, y).predict(X)
        bagging_regressor = BaggingRegressor(pipeline)
        y_hat = bagging_regressor.fit(X, y).predict(X)
        assert y.shape == y_hat.shape

        # Verify that exceptions can be raised by wrapper regressor
        regressor = DecisionTreeRegressor()
        pipeline = make_pipeline(regressor)
        with pytest.raises(ValueError):
            pipeline.fit(X, y)
        bagging_regressor = BaggingRegressor(pipeline)
        with pytest.raises(ValueError):
            bagging_regressor.fit(X, y)


def test_bagging_classifier_with_missing_inputs():
    # Check that BaggingClassifier can accept X with missing/infinite data
    X = np.array(
        [
            [1, 3, 5],
            [2, None, 6],
            [2, np.nan, 6],
            [2, np.inf, 6],
            [2, np.NINF, 6],
        ]
    )
    y = np.array([3, 6, 6, 6, 6])
    classifier = DecisionTreeClassifier()
    pipeline = make_pipeline(FunctionTransformer(replace), classifier)
    pipeline.fit(X, y).predict(X)
    bagging_classifier = BaggingClassifier(pipeline)
    bagging_classifier.fit(X, y)
    y_hat = bagging_classifier.predict(X)
    assert y.shape == y_hat.shape
    bagging_classifier.predict_log_proba(X)
    bagging_classifier.predict_proba(X)

    # Verify that exceptions can be raised by wrapper classifier
    classifier = DecisionTreeClassifier()
    pipeline = make_pipeline(classifier)
    with pytest.raises(ValueError):
        pipeline.fit(X, y)
    bagging_classifier = BaggingClassifier(pipeline)
    with pytest.raises(ValueError):
        bagging_classifier.fit(X, y)


def test_bagging_small_max_features():
    # Check that Bagging estimator can accept low fractional max_features

    X = np.array([[1, 2], [3, 4]])
    y = np.array([1, 0])

    bagging = BaggingClassifier(LogisticRegression(), max_features=0.3, random_state=1)
    bagging.fit(X, y)


def test_bagging_get_estimators_indices():
    # Check that Bagging estimator can generate sample indices properly
    # Non-regression test for:
    # https://github.com/scikit-learn/scikit-learn/issues/16436

    rng = np.random.RandomState(0)
    X = rng.randn(13, 4)
    y = np.arange(13)

    class MyEstimator(DecisionTreeRegressor):
        """An estimator which stores y indices information at fit."""

        def fit(self, X, y):
            self._sample_indices = y

    clf = BaggingRegressor(estimator=MyEstimator(), n_estimators=1, random_state=0)
    clf.fit(X, y)

    assert_array_equal(clf.estimators_[0]._sample_indices, clf.estimators_samples_[0])


# TODO(1.4): remove in 1.4
@pytest.mark.parametrize(
    "Bagging, Estimator",
    [
        (BaggingClassifier, DecisionTreeClassifier),
        (BaggingRegressor, DecisionTreeRegressor),
    ],
)
def test_base_estimator_argument_deprecated(Bagging, Estimator):
    X = np.array([[1, 2], [3, 4]])
    y = np.array([1, 0])
    model = Bagging(base_estimator=Estimator(), n_estimators=10)

    warn_msg = (
        "`base_estimator` was renamed to `estimator` in version 1.2 and "
        "will be removed in 1.4."
    )
    with pytest.warns(FutureWarning, match=warn_msg):
        model.fit(X, y)


# TODO(1.4): remove in 1.4
@pytest.mark.parametrize(
    "Bagging",
    [BaggingClassifier, BaggingClassifier],
)
def test_base_estimator_property_deprecated(Bagging):
    X = np.array([[1, 2], [3, 4]])
    y = np.array([1, 0])
    model = Bagging()
    model.fit(X, y)

    warn_msg = (
        "Attribute `base_estimator_` was deprecated in version 1.2 and "
        "will be removed in 1.4. Use `estimator_` instead."
    )
    with pytest.warns(FutureWarning, match=warn_msg):
        model.base_estimator_


# TODO(1.4): remove
def test_deprecated_base_estimator_has_decision_function():
    """Check that `BaggingClassifier` delegate to classifier with
    `decision_function`."""
    iris = load_iris()
    X, y = iris.data, iris.target
    clf = BaggingClassifier(base_estimator=SVC())
    assert hasattr(clf, "decision_function")
    warn_msg = (
        "`base_estimator` was renamed to `estimator` in version 1.2 and "
        "will be removed in 1.4."
    )
    with pytest.warns(FutureWarning, match=warn_msg):
        y_decision = clf.fit(X, y).decision_function(X)
    assert y_decision.shape == (150, 3)


# TODO(1.4): remove in 1.4
@pytest.mark.parametrize(
    "Bagging, Estimator",
    [
        (BaggingClassifier, DecisionTreeClassifier),
        (BaggingRegressor, DecisionTreeRegressor),
    ],
)
def test_base_estimator_argument_missing_fit_request_warns(Bagging, Estimator):
    # When using the deprecated base_estimator, there should still be a warning
    # when not setting the fit request
    X = np.array([[1, 2], [3, 4]])
    y = np.array([1, 0])
    sample_weight = np.array([1.0, 1.0])
    model = Bagging(base_estimator=Estimator(), n_estimators=10)

    warn_msg = (
        "You are passing metadata for which the request values are not "
        "explicitly set: sample_weight"
    )
    with pytest.warns(FutureWarning, match=warn_msg):
        model.fit(X, y, sample_weight=sample_weight)

    # When setting the fit request, the warning disappears. To check that, it's
    # not possible to check for the absence of FutureWarnings, since there will
    # still be the deprecation for base_estimator. Therefore, catch all warnings
    # and check the warning message.
    model = Bagging(base_estimator=_weighted(Estimator()), n_estimators=10)
    with pytest.warns(FutureWarning) as record:
        model.fit(X, y, sample_weight=sample_weight)
        assert not any(warn_msg in w.message.args[0] for w in record.list)


@pytest.mark.parametrize("estimator_cls", [BaggingClassifier, BaggingRegressor])
def test_base_estimator_with_aliased_sample_weight_raises(estimator_cls):
    # Aliasing sample weight on the base estimator is not supported because it
    # means that there are potentially two different sample weights and it's
    # unclear what to do with them. Since aliasing was not possible before
    # introducing metadata routing, this does not preclude something that was
    # possible before.

    # Using iris even for regression, which is fine for the purpose of this test
    X, y = iris.data, iris.target
    sample_weight = np.ones_like(y)

    estimator = SVC().set_fit_request(sample_weight="my_sample_weight")
    estimator = estimator_cls(estimator, max_samples=0.5)

    msg = "Aliasing sample_weight is not allowed when using Bagging"
    with pytest.raises(ValueError, match=msg):
        estimator.fit(X, y, my_sample_weight=sample_weight)


@pytest.mark.parametrize("estimator_cls", [BaggingClassifier, BaggingRegressor])
def test_base_estimator_metaestimator_with_aliased_sample_weight_raises(estimator_cls):
    # Same reasoning as previous test but using metaestimator as base estimator

    # Using iris even for regression, which is fine for the purpose of this test
    X, y = iris.data, iris.target
    sample_weight = np.ones_like(y)

    estimator = CalibratedClassifierCV(SVC()).set_fit_request(
        sample_weight="my_sample_weight"
    )
    estimator = estimator_cls(estimator, max_samples=0.5)

    msg = "Aliasing sample_weight is not allowed when using Bagging"
    with pytest.raises(ValueError, match=msg):
        estimator.fit(X, y, my_sample_weight=sample_weight)


@pytest.mark.parametrize("estimator_cls", [BaggingClassifier, BaggingRegressor])
def test_base_base_estimator_metaestimator_with_aliased_sample_weight_raises(
    estimator_cls,
):
    # Same reasoning as previous test but using metaestimator as base estimator
    # whose own base estimator aliases sample_weight.

    # Using iris even for regression, which is fine for the purpose of this test
    X, y = iris.data, iris.target
    sample_weight = np.ones_like(y)

    estimator = CalibratedClassifierCV(
        SVC().set_fit_request(sample_weight="my_sample_weight")
    )
    estimator = estimator_cls(estimator, max_samples=0.5)

    msg = "Aliasing sample_weight is not allowed when using Bagging"
    with pytest.raises(ValueError, match=msg):
        estimator.fit(X, y, my_sample_weight=sample_weight)

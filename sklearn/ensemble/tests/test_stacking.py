"""Test the stacking classifier and regressor."""

# Authors: Guillaume Lemaitre <g.lemaitre58@gmail.com>
# License: BSD 3 clause

import pytest
import numpy as np
from numpy.testing import assert_array_equal
import scipy.sparse as sparse

from sklearn.base import BaseEstimator
from sklearn.base import ClassifierMixin
from sklearn.base import RegressorMixin
from sklearn.base import clone

from sklearn.exceptions import ConvergenceWarning

from sklearn.datasets import load_iris
from sklearn.datasets import load_diabetes
from sklearn.datasets import load_breast_cancer
from sklearn.datasets import make_regression
from sklearn.datasets import make_classification
from sklearn.datasets import make_multilabel_classification

from sklearn.dummy import DummyClassifier
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
from sklearn.svm import LinearSVC
from sklearn.svm import LinearSVR
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import scale

from sklearn.ensemble import StackingClassifier
from sklearn.ensemble import StackingRegressor

from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import KFold

from sklearn.utils._mocking import CheckingClassifier
from sklearn.utils._testing import assert_allclose
from sklearn.utils._testing import assert_allclose_dense_sparse
from sklearn.utils._testing import ignore_warnings

from sklearn.exceptions import NotFittedError

from unittest.mock import Mock

diabetes = load_diabetes()
X_diabetes, y_diabetes = diabetes.data, diabetes.target
iris = load_iris()
X_iris, y_iris = iris.data, iris.target
X_multilabel_r, y_multilabel_r = make_multilabel_classification(
    n_classes=3, random_state=42
)
X_binary, y_binary = make_classification(n_classes=2, random_state=42)


@pytest.mark.parametrize(
    "cv", [3, StratifiedKFold(n_splits=3, shuffle=True, random_state=42)]
)
@pytest.mark.parametrize(
    "final_estimator", [None, RandomForestClassifier(random_state=42)]
)
@pytest.mark.parametrize("passthrough", [False, True])
def test_stacking_classifier_iris(cv, final_estimator, passthrough):
    # prescale the data to avoid convergence warning without using a pipeline
    # for later assert
    X_train, X_test, y_train, y_test = train_test_split(
        scale(X_iris), y_iris, stratify=y_iris, random_state=42
    )
    estimators = [("lr", LogisticRegression()), ("svc", LinearSVC())]
    clf = StackingClassifier(
        estimators=estimators,
        final_estimator=final_estimator,
        cv=cv,
        passthrough=passthrough,
    )
    clf.fit(X_train, y_train)
    clf.predict(X_test)
    clf.predict_proba(X_test)
    assert clf.score(X_test, y_test) > 0.8

    X_trans = clf.transform(X_test)
    expected_column_count = 10 if passthrough else 6
    assert X_trans.shape[1] == expected_column_count
    if passthrough:
        assert_allclose(X_test, X_trans[:, -4:])

    clf.set_params(lr="drop")
    clf.fit(X_train, y_train)
    clf.predict(X_test)
    clf.predict_proba(X_test)
    if final_estimator is None:
        # LogisticRegression has decision_function method
        clf.decision_function(X_test)

    X_trans = clf.transform(X_test)
    expected_column_count_drop = 7 if passthrough else 3
    assert X_trans.shape[1] == expected_column_count_drop
    if passthrough:
        assert_allclose(X_test, X_trans[:, -4:])


def test_stacking_classifier_drop_column_binary_classification():
    # check that a column is dropped in binary classification
    X, y = load_breast_cancer(return_X_y=True)
    X_train, X_test, y_train, _ = train_test_split(
        scale(X), y, stratify=y, random_state=42
    )

    # both classifiers implement 'predict_proba' and will both drop one column
    estimators = [
        ("lr", LogisticRegression()),
        ("rf", RandomForestClassifier(random_state=42)),
    ]
    clf = StackingClassifier(estimators=estimators, cv=3)

    clf.fit(X_train, y_train)
    X_trans = clf.transform(X_test)
    assert X_trans.shape[1] == 2

    # LinearSVC does not implement 'predict_proba' and will not drop one column
    estimators = [("lr", LogisticRegression()), ("svc", LinearSVC())]
    clf.set_params(estimators=estimators)

    clf.fit(X_train, y_train)
    X_trans = clf.transform(X_test)
    assert X_trans.shape[1] == 2


def test_stacking_classifier_drop_estimator():
    # prescale the data to avoid convergence warning without using a pipeline
    # for later assert
    X_train, X_test, y_train, _ = train_test_split(
        scale(X_iris), y_iris, stratify=y_iris, random_state=42
    )
    estimators = [("lr", "drop"), ("svc", LinearSVC(random_state=0))]
    rf = RandomForestClassifier(n_estimators=10, random_state=42)
    clf = StackingClassifier(
        estimators=[("svc", LinearSVC(random_state=0))], final_estimator=rf, cv=5
    )
    clf_drop = StackingClassifier(estimators=estimators, final_estimator=rf, cv=5)

    clf.fit(X_train, y_train)
    clf_drop.fit(X_train, y_train)
    assert_allclose(clf.predict(X_test), clf_drop.predict(X_test))
    assert_allclose(clf.predict_proba(X_test), clf_drop.predict_proba(X_test))
    assert_allclose(clf.transform(X_test), clf_drop.transform(X_test))


def test_stacking_regressor_drop_estimator():
    # prescale the data to avoid convergence warning without using a pipeline
    # for later assert
    X_train, X_test, y_train, _ = train_test_split(
        scale(X_diabetes), y_diabetes, random_state=42
    )
    estimators = [("lr", "drop"), ("svr", LinearSVR(random_state=0))]
    rf = RandomForestRegressor(n_estimators=10, random_state=42)
    reg = StackingRegressor(
        estimators=[("svr", LinearSVR(random_state=0))], final_estimator=rf, cv=5
    )
    reg_drop = StackingRegressor(estimators=estimators, final_estimator=rf, cv=5)

    reg.fit(X_train, y_train)
    reg_drop.fit(X_train, y_train)
    assert_allclose(reg.predict(X_test), reg_drop.predict(X_test))
    assert_allclose(reg.transform(X_test), reg_drop.transform(X_test))


@pytest.mark.parametrize("cv", [3, KFold(n_splits=3, shuffle=True, random_state=42)])
@pytest.mark.parametrize(
    "final_estimator, predict_params",
    [
        (None, {}),
        (RandomForestRegressor(random_state=42), {}),
        (DummyRegressor(), {"return_std": True}),
    ],
)
@pytest.mark.parametrize("passthrough", [False, True])
def test_stacking_regressor_diabetes(cv, final_estimator, predict_params, passthrough):
    # prescale the data to avoid convergence warning without using a pipeline
    # for later assert
    X_train, X_test, y_train, _ = train_test_split(
        scale(X_diabetes), y_diabetes, random_state=42
    )
    estimators = [("lr", LinearRegression()), ("svr", LinearSVR())]
    reg = StackingRegressor(
        estimators=estimators,
        final_estimator=final_estimator,
        cv=cv,
        passthrough=passthrough,
    )
    reg.fit(X_train, y_train)
    result = reg.predict(X_test, **predict_params)
    expected_result_length = 2 if predict_params else 1
    if predict_params:
        assert len(result) == expected_result_length

    X_trans = reg.transform(X_test)
    expected_column_count = 12 if passthrough else 2
    assert X_trans.shape[1] == expected_column_count
    if passthrough:
        assert_allclose(X_test, X_trans[:, -10:])

    reg.set_params(lr="drop")
    reg.fit(X_train, y_train)
    reg.predict(X_test)

    X_trans = reg.transform(X_test)
    expected_column_count_drop = 11 if passthrough else 1
    assert X_trans.shape[1] == expected_column_count_drop
    if passthrough:
        assert_allclose(X_test, X_trans[:, -10:])


@pytest.mark.parametrize("fmt", ["csc", "csr", "coo"])
def test_stacking_regressor_sparse_passthrough(fmt):
    # Check passthrough behavior on a sparse X matrix
    X_train, X_test, y_train, _ = train_test_split(
        sparse.coo_matrix(scale(X_diabetes)).asformat(fmt), y_diabetes, random_state=42
    )
    estimators = [("lr", LinearRegression()), ("svr", LinearSVR())]
    rf = RandomForestRegressor(n_estimators=10, random_state=42)
    clf = StackingRegressor(
        estimators=estimators, final_estimator=rf, cv=5, passthrough=True
    )
    clf.fit(X_train, y_train)
    X_trans = clf.transform(X_test)
    assert_allclose_dense_sparse(X_test, X_trans[:, -10:])
    assert sparse.issparse(X_trans)
    assert X_test.format == X_trans.format


@pytest.mark.parametrize("fmt", ["csc", "csr", "coo"])
def test_stacking_classifier_sparse_passthrough(fmt):
    # Check passthrough behavior on a sparse X matrix
    X_train, X_test, y_train, _ = train_test_split(
        sparse.coo_matrix(scale(X_iris)).asformat(fmt), y_iris, random_state=42
    )
    estimators = [("lr", LogisticRegression()), ("svc", LinearSVC())]
    rf = RandomForestClassifier(n_estimators=10, random_state=42)
    clf = StackingClassifier(
        estimators=estimators, final_estimator=rf, cv=5, passthrough=True
    )
    clf.fit(X_train, y_train)
    X_trans = clf.transform(X_test)
    assert_allclose_dense_sparse(X_test, X_trans[:, -4:])
    assert sparse.issparse(X_trans)
    assert X_test.format == X_trans.format


def test_stacking_classifier_drop_binary_prob():
    # check that classifier will drop one of the probability column for
    # binary classification problem

    # Select only the 2 first classes
    X_, y_ = scale(X_iris[:100]), y_iris[:100]

    estimators = [("lr", LogisticRegression()), ("rf", RandomForestClassifier())]
    clf = StackingClassifier(estimators=estimators)
    clf.fit(X_, y_)
    X_meta = clf.transform(X_)
    assert X_meta.shape[1] == 2


class NoWeightRegressor(RegressorMixin, BaseEstimator):
    def fit(self, X, y):
        self.reg = DummyRegressor()
        return self.reg.fit(X, y)

    def predict(self, X):
        return np.ones(X.shape[0])


class NoWeightClassifier(ClassifierMixin, BaseEstimator):
    def fit(self, X, y):
        self.clf = DummyClassifier(strategy="stratified")
        return self.clf.fit(X, y)


@pytest.mark.parametrize(
    "y, params, type_err, msg_err",
    [
        (y_iris, {"estimators": []}, ValueError, "Invalid 'estimators' attribute,"),
        (
            y_iris,
            {
                "estimators": [
                    ("lr", LogisticRegression()),
                    ("svm", SVC(max_iter=50_000)),
                ],
                "stack_method": "predict_proba",
            },
            ValueError,
            "does not implement the method predict_proba",
        ),
        (
            y_iris,
            {
                "estimators": [
                    ("lr", LogisticRegression()),
                    ("cor", NoWeightClassifier()),
                ]
            },
            TypeError,
            "does not support sample weight",
        ),
        (
            y_iris,
            {
                "estimators": [
                    ("lr", LogisticRegression()),
                    ("cor", LinearSVC(max_iter=50_000)),
                ],
                "final_estimator": NoWeightClassifier(),
            },
            TypeError,
            "does not support sample weight",
        ),
    ],
)
def test_stacking_classifier_error(y, params, type_err, msg_err):
    with pytest.raises(type_err, match=msg_err):
        clf = StackingClassifier(**params, cv=3)
        clf.fit(scale(X_iris), y, sample_weight=np.ones(X_iris.shape[0]))


@pytest.mark.parametrize(
    "y, params, type_err, msg_err",
    [
        (y_diabetes, {"estimators": []}, ValueError, "Invalid 'estimators' attribute,"),
        (
            y_diabetes,
            {"estimators": [("lr", LinearRegression()), ("cor", NoWeightRegressor())]},
            TypeError,
            "does not support sample weight",
        ),
        (
            y_diabetes,
            {
                "estimators": [("lr", LinearRegression()), ("cor", LinearSVR())],
                "final_estimator": NoWeightRegressor(),
            },
            TypeError,
            "does not support sample weight",
        ),
    ],
)
def test_stacking_regressor_error(y, params, type_err, msg_err):
    with pytest.raises(type_err, match=msg_err):
        reg = StackingRegressor(**params, cv=3)
        reg.fit(scale(X_diabetes), y, sample_weight=np.ones(X_diabetes.shape[0]))


@pytest.mark.parametrize(
    "estimator, X, y",
    [
        (
            StackingClassifier(
                estimators=[
                    ("lr", LogisticRegression(random_state=0)),
                    ("svm", LinearSVC(random_state=0)),
                ]
            ),
            X_iris[:100],
            y_iris[:100],
        ),  # keep only classes 0 and 1
        (
            StackingRegressor(
                estimators=[
                    ("lr", LinearRegression()),
                    ("svm", LinearSVR(random_state=0)),
                ]
            ),
            X_diabetes,
            y_diabetes,
        ),
    ],
    ids=["StackingClassifier", "StackingRegressor"],
)
def test_stacking_randomness(estimator, X, y):
    # checking that fixing the random state of the CV will lead to the same
    # results
    estimator_full = clone(estimator)
    estimator_full.set_params(
        cv=KFold(shuffle=True, random_state=np.random.RandomState(0))
    )

    estimator_drop = clone(estimator)
    estimator_drop.set_params(lr="drop")
    estimator_drop.set_params(
        cv=KFold(shuffle=True, random_state=np.random.RandomState(0))
    )

    assert_allclose(
        estimator_full.fit(X, y).transform(X)[:, 1:],
        estimator_drop.fit(X, y).transform(X),
    )


def test_stacking_classifier_stratify_default():
    # check that we stratify the classes for the default CV
    clf = StackingClassifier(
        estimators=[
            ("lr", LogisticRegression(max_iter=10_000)),
            ("svm", LinearSVC(max_iter=10_000)),
        ]
    )
    # since iris is not shuffled, a simple k-fold would not contain the
    # 3 classes during training
    clf.fit(X_iris, y_iris)


@pytest.mark.parametrize(
    "stacker, X, y",
    [
        (
            StackingClassifier(
                estimators=[
                    ("lr", LogisticRegression()),
                    ("svm", LinearSVC(random_state=42)),
                ],
                final_estimator=LogisticRegression(),
                cv=KFold(shuffle=True, random_state=42),
            ),
            *load_breast_cancer(return_X_y=True),
        ),
        (
            StackingRegressor(
                estimators=[
                    ("lr", LinearRegression()),
                    ("svm", LinearSVR(random_state=42)),
                ],
                final_estimator=LinearRegression(),
                cv=KFold(shuffle=True, random_state=42),
            ),
            X_diabetes,
            y_diabetes,
        ),
    ],
    ids=["StackingClassifier", "StackingRegressor"],
)
def test_stacking_with_sample_weight(stacker, X, y):
    # check that sample weights has an influence on the fitting
    # note: ConvergenceWarning are catch since we are not worrying about the
    # convergence here
    n_half_samples = len(y) // 2
    total_sample_weight = np.array(
        [0.1] * n_half_samples + [0.9] * (len(y) - n_half_samples)
    )
    X_train, X_test, y_train, _, sample_weight_train, _ = train_test_split(
        X, y, total_sample_weight, random_state=42
    )

    with ignore_warnings(category=ConvergenceWarning):
        stacker.fit(X_train, y_train)
    y_pred_no_weight = stacker.predict(X_test)

    with ignore_warnings(category=ConvergenceWarning):
        stacker.fit(X_train, y_train, sample_weight=np.ones(y_train.shape))
    y_pred_unit_weight = stacker.predict(X_test)

    assert_allclose(y_pred_no_weight, y_pred_unit_weight)

    with ignore_warnings(category=ConvergenceWarning):
        stacker.fit(X_train, y_train, sample_weight=sample_weight_train)
    y_pred_biased = stacker.predict(X_test)

    assert np.abs(y_pred_no_weight - y_pred_biased).sum() > 0


def test_stacking_classifier_sample_weight_fit_param():
    # check sample_weight is passed to all invocations of fit
    stacker = StackingClassifier(
        estimators=[("lr", CheckingClassifier(expected_sample_weight=True))],
        final_estimator=CheckingClassifier(expected_sample_weight=True),
    )
    stacker.fit(X_iris, y_iris, sample_weight=np.ones(X_iris.shape[0]))


@pytest.mark.filterwarnings("ignore::sklearn.exceptions.ConvergenceWarning")
@pytest.mark.parametrize(
    "stacker, X, y",
    [
        (
            StackingClassifier(
                estimators=[
                    ("lr", LogisticRegression()),
                    ("svm", LinearSVC(random_state=42)),
                ],
                final_estimator=LogisticRegression(),
            ),
            *load_breast_cancer(return_X_y=True),
        ),
        (
            StackingRegressor(
                estimators=[
                    ("lr", LinearRegression()),
                    ("svm", LinearSVR(random_state=42)),
                ],
                final_estimator=LinearRegression(),
            ),
            X_diabetes,
            y_diabetes,
        ),
    ],
    ids=["StackingClassifier", "StackingRegressor"],
)
def test_stacking_cv_influence(stacker, X, y):
    # check that the stacking affects the fit of the final estimator but not
    # the fit of the base estimators
    # note: ConvergenceWarning are catch since we are not worrying about the
    # convergence here
    stacker_cv_3 = clone(stacker)
    stacker_cv_5 = clone(stacker)

    stacker_cv_3.set_params(cv=3)
    stacker_cv_5.set_params(cv=5)

    stacker_cv_3.fit(X, y)
    stacker_cv_5.fit(X, y)

    # the base estimators should be identical
    for est_cv_3, est_cv_5 in zip(stacker_cv_3.estimators_, stacker_cv_5.estimators_):
        assert_allclose(est_cv_3.coef_, est_cv_5.coef_)

    # the final estimator should be different
    with pytest.raises(AssertionError, match="Not equal"):
        assert_allclose(
            stacker_cv_3.final_estimator_.coef_, stacker_cv_5.final_estimator_.coef_
        )


@pytest.mark.parametrize(
    "Stacker, Estimator, stack_method, final_estimator, X, y",
    [
        (
            StackingClassifier,
            DummyClassifier,
            "predict_proba",
            LogisticRegression(random_state=42),
            X_iris,
            y_iris,
        ),
        (
            StackingRegressor,
            DummyRegressor,
            "predict",
            LinearRegression(),
            X_diabetes,
            y_diabetes,
        ),
    ],
)
def test_stacking_prefit(Stacker, Estimator, stack_method, final_estimator, X, y):
    """Check the behaviour of stacking when `cv='prefit'`"""
    X_train1, X_train2, y_train1, y_train2 = train_test_split(
        X, y, random_state=42, test_size=0.5
    )
    estimators = [
        ("d0", Estimator().fit(X_train1, y_train1)),
        ("d1", Estimator().fit(X_train1, y_train1)),
    ]

    # mock out fit and stack_method to be asserted later
    for _, estimator in estimators:
        estimator.fit = Mock()
        stack_func = getattr(estimator, stack_method)
        setattr(estimator, stack_method, Mock(side_effect=stack_func))

    stacker = Stacker(
        estimators=estimators, cv="prefit", final_estimator=final_estimator
    )
    stacker.fit(X_train2, y_train2)

    assert stacker.estimators_ == [estimator for _, estimator in estimators]
    # fit was not called again
    assert all(estimator.fit.call_count == 0 for estimator in stacker.estimators_)

    # stack method is called with the proper inputs
    for estimator in stacker.estimators_:
        stack_func_mock = getattr(estimator, stack_method)
        stack_func_mock.assert_called_with(X_train2)


@pytest.mark.parametrize(
    "stacker, X, y",
    [
        (
            StackingClassifier(
                estimators=[("lr", LogisticRegression()), ("svm", SVC())],
                cv="prefit",
            ),
            X_iris,
            y_iris,
        ),
        (
            StackingRegressor(
                estimators=[
                    ("lr", LinearRegression()),
                    ("svm", LinearSVR()),
                ],
                cv="prefit",
            ),
            X_diabetes,
            y_diabetes,
        ),
    ],
)
def test_stacking_prefit_error(stacker, X, y):
    # check that NotFittedError is raised
    # if base estimators are not fitted when cv="prefit"
    with pytest.raises(NotFittedError):
        stacker.fit(X, y)


@pytest.mark.parametrize(
    "make_dataset, Stacking, Estimator",
    [
        (make_classification, StackingClassifier, LogisticRegression),
        (make_regression, StackingRegressor, LinearRegression),
    ],
)
def test_stacking_without_n_features_in(make_dataset, Stacking, Estimator):
    # Stacking supports estimators without `n_features_in_`. Regression test
    # for #17353

    class MyEstimator(Estimator):
        """Estimator without n_features_in_"""

        def fit(self, X, y):
            super().fit(X, y)
            del self.n_features_in_

    X, y = make_dataset(random_state=0, n_samples=100)
    stacker = Stacking(estimators=[("lr", MyEstimator())])

    msg = f"{Stacking.__name__} object has no attribute n_features_in_"
    with pytest.raises(AttributeError, match=msg):
        stacker.n_features_in_

    # Does not raise
    stacker.fit(X, y)

    msg = "'MyEstimator' object has no attribute 'n_features_in_'"
    with pytest.raises(AttributeError, match=msg):
        stacker.n_features_in_


@pytest.mark.parametrize("stack_method", ["auto", "predict", "predict_proba"])
@pytest.mark.parametrize("passthrough", [False, True])
@pytest.mark.parametrize(
    "X, y, is_multilabel",
    [(X_multilabel_r, y_multilabel_r, True), (X_binary, y_binary, False)],
    ids=["multilabel", "binary"],
)
def test_stacking_classifier_multilabel(stack_method, passthrough, X, y, is_multilabel):
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, stratify=y, random_state=42
    )

    estimators = [
        ("mlp", MLPClassifier()),
        ("knnc", RandomForestClassifier()),
        ("dcs", DummyClassifier(strategy="stratified", random_state=42)),
        ("dcu", DummyClassifier(strategy="uniform", random_state=42)),
    ]
    final_estimator = KNeighborsClassifier()

    clf = StackingClassifier(
        estimators=estimators,
        final_estimator=final_estimator,
        passthrough=passthrough,
        stack_method=stack_method,
    )
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    assert y_pred.shape == y_test.shape

    y_proba = clf.predict_proba(X_test)
    if is_multilabel:
        assert y_proba.shape == y_test.shape
    else:
        assert y_proba.shape == (y_test.shape[0], 2)

    X_trans = clf.transform(X_test)
    expected_n_features = len(estimators)
    if is_multilabel:
        expected_n_features *= y_test.shape[1]
    if passthrough:
        expected_n_features += X_test.shape[1]
    assert X_trans.shape[1] == expected_n_features

    if is_multilabel:
        assert len(clf.classes_) == y_test.shape[1]
    else:
        assert len(clf.classes_) == 2

    if passthrough:
        X_test_last_cols = X_trans[:, -X_test.shape[1] :]
        assert_array_equal(X_test, X_test_last_cols)


@pytest.mark.parametrize(
    "stacker, feature_names, X, y, expected_names",
    [
        (
            StackingClassifier(
                estimators=[
                    ("lr", LogisticRegression(random_state=0)),
                    ("svm", LinearSVC(random_state=0)),
                ]
            ),
            iris.feature_names,
            X_iris,
            y_iris,
            [
                "stackingclassifier_lr0",
                "stackingclassifier_lr1",
                "stackingclassifier_lr2",
                "stackingclassifier_svm0",
                "stackingclassifier_svm1",
                "stackingclassifier_svm2",
            ],
        ),
        (
            StackingClassifier(
                estimators=[
                    ("lr", LogisticRegression(random_state=0)),
                    ("other", "drop"),
                    ("svm", LinearSVC(random_state=0)),
                ]
            ),
            iris.feature_names,
            X_iris[:100],
            y_iris[:100],  # keep only classes 0 and 1
            [
                "stackingclassifier_lr",
                "stackingclassifier_svm",
            ],
        ),
        (
            StackingRegressor(
                estimators=[
                    ("lr", LinearRegression()),
                    ("svm", LinearSVR(random_state=0)),
                ]
            ),
            diabetes.feature_names,
            X_diabetes,
            y_diabetes,
            [
                "stackingregressor_lr",
                "stackingregressor_svm",
            ],
        ),
    ],
    ids=[
        "StackingClassifier_multiclass",
        "StackingClassifier_binary",
        "StackingRegressor",
    ],
)
@pytest.mark.parametrize("passthrough", [True, False])
def test_get_feature_names_out(
    stacker, feature_names, X, y, expected_names, passthrough
):
    """Check get_feature_names_out works for stacking."""

    stacker.set_params(passthrough=passthrough)
    stacker.fit(scale(X), y)

    if passthrough:
        expected_names = np.concatenate((expected_names, feature_names))

    names_out = stacker.get_feature_names_out(feature_names)
    assert_array_equal(names_out, expected_names)

# Author: Maria Telenczuk <https://github.com/maikia>
#
# License: BSD 3 clause

import pytest

import numpy as np

from sklearn.base import is_classifier
from sklearn.datasets import make_classification, make_regression
from sklearn.preprocessing import MinMaxScaler
from sklearn.utils import check_random_state
from sklearn.utils._testing import set_random_state

from sklearn.linear_model import (
    ARDRegression,
    BayesianRidge,
    ElasticNet,
    ElasticNetCV,
    GammaRegressor,
    HuberRegressor,
    Lars,
    LarsCV,
    Lasso,
    LassoCV,
    LassoLars,
    LassoLarsCV,
    LassoLarsIC,
    LinearRegression,
    LogisticRegression,
    LogisticRegressionCV,
    OrthogonalMatchingPursuit,
    OrthogonalMatchingPursuitCV,
    PassiveAggressiveClassifier,
    PassiveAggressiveRegressor,
    Perceptron,
    PoissonRegressor,
    Ridge,
    RidgeCV,
    RidgeClassifier,
    RidgeClassifierCV,
    SGDClassifier,
    SGDRegressor,
    TheilSenRegressor,
    TweedieRegressor,
)
from sklearn.svm import LinearSVC, LinearSVR


@pytest.mark.parametrize(
    "normalize, n_warnings, warning_category",
    [(True, 1, FutureWarning), (False, 1, FutureWarning), ("deprecated", 0, None)],
)
@pytest.mark.parametrize(
    "estimator",
    [
        LinearRegression,
        Ridge,
        RidgeCV,
        RidgeClassifier,
        RidgeClassifierCV,
        BayesianRidge,
        ARDRegression,
    ],
)
# FIXME remove test in 1.2
def test_linear_model_normalize_deprecation_message(
    estimator, normalize, n_warnings, warning_category
):
    # check that we issue a FutureWarning when normalize was set in
    # linear model
    rng = check_random_state(0)
    n_samples = 200
    n_features = 2
    X = rng.randn(n_samples, n_features)
    X[X < 0.1] = 0.0
    y = rng.rand(n_samples)
    if is_classifier(estimator):
        y = np.sign(y)

    model = estimator(normalize=normalize)
    with pytest.warns(warning_category) as record:
        model.fit(X, y)
    # Filter record in case other unrelated warnings are raised
    unwanted = [r for r in record if r.category != warning_category]
    if len(unwanted):
        msg = "unexpected warnings:\n"
        for w in unwanted:
            msg += str(w)
            msg += "\n"
        raise AssertionError(msg)
    wanted = [r for r in record if r.category == warning_category]
    if warning_category is not None:
        assert "'normalize' was deprecated" in str(wanted[0].message)
    assert len(wanted) == n_warnings


@pytest.mark.filterwarnings("ignore:The default of 'normalize'")
@pytest.mark.filterwarnings("ignore:lbfgs failed to converge")
@pytest.mark.parametrize(
    "Regressor",
    [
        ARDRegression,
        BayesianRidge,
        ElasticNet,
        ElasticNetCV,
        GammaRegressor,
        HuberRegressor,
        Lars,
        LarsCV,
        Lasso,
        LassoCV,
        LassoLars,
        LassoLarsCV,
        LassoLarsIC,
        LinearSVR,
        LinearRegression,
        OrthogonalMatchingPursuit,
        OrthogonalMatchingPursuitCV,
        PassiveAggressiveRegressor,
        PoissonRegressor,
        Ridge,
        RidgeCV,
        SGDRegressor,
        TheilSenRegressor,
        TweedieRegressor,
    ],
)
@pytest.mark.parametrize("ndim", [1, 2])
def test_linear_model_regressor_coef_shape(Regressor, ndim):
    """Check the consistency of linear models `coef` shape."""
    X, y = make_regression(random_state=0)
    y = MinMaxScaler().fit_transform(y.reshape(-1, 1))[:, 0] + 1
    y = y[:, np.newaxis] if ndim == 2 else y

    regressor = Regressor()
    set_random_state(regressor)
    regressor.fit(X, y)
    assert regressor.coef_.shape == (X.shape[1],)


@pytest.mark.parametrize(
    "Classifier",
    [
        LinearSVC,
        LogisticRegression,
        LogisticRegressionCV,
        PassiveAggressiveClassifier,
        Perceptron,
        RidgeClassifier,
        RidgeClassifierCV,
        SGDClassifier,
    ],
)
@pytest.mark.parametrize("n_classes", [2, 3])
def test_linear_model_classifier_coef_shape(Classifier, n_classes):
    X, y = make_classification(n_informative=10, n_classes=n_classes, random_state=0)
    n_features = X.shape[1]

    classifier = Classifier()
    set_random_state(classifier)
    classifier.fit(X, y)
    expected_shape = (1, n_features) if n_classes == 2 else (n_classes, n_features)
    assert classifier.coef_.shape == expected_shape

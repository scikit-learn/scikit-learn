# Author: Maria Telenczuk <https://github.com/maikia>
#
# License: BSD 3 clause

import inspect
import sys
import warnings

import numpy as np
import pytest

from sklearn.base import is_classifier
from sklearn.datasets import make_low_rank_matrix
from sklearn.linear_model import RidgeClassifier
from sklearn.linear_model import RidgeClassifierCV
from sklearn.linear_model import (
    ARDRegression,
    BayesianRidge,
    ElasticNet,
    ElasticNetCV,
    Lars,
    LarsCV,
    Lasso,
    LassoCV,
    LassoLarsCV,
    LassoLarsIC,
    LinearRegression,
    LogisticRegression,
    LogisticRegressionCV,
    MultiTaskElasticNet,
    MultiTaskElasticNetCV,
    MultiTaskLasso,
    MultiTaskLassoCV,
    OrthogonalMatchingPursuit,
    OrthogonalMatchingPursuitCV,
    PoissonRegressor,
    Ridge,
    RidgeCV,
    SGDRegressor,
    TweedieRegressor,
)

from sklearn.utils.fixes import np_version, parse_version
from sklearn.utils import check_random_state


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
@pytest.mark.xfail(
    sys.platform == "darwin" and np_version < parse_version("1.22"),
    reason="https://github.com/scikit-learn/scikit-learn/issues/21395",
)
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
    if warning_category is None:
        with warnings.catch_warnings():
            warnings.simplefilter("error", FutureWarning)
            model.fit(X, y)
        return

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
    assert "'normalize' was deprecated" in str(wanted[0].message)
    assert len(wanted) == n_warnings


# Note: GammaRegressor() and TweedieRegressor(power != 1) have a non-canonical link.
@pytest.mark.parametrize(
    "model",
    [
        ARDRegression(),
        BayesianRidge(),
        ElasticNet(),
        ElasticNetCV(),
        Lars(),
        LarsCV(),
        Lasso(),
        LassoCV(),
        LassoLarsCV(),
        LassoLarsIC(),
        LinearRegression(),
        # TODO: Why does this fail badly with sample_weights?
        pytest.param(
            LogisticRegression(
                penalty="elasticnet", solver="saga", l1_ratio=0.5, tol=1e-15
            ),
            marks=pytest.mark.xfail(reason="Unknown bug"),
        ),
        LogisticRegressionCV(),
        MultiTaskElasticNet(),
        MultiTaskElasticNetCV(),
        MultiTaskLasso(),
        MultiTaskLassoCV(),
        OrthogonalMatchingPursuit(normalize=False),
        OrthogonalMatchingPursuitCV(normalize=False),
        PoissonRegressor(),
        Ridge(),
        RidgeCV(),
        pytest.param(
            SGDRegressor(tol=1e-15),
            marks=pytest.mark.xfail(reason="Unsufficient precision."),
        ),
        SGDRegressor(penalty="elasticnet"),
        TweedieRegressor(power=0),  # same as Ridge
    ],
    ids=lambda x: x.__class__.__name__,
)
@pytest.mark.parametrize("with_sample_weight", [False, True])
def test_balance_property(model, with_sample_weight, global_random_seed):
    # Test that sum(y_predicted) == sum(y_observed) on the training set.
    # This must hold for all linear models with deviance of an exponential disperson
    # family as loss and the corresponding canonical link if fit_intercept=True.
    # Examples:
    #     - squared error and identity link (most linear models)
    #     - Poisson deviance with log link
    #     - log loss with logit link
    # This is known as balance property or unconditional calibration/unbiasedness.

    if (
        with_sample_weight
        and "sample_weight" not in inspect.signature(model.fit).parameters.keys()
    ):
        pytest.skip("Estimator does not support sample_weight.")

    rel = 1e-4  # test precision
    # TODO: Investigate why these models achieve less accurate results for the
    #       intercept.
    if isinstance(model, SGDRegressor) or (
        hasattr(model, "solver") and model.solver == "saga"
    ):
        rel = 1e-2

    rng = np.random.RandomState(global_random_seed)
    n_train, n_test, n_features, n_targets = 100, 50, 10, None
    if isinstance(
        model,
        (MultiTaskElasticNet, MultiTaskElasticNetCV, MultiTaskLasso, MultiTaskLassoCV),
    ):
        n_targets = 3
    X = make_low_rank_matrix(
        n_samples=n_train + n_test, n_features=n_features, random_state=rng
    )
    if n_targets:
        coef = (
            rng.uniform(low=-2, high=2, size=(n_features, n_targets))
            / np.max(X, axis=0)[:, None]
        )
    else:
        coef = rng.uniform(low=-2, high=2, size=n_features) / np.max(X, axis=0)

    expectation = np.exp(X @ coef + 0.5)
    y = rng.poisson(lam=expectation) + 1  # strict positive, i.e. y > 0
    if is_classifier(model):
        y = (y > expectation + 1).astype(np.float64)

    if with_sample_weight:
        sw = rng.uniform(low=1, high=10, size=y.shape[0])
    else:
        sw = None

    model.set_params(fit_intercept=True)  # to be sure
    if with_sample_weight:
        model.fit(X, y, sample_weight=sw)
    else:
        model.fit(X, y)

    if is_classifier(model):
        assert np.average(model.predict_proba(X)[:, 1], weights=sw) == pytest.approx(
            np.average(y, weights=sw), rel=rel
        )
    else:
        assert np.average(model.predict(X), weights=sw, axis=0) == pytest.approx(
            np.average(y, weights=sw, axis=0), rel=rel
        )

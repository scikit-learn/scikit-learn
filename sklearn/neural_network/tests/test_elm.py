import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn.base import clone
from sklearn.datasets import (
    fetch_openml,
    load_breast_cancer,
    load_digits,
    make_classification,
    make_regression,
)
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.neural_network import ExtremeLearningClassifier, ExtremeLearningRegressor
from sklearn.neural_network._base import ACTIVATIONS, WEIGHTS
from sklearn.preprocessing import MinMaxScaler, OneHotEncoder, StandardScaler
from sklearn.utils.estimator_checks import parametrize_with_checks

ACTIVATION_NAMES = [name for name in ACTIVATIONS.keys() if name != "exp"]
WEIGHT_NAMES = list(WEIGHTS.keys())


@pytest.mark.parametrize(
    "hidden_layer_sizes", [(10,), (10, 10), (5, 10, 15, 20), (100,)]
)
@pytest.mark.parametrize("n_classes", [2, 5])
@pytest.mark.parametrize("direct_links", [True, False])
@pytest.mark.parametrize("activation", ACTIVATION_NAMES)
@pytest.mark.parametrize("weight_init", WEIGHT_NAMES)
def test_params(hidden_layer_sizes, n_classes, activation, weight_init, direct_links):
    # Test shapes of internal params
    N, d = 60, 10
    X, y = make_classification(
        n_samples=N, n_features=d, n_classes=n_classes, n_informative=8, random_state=42
    )

    model = ExtremeLearningClassifier(
        hidden_layer_sizes, activation, weight_init, direct_links, 0
    )

    model.fit(X, y)

    assert len(model.W_) == len(hidden_layer_sizes)
    assert model.W_[0].T.shape == (d, hidden_layer_sizes[0])

    for layer, w, b, i in zip(
        hidden_layer_sizes[1:],
        model.W_[1:],
        model.b_[1:],
        range(len(model.W_) - 1),
        strict=False,
    ):
        assert w.T.shape == (hidden_layer_sizes[i], layer)
        assert b.shape == (layer,)

    if direct_links:
        assert model.coef_.shape == (
            sum(hidden_layer_sizes) + d,
            n_classes,
        )
    else:
        assert model.coef_.shape == (sum(hidden_layer_sizes), n_classes)

    pred = model.predict(X[:10])
    assert np.isin(pred, np.arange(n_classes)).all()

    P = model.predict_proba(X[:10])
    assert_allclose(P.sum(axis=1), 1.0, atol=1e-6)
    assert (P >= 0).all() and (P <= 1).all()
    np.testing.assert_array_equal(pred, model.classes_[np.argmax(P, axis=1)])


@pytest.mark.parametrize("weight_init", WEIGHT_NAMES)
@pytest.mark.parametrize(
    "hidden_layer_size", [(10,), (2, 3, 2, 1), (5, 10, 15, 20, 15, 10), (100,)]
)
def test_multilayer_math(weight_init, hidden_layer_size):
    # Validates construction of design matrix algebraically
    #
    # H_1 = X @ W_1.T + b_1
    # H_2 = H_1 @ W_2.T + b_2
    # ...
    # H_N = H_(N−1) @ W_N.T + b_N
    #
    # When activation is identity, we can collapse these layers
    # into a single linear equation:
    #
    # H_N = X @ (W_1.T @ ... @ W_N.T) + c_N
    #
    # where the accumulated bias (c_N) is:
    #
    # c_N =
    # b_1 @ W_2.T @ W_3.T @ ... @ W_N.T +
    # b_2 @ W_3.T @ W_4.T @ ... @ W_N.T +
    # ... +
    # b_(N-1) @ W_N.T +
    # b_N

    N, d = 60, 10
    X, y = make_classification(
        n_samples=N, n_features=d, n_classes=3, n_informative=8, random_state=42
    )

    model = ExtremeLearningClassifier(
        hidden_layer_sizes=hidden_layer_size,
        activation="identity",
        weight_init=weight_init,
        direct_links=False,
        random_state=0,
    )

    model.fit(X, y)

    enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    Y = enc.fit_transform(y.reshape(-1, 1))

    # collapsing weights and biases for representation as linear operation
    weights = [w.T for w in model.W_]
    Ts, cs = [], []
    T = np.eye(X.shape[1])
    c = np.zeros((X.shape[1],))

    for w, b in zip(weights, model.b_, strict=False):
        T = T @ w
        c = c @ w + b
        Ts.append(T)
        cs.append(c)

    # design matrix with ALL layers concatenated
    expected_phi = np.hstack([X @ T_l + c_l for T_l, c_l in zip(Ts, cs, strict=False)])

    expected_beta = np.linalg.pinv(expected_phi) @ Y

    assert_allclose(model.coef_, expected_beta)


@pytest.mark.parametrize("n_class", [3, 2])
@pytest.mark.parametrize("activation", ACTIVATION_NAMES)
@pytest.mark.parametrize("direct_links", [True, False])
def test_lls_classification(n_class, activation, direct_links):
    # Test linear least squares on classification.
    # It should achieve a score higher than 0.93 for the binary and multi-class
    # versions of the digits dataset.
    X, y = load_digits(n_class=n_class, return_X_y=True)

    X = MinMaxScaler().fit_transform(X[:200])
    y = y[:200]

    X_train = X[:150]
    y_train = y[:150]
    X_test = X[150:]
    expected_shape_dtype = (X_test.shape[0], y_train.dtype.kind)

    elm = ExtremeLearningClassifier(
        random_state=1,
        activation=activation,
        direct_links=direct_links,
    )
    elm.fit(X_train, y_train)
    y_predict = elm.predict(X_test)
    assert elm.score(X_train, y_train) > 0.93
    assert (y_predict.shape[0], y_predict.dtype.kind) == expected_shape_dtype


@pytest.mark.parametrize("activation", ACTIVATION_NAMES)
@pytest.mark.parametrize("direct_links", [True, False])
def test_lls_regression(activation, direct_links):
    # Test linear least squares on the regression dataset.

    X, y = make_regression(n_samples=200, n_features=10, random_state=1)
    elm = ExtremeLearningRegressor(
        random_state=1,
        activation=activation,
        direct_links=direct_links,
    )
    elm.fit(X, y)
    assert elm.score(X, y) > 0.73


@pytest.mark.parametrize("hidden_layer_sizes", [(10,), (5, 5)])
@pytest.mark.parametrize("direct_links", [True, False])
@pytest.mark.parametrize("n_classes", [2, 5])
@pytest.mark.parametrize("activation", ACTIVATION_NAMES)
@pytest.mark.parametrize("weight_init", WEIGHT_NAMES[1:])
@pytest.mark.parametrize("alpha", [None, 0.1, 0.5, 1])
def test_partial_fit_classifier(
    hidden_layer_sizes, direct_links, n_classes, activation, weight_init, alpha
):
    # Test coefficient equivalence between partial_fit() and fit() as long
    # as D.T @ D is well-conditioned

    # partial_fit is equivalent to accumulating normal equations
    X, y = make_classification(
        n_samples=600,
        n_features=20,
        n_informative=10,
        n_redundant=2,
        n_classes=n_classes,
        random_state=0,
    )

    ff_model = ExtremeLearningClassifier(
        hidden_layer_sizes=hidden_layer_sizes,
        activation=activation,
        weight_init=weight_init,
        direct_links=direct_links,
        random_state=0,
        ridge_alpha=alpha,
    )
    pf_model = clone(ff_model)

    ff_model.fit(X, y)

    classes = np.unique(y)
    batch = 25
    for start in range(0, len(X), batch):
        end = min(start + batch, len(X))
        Xb = X[start:end]
        yb = y[start:end]
        if start == 0:
            pf_model.partial_fit(Xb, yb, classes=classes)
        else:
            pf_model.partial_fit(Xb, yb)

    assert_allclose(ff_model.coef_, pf_model.coef_, rtol=1e-5, atol=1e-5)


@pytest.mark.parametrize("hidden_layer_sizes", [(10,), (5, 5)])
@pytest.mark.parametrize("direct_links", [True, False])
@pytest.mark.parametrize("activation", ACTIVATION_NAMES)
@pytest.mark.parametrize("weight_init", WEIGHT_NAMES[1:])
@pytest.mark.parametrize("alpha", [None, 0.1, 0.5, 1])
def test_partial_fit_regressor(
    hidden_layer_sizes, direct_links, activation, weight_init, alpha
):
    # Test coefficient equivalence between partial_fit() and fit() as long
    # as D.T@D is well-conditioned

    # partial_fit is equivalent to accumulating normal equations
    X, y = make_regression(
        n_samples=600,
        n_features=20,
        n_informative=10,
        random_state=0,
    )

    ff_model = ExtremeLearningRegressor(
        hidden_layer_sizes=hidden_layer_sizes,
        activation=activation,
        weight_init=weight_init,
        direct_links=direct_links,
        random_state=0,
        ridge_alpha=alpha,
    )
    pf_model = clone(ff_model)

    ff_model.fit(X, y)

    batch = 25
    for start in range(0, len(X), batch):
        end = min(start + batch, len(X))
        Xb = X[start:end]
        yb = y[start:end]
        pf_model.partial_fit(Xb, yb)

    assert_allclose(ff_model.coef_, pf_model.coef_, rtol=1e-5, atol=4e-5)


@pytest.mark.parametrize(
    "ridge_alpha, expected",
    [
        (0.1, 0.8327257206362838),
        # NOTE: for Moore-Penrose, a large singular value
        # cutoff (rcond) is required to achieve reasonable R2 with
        # the Boston Housing dataset
        # Without rtol, R2 ~= -59.615
        (None, 0.764164981634744),
    ],
)
def test_regression_boston(ridge_alpha, expected):
    # real-world data test with multi-layer RVFL
    X, y = fetch_openml(
        name="boston", version=1, return_X_y=True, as_frame=False, parser="liac-arff"
    )
    y = y.astype(float)
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.2,
        random_state=42,
        shuffle=True,
    )

    scaler = StandardScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    model = ExtremeLearningRegressor(
        hidden_layer_sizes=[800] * 10,
        activation="logistic",
        weight_init="uniform",
        direct_links=True,
        random_state=0,
        ridge_alpha=ridge_alpha,
        rtol=1e-3,  # has no effect for `Ridge`
    )
    model.fit(X_train, y_train)
    # RandomForestRegressor() with default params scores
    # 0.8733907 here; multi-layer RVFL with above params is a bit
    # worse, but certainly better than random chance:
    actual = model.score(X_test, y_test)
    assert_allclose(actual, expected)


@pytest.mark.parametrize(
    "ridge_alpha, expected",
    [
        (0.1, 0.9473684210526315),
        # NOTE: for Moore-Penrose, a large singular value
        # cutoff (rcond) is required to achieve reasonable accuracy
        # with the Wisconsin breast cancer dataset
        # Without rtol accuracy ~= 0.7105
        (None, 0.9473684210526315),
    ],
)
def test_classification_breast_cancer(ridge_alpha, expected):
    # real-world data test with multi-layer RVFL
    X, y = load_breast_cancer(return_X_y=True)

    X = MinMaxScaler().fit_transform(X)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, shuffle=True
    )

    scaler = StandardScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    model = ExtremeLearningClassifier(
        hidden_layer_sizes=[800] * 10,
        activation="logistic",
        weight_init="uniform",
        direct_links=True,
        random_state=0,
        ridge_alpha=ridge_alpha,
        rtol=1e-3,  # has no effect for `Ridge`
    )
    model.fit(X_train, y_train)
    # RandomForestRegressor() with default params scores 0.958 here
    # RVFL with above params scores comparatively:
    actual = model.score(X_test, y_test)
    assert_allclose(actual, expected)


@pytest.mark.parametrize(
    "ridge_alpha, expected",
    [
        (0.1, 0.8327240971436314),
        # NOTE: for Moore-Penrose, a large singular value
        # cutoff (rcond) is required to achieve reasonable R2
        # with the Boston Housing dataset
        # Without rtol R2 ~= -59.614
        (None, 0.8083442423720848),
    ],
)
def test_rtol_partial_fit(ridge_alpha, expected):
    # this test takes ~41.70s... should we simplify?
    X, y = fetch_openml(
        name="boston", version=1, return_X_y=True, as_frame=False, parser="liac-arff"
    )
    y = y.astype(float)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, shuffle=True
    )
    scaler = StandardScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    model = ExtremeLearningRegressor(
        hidden_layer_sizes=[800]
        * 2,  # partial_fit is slow so smaller network for speed
        activation="logistic",
        weight_init="uniform",
        direct_links=True,
        random_state=0,
        ridge_alpha=ridge_alpha,
        rtol=1e-6,
    )

    batch = 50
    for start in range(0, X_train.shape[0], batch):
        end = min(start + batch, X_train.shape[0])
        model.partial_fit(X_train[start:end], y_train[start:end])

    actual = model.score(X_test, y_test)
    assert_allclose(actual, expected)


def test_partial_fit_classes_error():
    # Test partial_fit error handling.
    X, y = make_classification(n_samples=50, n_features=10, n_classes=2, random_state=0)
    clf = ExtremeLearningClassifier(random_state=0)

    with pytest.raises(TypeError, match="Classes must not be None"):
        # classes parameter is required for first partial_fit call
        clf.partial_fit(X[:25], y[:25])

    clf.partial_fit(X[:25], y[:25], classes=np.array([0, 1]))
    y_bad = y[25:].copy()
    y_bad[0] = 2
    with pytest.raises(ValueError, match="Expected only labels in classes_"):
        # Raised when unseen classes are passed after initial partial_fit call
        clf.partial_fit(X[25:], y_bad)


def test_batch_order_invariance():
    # Order-invariance test for partial_fit
    X, y = make_classification(
        n_samples=400,
        n_features=20,
        n_informative=10,
        n_redundant=5,
        n_classes=2,
        random_state=0,
    )

    ff_model = ExtremeLearningClassifier(random_state=0, ridge_alpha=0.1)
    pf_model = clone(ff_model)

    ff_model.fit(X, y)

    classes = np.unique(y)
    batch = 10
    rng = np.random.default_rng(0)
    indices = rng.permutation(np.arange(0, len(X), batch))
    for start in indices:
        end = min(start + batch, len(X))
        Xb = X[start:end]
        yb = y[start:end]
        if start == indices[0]:
            pf_model.partial_fit(Xb, yb, classes=classes)
        else:
            pf_model.partial_fit(Xb, yb)

    assert_allclose(ff_model.coef_, pf_model.coef_, rtol=4e-8, atol=2e-10)


def test_batch_partition_invariance():
    # Partition-invariance test for partial_fit
    X, y = make_classification(
        n_samples=400,
        n_features=20,
        n_informative=10,
        n_redundant=5,
        n_classes=2,
        random_state=0,
    )
    pf1 = ExtremeLearningClassifier(random_state=0, ridge_alpha=0.1)
    pf2 = clone(pf1)

    classes = np.unique(y)
    for i in range(0, len(X), 10):
        pf1.partial_fit(X[i : i + 10], y[i : i + 10], classes=classes)

    cuts = [17, 61, 140]
    starts = [0] + cuts
    ends = cuts + [len(X)]
    for s, e in zip(starts, ends):
        pf2.partial_fit(X[s:e], y[s:e], classes=classes)

    assert_allclose(pf1.coef_, pf2.coef_, rtol=1e-7, atol=1e-9)


@pytest.mark.xfail(
    reason="partial_fit diverges from fit for ill-conditioned design matrices",
    strict=True,
)
def test_partial_fit_ill_conditioned():
    # For direct_links=True and certain activations and weight combinations,
    # the design matrix becomes rank-deficient and the exact
    # solve can diverge from fit()

    X, y = make_classification(
        n_samples=100,
        n_features=10,
        n_informative=4,
        n_classes=5,
        random_state=0,
    )

    ff_model = ExtremeLearningClassifier(
        hidden_layer_sizes=(50, 50),
        activation="softmax",
        weight_init="uniform",
        direct_links=True,
        random_state=0,
        ridge_alpha=None,
    )
    pf_model = clone(ff_model)

    ff_model.fit(X, y)

    classes = np.unique(y)
    batch = 10
    for start in range(0, len(X), batch):
        end = min(start + batch, len(X))
        if start == 0:
            pf_model.partial_fit(X[start:end], y[start:end], classes=classes)
        else:
            pf_model.partial_fit(X[start:end], y[start:end])

    # partial_fit() is expected to diverge from fit() given
    # these params
    assert_allclose(pf_model.coef_, ff_model.coef_, rtol=1e-3, atol=1e-3)


@pytest.mark.parametrize(
    "hidden_layer_sizes, activation, weight_init, exp_auc",
    [
        # when direct links are absent (ELM), we expect the
        # ROC AUC to increase with multi-layer network complexity
        # up to a reasonable degree, when the width of the layers is
        # quite small
        ((2,), "relu", "uniform", 0.5608145371438336),
        ((2, 2), "relu", "uniform", 0.5686688880345168),
        # start hitting diminishing returns here:
        ((2, 2, 2, 2), "relu", "uniform", 0.5686688880345168),
        ((2, 2, 2, 2, 2, 2, 2), "relu", "uniform", 0.5686688880345168),
        # effectively no improvement here:
        ((2, 2, 2, 2, 2, 2, 2, 2, 2), "relu", "uniform", 0.5686688880345168),
    ],
)
def test_multilayer_progression(weight_init, hidden_layer_sizes, activation, exp_auc):
    # Regression test for multilayer behavior
    X, y = make_classification(
        n_samples=400,
        n_features=100,
        n_classes=5,
        n_informative=26,
        random_state=42,
        class_sep=0.5,
    )
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=0
    )
    model = ExtremeLearningClassifier(
        hidden_layer_sizes=hidden_layer_sizes,
        activation=activation,
        weight_init=weight_init,
        direct_links=False,
        random_state=0,
    )
    model.fit(X_train, y_train)
    y_score = model.predict_proba(X_test)
    actual_auc = roc_auc_score(y_test, y_score, multi_class="ovo")

    # This test is quite brittle
    assert_allclose(actual_auc, exp_auc, rtol=1e-2, atol=5e-3)


@parametrize_with_checks([ExtremeLearningClassifier(), ExtremeLearningRegressor()])
def test_sklearn_api_conformance(estimator, check):
    check(estimator)

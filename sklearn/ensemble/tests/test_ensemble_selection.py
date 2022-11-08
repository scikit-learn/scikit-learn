"""Test the ensemble_selection classifier and regressor."""
import numpy as np
from scipy import sparse
import pytest

from sklearn.datasets import load_diabetes
from sklearn.feature_extraction.text import TfidfVectorizer

from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVR
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.neural_network import MLPClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn import preprocessing
from sklearn.preprocessing import scale
from sklearn.metrics import log_loss, mean_squared_error

from sklearn.cluster import KMeans, DBSCAN

from sklearn.ensemble import EnsembleSelection
from sklearn.pipeline import Pipeline


def _get_trained_evaluated_estimators(library_name, seed=1, X=None, y=None):
    # build estimators
    estimators = []
    if library_name == "regressor":
        estimators.append(("est1", GaussianProcessRegressor(random_state=seed)))
        estimators.append(("est2", RandomForestRegressor(random_state=seed)))
        estimators.append(("est3", KNeighborsRegressor()))
        estimators.append(("est4", MLPRegressor(random_state=seed)))
        estimators.append(("est5", LinearSVR(random_state=seed)))
    elif library_name == "classifier":
        estimators.append(("est1", GaussianProcessClassifier(random_state=seed)))
        estimators.append(("est2", SVC(random_state=seed, probability=True)))
        estimators.append(("est3", RandomForestClassifier(random_state=seed)))
        estimators.append(("est5", MLPClassifier(random_state=seed)))
        estimators.append(("est4", KNeighborsClassifier()))
    elif library_name == "short_regressor":
        estimators.append(("est1", MLPRegressor(random_state=seed)))
        estimators.append(("est2", LinearSVR(random_state=seed)))
    elif library_name == "short_classifier":
        estimators.append(("est1", LogisticRegression(C=0.1, random_state=seed)))
        estimators.append(("est2", LogisticRegression(C=1, random_state=seed)))

    # fit estimators if data is given
    if X is not None and y is not None:
        for _, base_model in estimators:
            base_model.fit(X, y)
    return estimators


def _get_simple_data():
    dtype = float
    X = np.array(
        [[0, 1], [0, 2], [0, 3], [0, 4], [1, 1], [1, 3], [1, 5], [1, 6]], dtype=dtype
    )
    y = np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=dtype)
    return X, y


def test_ensemble_selection_realistic_example():
    seed = 1

    diabetes = load_diabetes()
    X, y = diabetes.data, diabetes.target
    X = scale(X)
    y = scale(y)

    X_train, y_train = X[:200], y[:200]  # 200 samples
    X_valid, y_valid = X[200:300], y[200:300]  # 100 samples
    X_test, y_test = X[300:], y[300:]  # 142 samples

    estimators = _get_trained_evaluated_estimators("regressor", seed, X_train, y_train)

    base_score = []
    for _, est in estimators:
        base_score.append(est.score(X_test, y_test))

    # fit ensemble_selector
    es = EnsembleSelection(
        estimators=estimators,
        score_metric=mean_squared_error,
        score_direction="min",
        min_estimators=2,
        max_estimators=10,
        pruning_factor=0.2,
        with_replacement=True,
        verbose=True,
    )
    es.fit(X_valid, y_valid)

    # predict
    score = es.score(X_test, y_test)

    # at least better than the worst one + margin
    assert score < np.max(base_score) + 0.1


def test_ensemble_selection_predict_regressor():
    seed = 1
    X, y = _get_simple_data()

    estimators = _get_trained_evaluated_estimators("regressor", seed, X, y)

    # fit ensemble_selector
    es = EnsembleSelection(
        estimators=estimators,
        score_metric=mean_squared_error,
        score_direction="min",
        is_base_estimator_proba=False,
        min_estimators=2,
        max_estimators=10,
        with_replacement=True,
        verbose=True,
    )
    es.fit(X, y)

    # predict
    y_pred = es.predict(X)
    loss1 = mean_squared_error(y, y_pred)
    assert loss1 < 5

    # predict_proba is implemented like predict for regressors
    y_pred = es.predict(X)
    loss2 = mean_squared_error(y, y_pred)
    assert loss2 == loss1


def test_ensemble_selection_predict_classifier_yone_hot():
    seed = 3

    X, y = _get_simple_data()

    # convert label id to one-hot vector
    y_one_hot = (
        preprocessing.OneHotEncoder().fit_transform(y.reshape((-1, 1))).toarray()
    )

    estimators = _get_trained_evaluated_estimators("classifier", seed, X, y)

    es = EnsembleSelection(
        estimators=estimators,
        score_metric=log_loss,
        score_direction="min",
        is_base_estimator_proba=True,
        min_estimators=2,
        max_estimators=10,
        with_replacement=True,
        verbose=True,
    )
    es.fit(X, y_one_hot)
    y_pred = es.predict_proba(X)
    loss1 = log_loss(y, y_pred)
    assert loss1 < 0.5  # not so bad

    y_pred = es.predict(X)
    loss2 = np.mean(y_pred == y)
    assert loss2 > 0.5  # not so bad


def test_ensemble_selection_classifier_without_replacement():
    seed = 3

    X, y = _get_simple_data()

    estimators = _get_trained_evaluated_estimators("classifier", seed, X, y)

    es = EnsembleSelection(
        estimators=estimators,
        score_direction="min",
        is_base_estimator_proba=True,
        min_estimators=1,
        max_estimators=5,
        with_replacement=False,
        verbose=False,
    )
    es.fit(X, y)
    y_pred = es.predict_proba(X)
    loss = log_loss(y, y_pred)
    assert loss < 0.2  # not so bad


def test_ensemble_selection_regressor_weights():
    seed = 1
    X, y = _get_simple_data()

    estimators = _get_trained_evaluated_estimators("regressor", seed, X, y)

    # fit ensemble_selector
    es = EnsembleSelection(estimators=estimators)
    weights1 = np.ones((y.shape[0],))
    es.fit(X, y, weights1)
    y_pred = es.predict(X)

    loss = mean_squared_error(y, y_pred)
    assert loss < 0.5


def test_ensemble_selection_regressor_weights_gaussian():
    seed = 1
    X, y = _get_simple_data()
    weights = np.ones((len(y),))

    estimators = []
    estimators.append(("est1", GaussianProcessRegressor(random_state=seed)))
    estimators.append(("est2", GaussianProcessRegressor(random_state=seed + 1)))

    # fit ensemble_selector
    es = EnsembleSelection(estimators=estimators)

    # Gaussian are not compatible with weigths, so weights are ignored
    es.fit(X, y, weights)

    assert es.score(X, y) < 0.5


def test_ensemble_selection_regressor_unfitted_estimators():
    seed = 1
    X, y = _get_simple_data()

    estimators = _get_trained_evaluated_estimators("regressor", seed, X, y)

    # fit ensemble_selector
    es = EnsembleSelection(
        estimators=estimators,
        score_metric=mean_squared_error,
        score_direction="min",
        is_base_estimator_proba=False,
        min_estimators=2,
        max_estimators=10,
        with_replacement=True,
        verbose=True,
    )
    es.fit(X, y)
    y_pred = es.predict(X)

    loss = mean_squared_error(y, y_pred)
    assert loss < 0.5  # not so bad


def test_ensemble_selection_classifier_unamed_estimators():
    X, y = _get_simple_data()

    # instantiate estimators
    estimators = _get_trained_evaluated_estimators("short_classifier")
    unamed_est = [named_est[1] for named_est in estimators]  # remove name

    es = EnsembleSelection(estimators=unamed_est)
    es.fit(X, y)
    loss = es.score(X, y)

    assert loss < 0.5  # not so bad


def test_ensemble_selection_fit_predict_classifier():
    X, y = _get_simple_data()

    # instantiate estimators
    estimators = _get_trained_evaluated_estimators("short_classifier")
    es = EnsembleSelection(estimators=estimators)
    es.fit(X, y)
    pred1 = es.predict(X)

    es = EnsembleSelection(estimators=estimators)
    pred2 = es.fit_predict(X, y)
    assert np.mean(np.abs(pred1 - pred2)) < 0.001
    assert es.score(X, y) < 0.5  # not so bad


def test_ensemble_selection_fit_predict_regressor():
    X, y = _get_simple_data()

    # instantiate estimators
    estimators = _get_trained_evaluated_estimators("short_regressor")
    es = EnsembleSelection(estimators=estimators)
    es.fit(X, y)
    pred1 = es.predict(X)

    es = EnsembleSelection(estimators=estimators)
    pred2 = es.fit_predict(X, y)
    assert np.mean(np.abs(pred1 - pred2)) < 0.001
    assert es.score(X, y) < 0.5  # not so bad


def test_ensemble_selection_custom_metric():
    X, y = _get_simple_data()

    # custom metrics
    def acc(y, y_pred):
        return np.mean(y == np.argmax(y_pred, axis=1))

    def mae(y, y_pred):
        return 1.0 - np.mean(np.abs(y - y_pred))

    # classifier
    estimators = _get_trained_evaluated_estimators("short_classifier")
    es = EnsembleSelection(
        estimators=estimators, score_metric=acc, score_direction="max"
    )
    loss = es.fit(X, y).score(X, y)
    assert loss > 0.5  # not so bad

    # regressor custom metric
    estimators = _get_trained_evaluated_estimators("short_regressor")
    es = EnsembleSelection(
        estimators=estimators, score_metric=mae, score_direction="max"
    )
    loss = es.fit(X, y).score(X, y)
    assert loss > 0.5  # not so bad


def test_ensemble_selection_classifier_not_proba():
    seed = 3

    X, y = _get_simple_data()

    # base estimators may predict proba or directly the class id.
    estimators = []
    estimators.append(("est1", MLPClassifier(random_state=seed)))
    estimators.append(("est2", SVC(random_state=seed, probability=True)))

    # Not using proba is forced
    es = EnsembleSelection(estimators=estimators, is_base_estimator_proba=False)
    loss = es.fit(X, y).score(X, y)

    assert loss < 0.2  # not so bad


def test_ensemble_selection_sample_weights_shape():
    X, y = _get_simple_data()
    estimators = _get_trained_evaluated_estimators("short_classifier")
    estimator = EnsembleSelection(estimators)

    estimator.fit(X, y, sample_weight=np.ones(len(y)))

    with pytest.raises(ValueError):
        estimator.fit(X, y, sample_weight=np.ones(2 * len(y)))

    with pytest.raises(ValueError):
        estimator.fit(X, y, sample_weight=np.ones((len(y), 2)))


def test_ensemble_selection_sparse_data():
    seed = 1
    estimators = _get_trained_evaluated_estimators("short_classifier")
    es = EnsembleSelection(estimators)
    rng = np.random.RandomState(seed)
    X = rng.uniform(size=(40, 3))
    X[X < 0.8] = 0

    X_csr = sparse.csr_matrix(X)
    y = (4 * rng.uniform(size=40)).astype(int)

    es.fit(X_csr, y)


def test_ensemble_selection_panda():
    seed = 1
    try:
        import pandas as pd

        estimators = _get_trained_evaluated_estimators("short_classifier", seed)
        es = EnsembleSelection(estimators)
        X, y = _get_simple_data()

        col_names = ["col_" + str(i) for i in range(X.shape[1])]
        pandas_x = pd.DataFrame(X, columns=col_names)
        pandas_y = pd.DataFrame(y)

        es.fit(pandas_x, pandas_y)

    except ImportError:
        print("Pandas compatibility not tested")


def test_ensemble_selection_3D():
    X = np.zeros((5, 5, 5))
    y = np.zeros((5, 5, 5))
    estimators = _get_trained_evaluated_estimators("short_classifier")
    es = EnsembleSelection(estimators)
    with pytest.raises(ValueError):
        es.fit(X, y)


def test_ensemble_selection_empty():
    X = np.zeros((0,))
    y = np.zeros((0,))
    estimators = _get_trained_evaluated_estimators("short_regressor")
    es = EnsembleSelection(estimators)
    with pytest.raises(ValueError):
        es.fit(X, y)


def test_ensemble_selection_regressor_transform():
    X, y = _get_simple_data()

    # regressor case
    estimators = _get_trained_evaluated_estimators("short_regressor")
    es = EnsembleSelection(estimators)
    es.fit(X, y)
    df = es.decision_function(X)
    shape_if_one_selected = (y.shape[0], 1)
    shape_if_two_selected = (y.shape[0], 2)
    assert df.shape == shape_if_one_selected or df.shape == shape_if_two_selected

    # classifier case
    estimators = _get_trained_evaluated_estimators("short_classifier")
    es = EnsembleSelection(estimators)
    es.fit(X, y)
    df = es.decision_function(X)
    nb_classes = np.max(y) + 1
    shape_if_one_selected = (y.shape[0], 1, nb_classes)
    shape_if_two_selected = (y.shape[0], 2, nb_classes)
    assert df.shape == shape_if_one_selected or df.shape == shape_if_two_selected


def test_ensemble_selection_text_pipeline_classifier():
    seed = 1
    rng = np.random.RandomState(seed)

    n_samples = 30
    X = rng.choice(np.array(["aa", "bb", "cc"], dtype=object), size=n_samples)

    # y = rng.normal(size=n_samples) # regression version
    y = rng.randint(3, size=n_samples)  # classification version

    pip1 = Pipeline(
        steps=[
            ("tfidfvectorizer", TfidfVectorizer()),
            ("logisticregression", LogisticRegression(C=0.1)),
        ]
    )
    pip2 = Pipeline(
        steps=[
            ("tfidfvectorizer", TfidfVectorizer()),
            ("logisticregression", LogisticRegression(C=1.0)),
        ]
    )
    estimators = [("est1", pip1), ("est2", pip2)]

    estimator = EnsembleSelection(estimators=estimators)
    estimator.fit(X, y)


def test_ensemble_selection_no_estimators():
    X, y = _get_simple_data()

    # fit ensemble_selector
    es = EnsembleSelection([])
    with pytest.raises(ValueError):
        es.fit(X, y)


def test_ensemble_selection_clustering():
    X, y = _get_simple_data()
    estimators = []
    estimators.append(("est1", KMeans(n_clusters=2)))
    estimators.append(("est2", DBSCAN(eps=0.5)))
    es = EnsembleSelection(estimators=estimators)

    # Not compatible compatible with clustering
    with pytest.raises(ValueError):
        es.fit(X, y)


def test_ensemble_selection_sequential():
    X, y = _get_simple_data()
    estimators = _get_trained_evaluated_estimators("short_regressor")
    es = EnsembleSelection(estimators=estimators, n_jobs=1)
    es.fit(X, y)
    score = es.score(X, y)
    assert score < 0.5

""" test the label propagation module """

import numpy as np
import pytest
import scipy.sparse as sp

from scipy.sparse import issparse
from sklearn.semi_supervised import _label_propagation as label_propagation
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.model_selection import train_test_split
from sklearn.neighbors import NearestNeighbors
from sklearn.datasets import make_classification
from sklearn.exceptions import ConvergenceWarning
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal
from sklearn.utils._testing import _convert_container

ESTIMATOR_CLS = [label_propagation.LabelPropagation,
                 label_propagation.LabelSpreading]


PARAMETERS = [{'kernel': 'rbf'}, {'kernel': 'knn', 'n_neighbors': 2},
              {'kernel': lambda x, y: rbf_kernel(x, y, gamma=20)}]


CONSTRUCTOR_TYPES = ("array", "sparse_csr", "sparse_csc")


@pytest.fixture
def estimators():
    estimators_ = []
    for estimator in ESTIMATOR_CLS:
        for params in PARAMETERS:
            estimators_.append(estimator(**params))
    return estimators_


@pytest.mark.parametrize("constructor_type", CONSTRUCTOR_TYPES)
def test_fit_transduction(constructor_type, estimators):
    for estimator in estimators:
        samples = _convert_container(
            [[1., 0.], [0., 2.], [1., 3.]], constructor_type)
        labels = [0, 1, -1]
        clf = estimator.fit(samples, labels)
        assert clf.transduction_[2] == 1


@pytest.mark.parametrize("constructor_type", CONSTRUCTOR_TYPES)
def test_distribution(constructor_type, estimators):
    for estimator in estimators:
        samples = _convert_container(
            [[1., 0.], [0., 1.], [1., 1.]], constructor_type)
        labels = [0, 1, -1]
        clf = estimator.fit(samples, labels)
        if estimator.kernel != 'knn':
            # unstable test; changes in k-NN ordering break it
            assert_array_almost_equal(np.asarray(clf.label_distributions_[2]),
                                      np.array([.5, .5]), 2)


@pytest.mark.parametrize("constructor_type", CONSTRUCTOR_TYPES)
def test_predict(constructor_type, estimators):
    for estimator in estimators:
        samples = _convert_container(
            [[1., 0.], [0., 2.], [1., 3.]], constructor_type)
        labels = [0, 1, -1]
        clf = estimator.fit(samples, labels)
        assert_array_equal(clf.predict([[0.5, 2.5]]), np.array([1]))


@pytest.mark.parametrize("constructor_type", CONSTRUCTOR_TYPES)
def test_predict_proba(constructor_type, estimators):
    for estimator in estimators:
        samples = _convert_container(
            [[1., 0.], [0., 1.], [1., 2.5]], constructor_type)
        labels = [0, 1, -1]
        clf = estimator.fit(samples, labels)
        assert_array_almost_equal(clf.predict_proba([[1., 1.]]),
                                np.array([[0.5, 0.5]]))


@pytest.mark.parametrize("constructor_type", CONSTRUCTOR_TYPES)
def test_label_spreading_closed_form(constructor_type):
    n_classes = 2
    X, y = make_classification(n_classes=n_classes, n_samples=200,
                               random_state=0)
    X = _convert_container(X, constructor_type)
    y[::3] = -1
    clf = label_propagation.LabelSpreading().fit(X, y)
    # adopting notation from Zhou et al (2004):
    S = clf._build_graph()
    Y = np.zeros((len(y), n_classes + 1))
    Y[np.arange(len(y)), y] = 1
    Y = Y[:, :-1]
    for alpha in [0.1, 0.3, 0.5, 0.7, 0.9]:
        expected = np.dot(np.linalg.inv(np.eye(len(S)) - alpha * S), Y)
        expected /= expected.sum(axis=1)[:, np.newaxis]
        clf = label_propagation.LabelSpreading(max_iter=10000, alpha=alpha)
        clf.fit(X, y)
        assert_array_almost_equal(expected, clf.label_distributions_, 4)


@pytest.mark.parametrize("constructor_type", CONSTRUCTOR_TYPES)
def test_label_propagation_closed_form(constructor_type):
    n_classes = 2
    X, y = make_classification(n_classes=n_classes, n_samples=200,
                               random_state=0)
    X = _convert_container(X, constructor_type)
    y[::3] = -1
    Y = np.zeros((len(y), n_classes + 1))
    Y[np.arange(len(y)), y] = 1
    unlabelled_idx = Y[:, (-1,)].nonzero()[0]
    labelled_idx = (Y[:, (-1,)] == 0).nonzero()[0]

    clf = label_propagation.LabelPropagation(max_iter=10000,
                                             gamma=0.1)
    clf.fit(X, y)
    # adopting notation from Zhu et al 2002
    T_bar = clf._build_graph()
    Tuu = T_bar[tuple(np.meshgrid(unlabelled_idx, unlabelled_idx,
                                  indexing='ij'))]
    Tul = T_bar[tuple(np.meshgrid(unlabelled_idx, labelled_idx,
                                  indexing='ij'))]
    Y = Y[:, :-1]
    Y_l = Y[labelled_idx, :]
    Y_u = np.dot(np.dot(np.linalg.inv(np.eye(Tuu.shape[0]) - Tuu), Tul), Y_l)

    expected = Y.copy()
    expected[unlabelled_idx, :] = Y_u
    expected /= expected.sum(axis=1)[:, np.newaxis]

    assert_array_almost_equal(expected, clf.label_distributions_, 4)


@pytest.mark.parametrize("constructor_type", CONSTRUCTOR_TYPES)
def test_valid_alpha(constructor_type):
    n_classes = 2
    X, y = make_classification(n_classes=n_classes, n_samples=200,
                               random_state=0)
    X = _convert_container(X, constructor_type)
    for alpha in [-0.1, 0, 1, 1.1, None]:
        with pytest.raises(ValueError):
            label_propagation.LabelSpreading(alpha=alpha).fit(X, y)


@pytest.mark.parametrize("constructor_type", CONSTRUCTOR_TYPES)
def test_convergence_speed(constructor_type):
    # This is a non-regression test for #5774
    X = _convert_container([[1., 0.], [0., 1.], [1., 2.5]], constructor_type)
    y = np.array([0, 1, -1])
    mdl = label_propagation.LabelSpreading(kernel='rbf', max_iter=5000)
    mdl.fit(X, y)

    # this should converge quickly:
    assert mdl.n_iter_ < 10
    assert_array_equal(mdl.predict(X), [0, 1, 1])


@pytest.mark.parametrize("constructor_type", CONSTRUCTOR_TYPES)
def test_convergence_warning(constructor_type):
    # This is a non-regression test for #5774
    X = _convert_container([[1., 0.], [0., 1.], [1., 2.5]], constructor_type)
    y = np.array([0, 1, -1])
    mdl = label_propagation.LabelSpreading(kernel='rbf', max_iter=1)
    warn_msg = ('max_iter=1 was reached without convergence.')
    with pytest.warns(ConvergenceWarning, match=warn_msg):
        mdl.fit(X, y)
    assert mdl.n_iter_ == mdl.max_iter

    mdl = label_propagation.LabelPropagation(kernel='rbf', max_iter=1)
    with pytest.warns(ConvergenceWarning, match=warn_msg):
        mdl.fit(X, y)
    assert mdl.n_iter_ == mdl.max_iter

    mdl = label_propagation.LabelSpreading(kernel='rbf', max_iter=500)
    with pytest.warns(None) as record:
        mdl.fit(X, y)
    assert len(record) == 0

    mdl = label_propagation.LabelPropagation(kernel='rbf', max_iter=500)
    with pytest.warns(None) as record:
        mdl.fit(X, y)

    assert len(record) == 0


@pytest.mark.parametrize("constructor_type", CONSTRUCTOR_TYPES)
@pytest.mark.parametrize("LabelPropagationCls", ESTIMATOR_CLS)
def test_label_propagation_non_zero_normalizer(
    LabelPropagationCls, constructor_type
):
    # check that we don't divide by zero in case of null normalizer
    # non-regression test for
    # https://github.com/scikit-learn/scikit-learn/pull/15946
    # https://github.com/scikit-learn/scikit-learn/issues/9292
    X = _convert_container([[100., 100.], [100., 100.], [
                    0., 0.], [0., 0.]], constructor_type)
    y = np.array([0, 1, -1, -1])
    mdl = LabelPropagationCls(kernel='knn',
                              max_iter=100,
                              n_neighbors=1)
    with pytest.warns(None) as record:
        mdl.fit(X, y)
    if constructor_type is not sp.dok_matrix:  # fails for dok_matrix
        assert len(record) == 0


@pytest.mark.parametrize("constructor_type", CONSTRUCTOR_TYPES)
def test_predict_sparse_callable_kernel(constructor_type):
    # This is a non-regression test for #15866

    # Custom sparse kernel (top-K RBF)
    def topk_rbf(X, Y=None, n_neighbors=10, gamma=1e-5):
        nn = NearestNeighbors(n_neighbors=10, metric='euclidean', n_jobs=-1)
        nn.fit(X)
        W = -1 * nn.kneighbors_graph(Y, mode='distance').power(2) * gamma
        np.exp(W.data, out=W.data)
        assert issparse(W)
        return W.T

    n_classes = 4
    n_samples = 500
    n_test = 10
    X, y = make_classification(n_classes=n_classes,
                               n_samples=n_samples,
                               n_features=20,
                               n_informative=20,
                               n_redundant=0,
                               n_repeated=0,
                               random_state=0)
    X = _convert_container(X, constructor_type)
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        test_size=n_test,
                                                        random_state=0)

    model = label_propagation.LabelSpreading(kernel=topk_rbf)
    model.fit(X_train, y_train)
    assert model.score(X_test, y_test) >= 0.9

    model = label_propagation.LabelPropagation(kernel=topk_rbf)
    model.fit(X_train, y_train)
    assert model.score(X_test, y_test) >= 0.9

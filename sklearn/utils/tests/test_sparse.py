import numpy as np
import pytest
from scipy.sparse import csc_array, csc_matrix, csr_array, csr_matrix

import sklearn


@pytest.mark.parametrize(
    ["sparse_interface", "x", "result_type"],
    [
        ("sparray", csr_array([[1, 2, 3]]), csr_array),
        ("sparray", csr_matrix([[1, 2, 3]]), csr_array),
        ("sparray", csc_array([[1, 2, 3]]), csc_array),
        ("sparray", csc_matrix([[1, 2, 3]]), csc_array),
        ("spmatrix", csr_array([[1, 2, 3]]), csr_matrix),
        ("spmatrix", csr_matrix([[1, 2, 3]]), csr_matrix),
        ("spmatrix", csc_array([[1, 2, 3]]), csc_matrix),
        ("spmatrix", csc_matrix([[1, 2, 3]]), csc_matrix),
    ],
)
def test_interface_if_sparse(sparse_interface, x, result_type):
    with sklearn.config_context(sparse_interface=sparse_interface):
        result = sklearn.utils._select_interface_if_sparse(x)
        result = sklearn.utils._as_sparse(x)
        assert isinstance(result, result_type)


@pytest.mark.parametrize(
    ["sparse_interface", "x", "result_type"],
    [
        ("sparray", np.array([[1, 2, 3]]), np.ndarray),
        ("spmatrix", np.array([[1, 2, 3]]), np.ndarray),
    ],
)
def test_ndarray_select_interface_if_sparse(sparse_interface, x, result_type):
    with sklearn.config_context(sparse_interface=sparse_interface):
        result = sklearn.utils._select_interface_if_sparse(x)
        assert isinstance(result, result_type)


@pytest.mark.parametrize("sparse_interface", ["sparray", "spmatrix"])
def test_ndarray_as_sparse_raises(sparse_interface):
    with sklearn.config_context(sparse_interface=sparse_interface):
        with pytest.raises(TypeError, match="sparse container"):
            sklearn.utils._as_sparse(np.array([[1, 2, 3]]))


@pytest.mark.parametrize(
    ["sparse_interface", "result_type"],
    [("sparray", csr_array), ("spmatrix", csr_matrix)],
)
def test_transform_returns_sparse(sparse_interface, result_type):
    corpus = [
        "This is the first document.",
        "This document is the second document.",
        "And this is the third one.",
        "Is this the first document?",
    ]
    with sklearn.config_context(sparse_interface=sparse_interface):
        vectorizer = sklearn.feature_extraction.text.CountVectorizer()
        X = vectorizer.fit_transform(corpus)
        assert isinstance(X, result_type)


@pytest.mark.parametrize(
    ["sparse_interface", "result_type"],
    [("sparray", csr_array), ("spmatrix", csr_matrix)],
)
def test_function_returns_sparse(sparse_interface, result_type):
    with sklearn.config_context(sparse_interface=sparse_interface):
        X, y = sklearn.datasets.make_regression(n_features=2, random_state=0)
        X = sklearn.manifold._locally_linear.barycenter_kneighbors_graph(X, 1)
        assert isinstance(X, result_type)


@pytest.mark.parametrize(
    ["sparse_interface", "result_type"],
    [("sparray", csr_array), ("spmatrix", csr_matrix)],
)
def test_estimator_property_sparse(sparse_interface, result_type):
    with sklearn.config_context(sparse_interface=sparse_interface):
        X, y = sklearn.datasets.make_regression(n_features=2, random_state=0)
        regr = sklearn.linear_model.ElasticNet(random_state=0)
        regr.fit(X, y)
        # check spec_coeff property
        assert isinstance(regr.sparse_coef_, result_type)

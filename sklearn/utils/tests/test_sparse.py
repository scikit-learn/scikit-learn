import numpy as np
import pytest
import scipy as sp
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
def test_align_api_if_sparse(sparse_interface, x, result_type):
    with sklearn.config_context(sparse_interface=sparse_interface):
        result = sklearn.utils._align_api_if_sparse(x)
        assert isinstance(result, result_type)


@pytest.mark.parametrize(
    ["sparse_interface", "x", "result_type"],
    [
        ("sparray", np.array([[1, 2, 3]]), np.ndarray),
        ("spmatrix", np.array([[1, 2, 3]]), np.ndarray),
    ],
)
def test_ndarray_align_api_if_sparse(sparse_interface, x, result_type):
    with sklearn.config_context(sparse_interface=sparse_interface):
        result = sklearn.utils._align_api_if_sparse(x)
        assert isinstance(result, result_type)


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


@pytest.mark.parametrize(
    "constructor",
    [
        sp.sparse.bsr_array, sp.sparse.bsr_matrix,
        sp.sparse.csc_array, sp.sparse.csc_matrix,
        sp.sparse.csr_array, sp.sparse.csr_matrix,
        sp.sparse.coo_array, sp.sparse.coo_matrix,
        sp.sparse.dia_array, sp.sparse.dia_matrix,
        sp.sparse.dok_array, sp.sparse.dok_matrix,
        sp.sparse.lil_array, sp.sparse.lil_matrix,
    ]
)
def test_estimator_property_sparse(constructor):
    A = constructor(np.array([[1., 2., 3.], [3., 2., 1.]]))
    sklearn.utils._sparse._ensure_sparse_index_int32(A)

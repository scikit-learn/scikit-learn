import scipy.sparse as sp
import numpy as np

from sklearn.metrics import euclidean_distances
from sklearn.datasets import make_low_rank_matrix
from sklearn.random_projection import SparseRandomProjection
from sklearn.random_projection import johnson_lindenstrauss_min_dim
from sklearn.random_projection import hashing_dot

from sklearn.utils.testing import assert_raise_message
from numpy.testing import assert_array_equal
from nose.tools import assert_equal
from nose.tools import assert_almost_equal
from nose.tools import assert_raises


def assert_lower(a, b, details=None):
    message = "%r is not lower than %r" % (a, b)
    if details is not None:
        message += ": " + details
    assert a < b, message


n_samples, n_features = (10, 10000)
data = make_low_rank_matrix(n_samples=n_samples, n_features=n_features,
                            random_state=0)


def test_invalid_jl_domain():
    assert_raises(ValueError, johnson_lindenstrauss_min_dim, 100, 1.1)
    assert_raises(ValueError, johnson_lindenstrauss_min_dim, 100, 0.0)
    assert_raises(ValueError, johnson_lindenstrauss_min_dim, 100, -0.1)


def test_hashing_dot():
    n_components = 100
    density = 0.01
    projected = hashing_dot(data, n_components, density, seed=0,
                            dense_output=True)

    assert_equal(projected.shape, (n_samples, n_components))
    assert_almost_equal(projected.mean(), 0.0, 2)


def test_sparse_random_project_invalid_input():
    assert_raises(ValueError, SparseRandomProjection(density=0.8).fit, data)
    assert_raises(ValueError, SparseRandomProjection().fit, [0, 1, 2])


def test_input_dimension_inconsistency():
    expected_msg = ("n_components=20000 should be smaller "
                    "than n_features=10000")
    rp = SparseRandomProjection(n_components=20000)
    assert_raise_message(ValueError, expected_msg, rp.fit, data)


def test_too_many_samples_to_find_a_safe_embedding():
    data = make_low_rank_matrix(n_samples=1000, n_features=100)
    expected_msg = (
        'eps=0.100000 and n_samples=1000 lead to a target dimension'
        ' of 5920 which is larger than the original space with'
        ' n_features=100')
    rp = SparseRandomProjection(n_components='auto')
    assert_raise_message(ValueError, expected_msg, rp.fit, data)


def test_sparse_random_projection_dimensions():
    rp = SparseRandomProjection(random_state=0).fit(data)

    # the number of components is adjusted from the shape of the training set
    assert_equal(rp.n_components, 'auto')
    assert_equal(rp.n_components_, 1973)
    assert_equal(rp.density, 'auto')
    assert_equal(rp.density_, 0.01)
    assert_equal(rp.components_.shape, (1973, n_features))

    projected_1 = rp.transform(data)
    assert_equal(projected_1.shape, (n_samples, 1973))

    # once the RP is 'fitted' the projection is always the same
    projected_2 = rp.transform(data)
    assert_array_equal(projected_1, projected_2)

    # fit transform with same random seed will lead to the same results
    rp2 = SparseRandomProjection(random_state=0)
    projected_3 = rp2.fit_transform(data)
    assert_array_equal(projected_1, projected_3)

    # it also possible to fix the number of components and the density level
    rp = SparseRandomProjection(n_components=100, density=0.001,
                                random_state=0)
    projected = rp.fit_transform(data)
    assert_equal(projected.shape, (n_samples, 100))
    assert_equal(rp.components_.shape, (100, n_features))
    assert_lower(rp.components_.nnz, 1100)  # close to 1% density
    assert_lower(900, rp.components_.nnz)  # close to 1% density


def test_sparse_projection_embedding_quality():
    eps = 0.1
    rp = SparseRandomProjection(n_components='auto', density='auto', eps=eps,
                                random_state=0)
    projected = rp.fit_transform(data)

    original_distances = euclidean_distances(data, squared=True).ravel()
    projected_distances = euclidean_distances(projected, squared=True).ravel()
    non_identical = original_distances != 0.0

    # remove 0 distances to avoid division by 0
    original_distances = original_distances[non_identical]
    projected_distances = projected_distances[non_identical]

    distances_ratio = projected_distances / original_distances

    # check that the automatically tuned values for the density respect the
    # contract for eps: pairwise distances are preserved according to the
    # Johnson Lindenstrauss bound
    assert_lower(distances_ratio.max(), 1 + eps)
    assert_lower(1 - eps, distances_ratio.min())


def test_output_representation():
    # by default dense ndarray output is enforced
    rp = SparseRandomProjection(n_components=100, random_state=0).fit(data)
    assert isinstance(rp.transform(data), np.ndarray)

    sparse_data = sp.csr_matrix(data)
    assert isinstance(rp.transform(sparse_data), np.ndarray)

    # this behavior can be disabled:
    rp = SparseRandomProjection(n_components=100, dense_output=False,
                                random_state=0).fit(data)
    # output for dense input will stay dense:
    assert isinstance(rp.transform(data), np.ndarray)

    # ouput for sparse output will be sparse:
    assert sp.isspmatrix_csr(rp.transform(sparse_data))

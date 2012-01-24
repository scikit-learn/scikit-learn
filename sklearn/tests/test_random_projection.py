import unittest
import scipy.sparse as sp
import numpy as np

from sklearn.metrics import euclidean_distances
from sklearn.random_projection import SparseRandomProjection
from sklearn.random_projection import johnson_lindenstrauss_min_dim
from sklearn.random_projection import hashing_dot

from sklearn.utils.testing import assert_raise_message
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal
from nose.tools import assert_almost_equal
from nose.tools import assert_raises
from nose.tools import assert_false


def assert_lower(a, b, details=None):
    message = "%r is not lower than %r" % (a, b)
    if details is not None:
        message += ": " + details
    assert a < b, message


# Make a some random data with uniformly located non zero entries with
# gaussian distributed values
def make_sparse_random_data(n_samples, n_features, n_nonzeros):
    rng = np.random.RandomState(0)
    data_coo = sp.coo_matrix(
        (rng.randn(n_nonzeros),
         (rng.randint(n_samples, size=n_nonzeros),
          rng.randint(n_features, size=n_nonzeros))),
        shape=(n_samples, n_features))
    return data_coo.toarray(), data_coo.tocsr()


n_samples, n_features = (10, 1000)
n_nonzeros = n_samples * n_features / 100.
data, data_csr = make_sparse_random_data(n_samples, n_features, n_nonzeros)


def test_invalid_jl_domain():
    assert_raises(ValueError, johnson_lindenstrauss_min_dim, 100, 1.1)
    assert_raises(ValueError, johnson_lindenstrauss_min_dim, 100, 0.0)
    assert_raises(ValueError, johnson_lindenstrauss_min_dim, 100, -0.1)


def test_hashing_dot():
    n_components = 300
    density = 0.01

    # random dot with dense array input (and output)
    projected_array_array = hashing_dot(data, n_components, density, seed=0)
    assert_equal(projected_array_array.shape, (n_samples, n_components))
    assert_almost_equal(projected_array_array.mean(), 0.0, 2)

    # random dot with CSR input and dense array output
    projected_csr_array = hashing_dot(data_csr, n_components, density, seed=0,
                                      dense_output=True)
    assert_equal(projected_csr_array.shape, (n_samples, n_components))

    # random dot with CSR input and sparse matrix output
    projected_csr_coo = hashing_dot(data_csr, n_components, density, seed=0,
                                      dense_output=False)
    assert_equal(projected_csr_coo.shape, (n_samples, n_components))

    # check that the projections are identical
    assert_array_almost_equal(projected_array_array, projected_csr_array)
    assert_array_almost_equal(projected_array_array,
                              projected_csr_coo.toarray())


class MaterializedRandomProjection(unittest.TestCase):

    materialize = True

    def setUp(self):
        self.rp = SparseRandomProjection(n_components='auto',
                                         materialize=self.materialize,
                                         random_state=0)

    def test_sparse_random_project_invalid_input(self):
        rp_08 = SparseRandomProjection(
            density=0.8, materialize=self.materialize)
        assert_raises(ValueError, rp_08.fit, data)
        assert_raises(ValueError, self.rp.fit, [0, 1, 2])

    def test_input_dimension_inconsistency(self):
        expected_msg = ("n_components=20000 should be smaller "
                        "than n_features=1000")
        rp = SparseRandomProjection(
            n_components=20000, materialize=self.materialize)
        assert_raise_message(ValueError, expected_msg, rp.fit, data)

    def test_too_many_samples_to_find_a_safe_embedding(self):
        data, _ = make_sparse_random_data(1000, 100, 1000)
        rp = SparseRandomProjection(
            n_components='auto', eps=0.1, materialize=self.materialize)
        expected_msg = (
            'eps=0.100000 and n_samples=1000 lead to a target dimension'
            ' of 5920 which is larger than the original space with'
            ' n_features=100')
        assert_raise_message(ValueError, expected_msg, rp.fit, data)

    def test_sparse_random_projection_dimensions(self):
        rp = self.rp.fit(data)

        # the number of components is adjusted from the shape of the training
        # set
        assert_equal(rp.n_components, 'auto')
        assert_equal(rp.n_components_, 110)
        assert_equal(rp.density, 'auto')
        assert_almost_equal(rp.density_, 0.03, 2)
        if self.materialize:
            assert_equal(rp.components_.shape, (110, n_features))

        projected_1 = rp.transform(data)
        assert_equal(projected_1.shape, (n_samples, 110))

        # once the RP is 'fitted' the projection is always the same
        projected_2 = rp.transform(data)
        assert_array_equal(projected_1, projected_2)

        # fit transform with same random seed will lead to the same results
        rp2 = SparseRandomProjection(materialize=self.materialize,
                                     random_state=0)
        projected_3 = rp2.fit_transform(data)
        assert_array_equal(projected_1, projected_3)

        # it is also possible to fix the number of components and the density
        # level
        rp = SparseRandomProjection(n_components=100, density=0.001,
                                    random_state=0,
                                    materialize=self.materialize)
        projected = rp.fit_transform(data)
        assert_equal(projected.shape, (n_samples, 100))
        if self.materialize:
            assert_equal(rp.components_.shape, (100, n_features))
            assert_lower(rp.components_.nnz, 110)  # close to 1% density
            assert_lower(90, rp.components_.nnz)  # close to 1% density
        else:
            assert_false(hasattr(rp, 'components_'))

    def test_sparse_projection_embedding_quality(self):
        eps = 0.1
        data, _ = make_sparse_random_data(10, 5000, 10000)
        rp = SparseRandomProjection(n_components='auto', density='auto',
                                    materialize=self.materialize, eps=eps,
                                    random_state=0)
        projected = rp.fit_transform(data)

        original_distances = euclidean_distances(data, squared=True)
        original_distances = original_distances.ravel()

        projected_distances = euclidean_distances(projected, squared=True)
        projected_distances = projected_distances.ravel()

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

    def test_output_representation(self):
        # by default dense ndarray output is enforced
        self.rp.fit(data)
        assert isinstance(self.rp.transform(data), np.ndarray)

        sparse_data = sp.csr_matrix(data)
        assert isinstance(self.rp.transform(sparse_data), np.ndarray)

        # this behavior can be disabled:
        rp = SparseRandomProjection(n_components=100, dense_output=False,
                                    materialize=self.materialize,
                                    random_state=0)
        rp = rp.fit(data)
        # output for dense input will stay dense:
        assert isinstance(rp.transform(data), np.ndarray)

        # ouput for sparse output will be sparse:
        assert sp.issparse(rp.transform(sparse_data))


class HashingSparseRandomProjection(MaterializedRandomProjection):
    """Same tests as above but without allocating the projection matrix"""

    materialize = False

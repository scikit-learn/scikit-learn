from sklearn.metrics import euclidean_distances
from sklearn.datasets import make_low_rank_matrix
from sklearn.random_projection import SparseRandomProjection

from numpy.testing import assert_array_equal
from nose.tools import assert_equal
from nose.tools import assert_raises
from nose.tools import assert_almost_equal


def assert_lower(a, b, details=None):
    message = "%r is not lower than %r" % (a, b)
    if details is not None:
        message += ": " + details
    assert a < b, message


n_samples, n_features = (10, 10000)
data = make_low_rank_matrix(n_samples=n_samples, n_features=n_features,
                            random_state=0)


def test_sparse_random_project_invalid_input():
    assert_raises(ValueError, SparseRandomProjection(density=0.8).fit, data)
    assert_raises(ValueError, SparseRandomProjection().fit, [0, 1, 2])


def test_too_many_samples_to_find_a_safe_embedding():
    # check that a warning is raised when no safe dim reduction is possible

    # install a monkey patch to collect the warnings
    collected_warnings = set()
    def collecting_warn(msg):
        collected_warnings.add(msg)
    import warnings
    warn_original = warnings.warn
    warnings.warn = collecting_warn

    try:

      data = make_low_rank_matrix(n_samples=1000, n_features=100)
      rp = SparseRandomProjection().fit(data)
      assert_equal(rp.n_components, 100)  # no safe dim reduction possible

      expected_warnings = set([
          'eps=0.100000 and n_samples=1000 lead to a target dimension'
          ' of 5920 which is larger than the original space with'
          ' n_features=100'])
      assert_equal(collected_warnings, expected_warnings)

    finally:
        # restore the warn function
        warnings.warn = warn_original


def test_sparse_random_projection_dimensions():
    rp = SparseRandomProjection(random_state=0).fit(data)

    # the number of components is adjusted from the shape of the training set
    assert_equal(rp.n_components, 1973)
    assert_equal(rp.density, 0.01)
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
    assert_equal(rp.components_.nnz, 960)  # close to 1% density


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
    # contract for eps: pairwise distances are preserved
    assert_almost_equal(distances_ratio.mean(), 1.00, 2)
    assert_almost_equal(distances_ratio.std(), 0.03, 2)

    # check the Johnson Lindenstrauss bound
    assert_lower(distances_ratio.max(), 1 + eps)
    assert_lower(1 - eps, distances_ratio.min())

from sklearn.datasets import make_low_rank_matrix
from sklearn.random_projection import SparseRandomProjection

from numpy.testing import assert_array_equal
from nose.tools import assert_equals


def test_sparse_random_projection():
    data = make_low_rank_matrix(n_samples=10, n_features=10000)

    rp = SparseRandomProjection(random_state=0).fit(data)

    # the number of components is adjusted from the shape of the training set
    assert_equals(rp.n_components, 1973)
    assert_equals(rp.density, 0.01)
    assert_equals(rp.components_.shape, (1973, 10000))

    # once the RP is 'fitted' the projection is always the same
    projected_1 = rp.transform(data)
    projected_2 = rp.transform(data)
    assert_array_equal(projected_1, projected_2)

    # TODO: check the Johnson Lindenstrauss embedding

    # fit transform with same random seed will lead to the same results
    rp2 = SparseRandomProjection(random_state=0)
    projected_3 = rp2.fit_transform(data)
    assert_array_equal(projected_1, projected_3)

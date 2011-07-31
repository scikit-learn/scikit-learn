import numpy as np
from numpy.testing import assert_equal

from .. import generate_sparse_coded_signal


def test_generate_sparse_coded_signal():
    Y, D, X = generate_sparse_coded_signal(n_samples=5, n_components=8,
                                           n_features=10, n_nonzero_coefs=3,
                                           random_state=0)
    assert_equal(Y.shape, (10, 5), 'Data shape mismatch')
    assert_equal(D.shape, (10, 8), 'Dictionary shape mismatch')
    assert_equal(X.shape, (8, 5), 'Code shape mismatch')
    for col in X.T:
        assert_equal(len(np.flatnonzero(col)), 3, 'Non-zero coefs mismatch')
    assert_equal(np.dot(D, X), Y)

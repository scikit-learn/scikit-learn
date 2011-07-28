import numpy as np
from numpy.testing import assert_equal

from .. import generate_sparse_coded_signal


def test_generate_sparse_coded_signal():
    X, W, H = generate_sparse_coded_signal(n_samples=5, n_components=8,
                                           n_features=10, n_nonzero_coefs=3,
                                           random_state=0)
    assert_equal(X.shape, (5, 10), 'Data shape mismatch')
    assert_equal(W.shape, (5, 8), 'Code shape mismatch')
    assert_equal(H.shape, (8, 10), 'Dictionary shape mismatch')
    for line in W:
        assert_equal(len(line.nonzero()[0]), 3, 'Non-zero coefs mismatch')
    assert_equal(np.dot(W, H), X)

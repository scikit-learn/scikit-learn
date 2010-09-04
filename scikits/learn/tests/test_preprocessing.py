import numpy as np

from numpy.testing import assert_array_almost_equal, assert_array_equal, assert_

from scikits.learn.preprocessing import Scaler

def test_scaler():
    """Test scaling of dataset along all axis
    """

    X = np.random.randn(4, 5)

    scaler = Scaler(axis=1)
    X_scaled = scaler.fit(X).transform(X, copy=False)
    assert_array_almost_equal(X_scaled.mean(axis=1), 4*[0.0])
    assert_array_almost_equal(X_scaled.std(axis=1), 4*[1.0])
    # Check that X has not been copied
    assert_(X_scaled is X)
    

    scaler = Scaler(axis=0, with_std=False)
    X_orig = X.copy()
    X_scaled = scaler.fit(X).transform(X, copy=True)
    assert_array_almost_equal(X_scaled.mean(axis=0), 5*[0.0])
    # Check that the data hasn't been modified
    assert_array_equal(X, X_orig)

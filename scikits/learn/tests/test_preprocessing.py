import numpy as np

from numpy.testing import assert_array_almost_equal, assert_array_equal, assert_

#from ..preprocessing import Scaler
from scikits.learn.preprocessing import Scaler, scale

def test_scaler():
    """Test scaling of dataset along all axis
    """
    X = np.random.randn(4, 5)

    scaler = Scaler()
    X_scaled = scaler.fit(X).transform(X, copy=False)
    assert_array_almost_equal(X_scaled.mean(axis=0), 5*[0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), 5*[1.0])
    # Check that X has not been copied
    assert_(X_scaled is X)

    X_scaled = scaler.fit(X).transform(X, copy=True)
    assert_array_almost_equal(X_scaled.mean(axis=0), 5*[0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), 5*[1.0])
    # Check that X has not been copied
    assert_(X_scaled is not X)

    X_scaled = scale(X, axis=1, with_std=False)
    assert_array_almost_equal(X_scaled.mean(axis=1), 4*[0.0])
    X_scaled = scale(X, axis=1, with_std=True)
    assert_array_almost_equal(X_scaled.std(axis=1), 4*[1.0])
    # Check that the data hasn't been modified


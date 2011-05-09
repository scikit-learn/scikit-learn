from .. import datasets, perceptron
import numpy as np
from numpy.testing import assert_array_equal

digits = datasets.load_digits()

def test_perceptron():
    '''
    Test perceptrons on digit recognition task.
    '''

    X = digits.data
    y = digits.target

    # Test self-returning from fit etc. while we're at it
    clf_batch = perceptron.Perceptron() \
              . fit(X[1:], y[1:], averaged=True)

    n_features = X.shape[1]
    n_samples  = X.shape[0] - 1
    half = n_samples // 2
    clf_online = perceptron.Perceptron() \
               . partial_setup(n_features, len(np.unique(y)),
                               averaged=True) \
               . partial_fit(X[: half], y[: half]) \
               . partial_fit(X[n_samples - half : -1],
                             y[n_samples - half : -1])

    assert_array_equal(clf_batch.predict(y[-1]),
                       clf_online.predict(y[-1]))


if __name__ == '__main__':
    import nose
    nose.runmodule()

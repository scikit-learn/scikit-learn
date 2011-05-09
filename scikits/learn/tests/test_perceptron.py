from .. import datasets, perceptron
import numpy as np
from numpy.testing import assert_array_equal

digits = datasets.load_digits()

def test_perceptron():
    '''
    Test perceptrons on digit recognition task.
    '''

    # Test self-returning from fit etc. while we're at it
    clf_batch = perceptron.Perceptron() \
              . fit(digits.data[1:], digits.target[1:], averaged=True)

    n_features = digits.data.shape[1]
    n_samples = digits.data.shape[0] - 1
    half = n_samples // 2
    clf_online = perceptron.Perceptron() \
               . partial_setup(n_features, len(np.unique(digits.target)),
                               averaged=True) \
               . partial_fit(digits.data[: half], digits.target[: half]) \
               . partial_fit(digits.data[n_samples - half : -1],
                             digits.target[n_samples - half : -1])

    assert_array_equal(clf_batch.predict(digits.data[-1]),
                       clf_online.predict(digits.data[-1]))


if __name__ == '__main__':
    import nose
    nose.runmodule()

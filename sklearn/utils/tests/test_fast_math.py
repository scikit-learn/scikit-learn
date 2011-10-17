import numpy as np

from sklearn.utils import fast_math
from sklearn.utils import check_random_state

def test_fast_log(random_state=1):
    """ Check that the fast_log gives reasonnable answers. """
    random_state = check_random_state(random_state)
    for i in range(10):
        l = random_state.normal()
        e = np.exp(l)
        np.testing.assert_array_almost_equal(l, fast_math.py_fast_log(e),
                    decimal=1)


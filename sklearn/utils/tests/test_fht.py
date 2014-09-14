import numpy as np
import numpy.testing as npt
import nose.tools as nt
from scipy.linalg import hadamard
from sklearn.utils.cyfht import fht as cyfht


def test_wikipedia_example():
    input_ = np.array([1, 0, 1, 0, 0, 1, 1, 0], dtype=np.float64)
    copy = input_.copy()
    H = hadamard(8)
    cyfht(input_)
    npt.assert_array_equal(np.dot(copy, H), input_)


def test_exception_when_input_not_power_two():
    nt.assert_raises(ValueError, cyfht, np.zeros(9, dtype=np.float64))

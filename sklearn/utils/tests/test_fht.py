import numpy as np
import numpy.testing as npt
import nose.tools as nt
from scipy.linalg import hadamard
from sklearn.utils.cyfht import fht as cyfht
from sklearn.utils.cyfht import fht2 as cyfht2


def test_wikipedia_example():
    input_ = np.array([1, 0, 1, 0, 0, 1, 1, 0], dtype=np.float64)
    copy = input_.copy()
    H = hadamard(8)
    cyfht(input_)
    npt.assert_array_equal(np.dot(copy, H), input_)


def test_numerical_fuzzing_fht():
    for length in [2, 4, 8, 16, 32, 64]:
        input_ = np.random.normal(size=length)
        copy = input_.copy()
        H = hadamard(length)
        cyfht(input_)
        npt.assert_array_almost_equal(np.dot(copy, H), input_)


def test_numerical_fuzzing_fht2():
    for length in [2, 4, 8, 16, 32, 64]:
        for rows in [1, 2, 3, 4, 5]:
            input_ = np.random.normal(size=(rows, length))
            copy = input_.copy()
            H = hadamard(length)
            cyfht2(input_)
            npt.assert_array_almost_equal(np.dot(copy, H), input_)


def test_exception_when_input_not_power_two():
    nt.assert_raises(ValueError, cyfht, np.zeros(9, dtype=np.float64))
    nt.assert_raises(ValueError, cyfht2, np.zeros((2, 9), dtype=np.float64))

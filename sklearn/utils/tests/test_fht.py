import numpy as np
import numpy.testing as npt
from scipy.linalg import hadamard
from sklearn.utils.fht import fht

def test():
    input_ = np.array([1, 0, 1, 0, 0, 1, 1, 0])
    copy = input_.copy()
    H = hadamard(8)
    fht(input_)
    npt.assert_array_equal(np.dot(copy, H), input_)

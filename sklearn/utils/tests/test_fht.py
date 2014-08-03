import numpy as np
import numpy.testing as npt
from scipy.linalg import hadamard
from sklearn.utils.fht import fht
from sklearn.utils.cyfht import fht as cyfht

def single(fht_type):
    input_ = np.array([1, 0, 1, 0, 0, 1, 1, 0], dtype=np.float64)
    copy = input_.copy()
    H = hadamard(8)
    fht_type(input_)
    npt.assert_array_equal(np.dot(copy, H), input_)

def test_all():
    single(fht)
    single(cyfht)

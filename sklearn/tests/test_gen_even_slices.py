""" test gen_even_slices"""

from sklearn.utils import gen_even_slices
from scipy.sparse import csr_matrix
import numpy as np


def check_gen_even_slices():

    even = csr_matrix((1032,1030))
    odd = csr_matrix((1033,1033))
    batch_size = 100

    even_batches = int(np.ceil(float(even.shape[0]) / batch_size))
    odd_batches = int(np.ceil(float(odd.shape[0]) / batch_size))

    odd_slices = list(gen_even_slices(odd_batches * batch_size,
                                            odd_batches, odd.shape[0]))

    even_slices = list(gen_even_slices(even_batches * batch_size,
                                            even_batches, even.shape[0]))

    assert slices_bounds_check(even,even_slices) == "passes", "Fails on Even number of rows"
    assert slices_bounds_check(odd,odd_slices) == "passes", "Fails on Odd number of rows" 

    print("OK")




def slices_bounds_check(matrix, slices):
    try:
        for batch_slice in slices:
            _= matrix[batch_slice]
        return "passes"

    except IndexError:
        return "oob"







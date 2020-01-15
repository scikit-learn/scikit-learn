""" Test simple_heap
"""
import numpy as np

from sklearn.utils._simple_heap import SimpleIntValueHeap


def test_heap():
    rng = np.random.RandomState(0)
    n = 100
    values = rng.rand(n)
    sorted_values = sorted(values)
    res = []
    heap = SimpleIntValueHeap(n)
    for i in range(n):
        heap.update(i, sorted_values[i])

    for i in range(n):
        v = heap.pop()
        res.append(v)

    for i in range(n):
        assert (res[i] == (sorted_values[i], i))

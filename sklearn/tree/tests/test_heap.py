import random
from heapq import heappop, heappush

import pytest

from sklearn.tree._utils import WeightedHeap


@pytest.mark.parametrize("min_heap", [True, False])
def test_weighted_heap(min_heap):
    n = 200
    w_heap = WeightedHeap(n, min_heap=min_heap)
    py_heap = []

    def pop_from_heaps_and_compare():
        top, top_w = w_heap._py_pop()
        top_, top_w_ = heappop(py_heap)
        if not min_heap:
            top_ = -top_
        assert top == top_
        assert top_w == top_w_

    for _ in range(n):
        if len(py_heap) > 0 and random.random() < 1 / 3:
            pop_from_heaps_and_compare()
        else:
            y = random.random()
            w = random.random()
            heappush(py_heap, (y if min_heap else -y, w))
            w_heap._py_push(y, w)

    for _ in range(len(py_heap)):
        pop_from_heaps_and_compare()

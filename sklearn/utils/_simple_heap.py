"""
Simple heap for storing pairs (integer, value) 
and allow decreasing value for given integer.
"""
import numpy as np
import heapq

class SimpleIntValueHeap(object):
    def __init__(self, n):
       self.data = [(np.inf, ordering_idx) for ordering_idx in range(n)]
       self.processed = np.zeros(n, dtype=bool)
    def update(self, idx, value):
        heapq.heappush(self.data, (value, idx))
    def pop(self):
        val, point = heapq.heappop(self.data)
        while self.processed[point]:
            val, point = heapq.heappop(self.data)
        self.processed[point] = True
        return val, point


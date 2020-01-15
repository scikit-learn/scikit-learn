"""
Binary heap implementation for storing pairs (integer, value) 
and allows decreasing value for given integer.
"""
import numpy as np

class UpdateableBinaryHeap(object):
    def __init__(self, n):
       self.data = [(np.inf, ordering_idx) for ordering_idx in range(n)]
       self.tab = [ordering_idx for ordering_idx in range(n)]
       self.processed = np.zeros(n, dtype=bool)
    
    def _siftup(self, pos):
        endpos = len(self.data)
        startpos = pos
        newitem = self.data[pos]
        # Bubble up the smaller child until hitting a leaf.
        childpos = 2*pos + 1    # leftmost child position
        while childpos < endpos:
            # Set childpos to index of smaller child.
            rightpos = childpos + 1
            if rightpos < endpos and self.data[childpos] > self.data[rightpos]:
                childpos = rightpos
            # Move the smaller child up.
            self.tab[self.data[childpos][1]] = pos
            self.data[pos] = self.data[childpos]
            pos = childpos
            childpos = 2*pos + 1
        # The leaf at pos is empty now.  Put newitem there, and bubble it up
        # to its final resting place (by sifting its parents down).
        self.data[pos] = newitem
        self.tab[newitem[1]] = pos
        self._siftdown(startpos, pos)

    def _siftdown(self, startpos, pos):
        newitem = self.data[pos]
        # Follow the path to the root, moving parents down until finding a place
        # newitem fits.
        while pos > startpos:
            parentpos = (pos - 1) >> 1
            parent = self.data[parentpos]
            if newitem < parent:
                self.tab[parent[1]] = pos
                self.data[pos] = parent
                pos = parentpos
                continue
            break
        self.data[pos] = newitem
        self.tab[newitem[1]] = pos
                
    def heappop(self):
        """Pop the smallest item off the heap, maintaining the heap invariant."""
        lastelt = self.data.pop()    # raises appropriate IndexError if heap is empty
        if self.data:
            returnitem = self.data[0]
            self.data[0] = lastelt
            self._siftup(0)
        else:
            returnitem = lastelt
        return returnitem

    def pop(self):
        val, point = self.heappop()
        return val, point
    
    def update(self, idx, value):
        val = (value, idx)
        self.data[self.tab[idx]] = val
        self._siftdown(0,self.tab[idx])




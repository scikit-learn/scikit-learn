import numpy as np
cimport numpy as np

def fht(array_):
    bit = length = len(array_)
    for _ in xrange(int(np.log2(length))):
        bit >>= 1
        for i in xrange(length):
            if i & bit == 0:
                j = i | bit
                temp = array_[i]
                array_[i] += array_[j]
                array_[j] = temp - array_[j]


""" Test fast_dict. So far only smoke test
"""
import numpy as np

from fast_dict import IntFloatDict#, max_merge

def test_int_float_dict():
    keys = np.random.randint(100, size=10).astype(np.int32)
    values = np.random.random(size=10)

    d = IntFloatDict(keys, values)

    #other_keys = np.arange(50).astype(np.int32)[::2]
    #other_values = 0.5*np.ones(50)[::2]
    #other = IntFloatDict(other_keys, other_values)
    #c = max_merge(d, other, mask=np.ones(100, dtype=np.int32))

    d.append(120, 3.)
    [d.append(i+1000, 3.0) for i in xrange(2000)]


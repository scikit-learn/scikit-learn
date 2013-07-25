""" Test fast_dict. So far only smoke test
"""
import numpy as np

from sklearn.utils.fast_dict import IntFloatDict#, max_merge

def test_int_float_dict():
    rng = np.random.RandomState(0)
    keys = np.unique(rng.randint(100, size=10).astype(np.int32))
    values = rng.rand(len(keys))

    d = IntFloatDict(keys, values)
    for key, value in zip(keys, values):
        assert d[key] == value

    #other_keys = np.arange(50).astype(np.int32)[::2]
    #other_values = 0.5*np.ones(50)[::2]
    #other = IntFloatDict(other_keys, other_values)
    #c = max_merge(d, other, mask=np.ones(100, dtype=np.int32))

    d.append(120, 3.)
    assert d[120] == 3.0
    [d.append(i+1000, 4.0) for i in xrange(2000)]
    assert d[1100] == 4.0


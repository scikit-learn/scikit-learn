""" Test fast_dict. So far only smoke test
"""
import numpy as np

from fast_dict import IntFloatDict, MinMergeDict

keys = np.random.randint(100, size=10).astype(np.int32)
values = np.random.random(size=10)

other_keys = np.arange(50).astype(np.int32)
other_values = .5*np.ones(50)

d = IntFloatDict(keys, values)
other = MinMergeDict(other_keys, other_values)

other.update_min(d)


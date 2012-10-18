""" Test fast_dict. So far only smoke test
"""
import numpy as np

from fast_dict import IntFloatDict, min_merge

keys = np.random.randint(100, size=10).astype(np.int32)
values = np.random.random(size=10)

#other_keys = np.arange(50).astype(np.int32)
#other_values = .5*np.ones(50)

d = IntFloatDict(keys, values)
#other = IntFloatDict(other_keys, other_values)
#c = min_merge(d, other, mask=np.ones(100, dtype=np.int32))
#d.append(120, 3.)
[d.append(i+1000, 3.0) for i in xrange(2000)]

import numpy as np
from scipy.stats import rankdata as _sp_rankdata

# To remove when we support scipy 0.13
try:
    _sp_rankdata([1.], 'max')
    rankdata = _sp_rankdata
except TypeError as e:
    def rankdata(a, method="average"):
        if method != "max":
            raise NotImplementedError()

        unique_all, inverse = np.unique(a, return_inverse=True)
        count = np.bincount(inverse, minlength=unique_all.size)
        cum_count = count.cumsum()
        rank = cum_count[inverse]
        return rank


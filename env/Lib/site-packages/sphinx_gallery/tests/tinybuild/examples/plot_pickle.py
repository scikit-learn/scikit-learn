"""
Pickling
--------

This example pickles a function.
"""

import pickle
from math import sqrt

from joblib import Parallel, delayed

assert __name__ == "__main__"
assert "__file__" not in globals()


def function(x):
    """Square root function."""
    return sqrt(x)


pickle.loads(pickle.dumps(function))

# Now with joblib
print(Parallel(n_jobs=2)(delayed(sqrt)(i**2) for i in range(10)))
print(Parallel(n_jobs=2)(delayed(function)(i**2) for i in range(10)))

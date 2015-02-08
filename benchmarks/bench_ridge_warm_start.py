import numpy as np
import time

from sklearn.datasets import fetch_20newsgroups_vectorized
from sklearn.linear_model import Ridge
from sklearn.random_projection import SparseRandomProjection

# Load News20 dataset from scikit-learn.
bunch = fetch_20newsgroups_vectorized(subset="train")
X = bunch.data
y = bunch.target
y[y >= 1] = 1
X = SparseRandomProjection(n_components=1000).fit_transform(X)

alphas = np.logspace(-2, -1, 10)[::-1]

for warm_start in (False, True):
    start = time.clock()
    ridge = Ridge(warm_start=warm_start, fit_intercept=False, solver="sparse_cg")
    for alpha in alphas:
        ridge.set_params(alpha=alpha)
        ridge.fit(X, y)

    print "Warm start", warm_start
    print time.clock() - start
    print

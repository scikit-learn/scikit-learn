"""
Benchmark scikit-learn's Ward implement compared to SciPy's
"""

import time

import numpy as np
from scipy.cluster import hierarchy
import pylab as pl

from sklearn.cluster import AgglomerativeClustering

ward = AgglomerativeClustering(n_clusters=3, linkage='ward')

n_samples = np.logspace(.5, 3, 9)
n_features = np.logspace(1, 3.5, 7)
N_samples, N_features = np.meshgrid(n_samples,
                                    n_features)
scikits_time = np.zeros(N_samples.shape)
scipy_time = np.zeros(N_samples.shape)

for i, n in enumerate(n_samples):
    for j, p in enumerate(n_features):
        X = np.random.normal(size=(n, p))
        t0 = time.time()
        ward.fit(X)
        scikits_time[j, i] = time.time() - t0
        t0 = time.time()
        hierarchy.ward(X)
        scipy_time[j, i] = time.time() - t0

ratio = scikits_time / scipy_time

pl.figure("scikit-learn Ward's method benchmark results")
pl.imshow(np.log(ratio), aspect='auto', origin="lower")
pl.colorbar()
pl.contour(ratio, levels=[1, ], colors='k')
pl.yticks(range(len(n_features)), n_features.astype(np.int))
pl.ylabel('N features')
pl.xticks(range(len(n_samples)), n_samples.astype(np.int))
pl.xlabel('N samples')
pl.title("Scikit's time, in units of scipy time (log)")
pl.show()

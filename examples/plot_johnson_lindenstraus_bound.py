"""
=====================================================================
The Johnson-Lindenstrauss bound for embedding with random projections
=====================================================================


The Johnson-Lindenstrauss states that any high dimensional dataset can
be randomly projected into a lower dimensional Euclidean space while
controlling the distortion in the pairwise distances.

The distortion introduced by a random projection `p` is asserted by
the fact that `p` is defining an eps-embedding with good probability
as defined by:

  (1 - eps) ||u - v||^2 < ||p(u) - p(v)||^2 < (1 + eps) ||u - v||^2

Where u and v are any rows taken from a dataset of shape [n_samples,
n_features] and p is a projection by a random gaussian N(0, 1) matrix
with shape [n_components, n_features] (or a sparse Achlioptas matrix).

The minimum number of components to guarantees the eps-embedding is
given by:

  n_components >= 4 log(n_samples) / (eps^2 / 2 - eps^3 / 3)

"""
import numpy as np
import pylab as pl
from sklearn.random_projection import johnson_lindenstrauss_bound

# range of admissible distortions
eps_range = np.linspace(0.1, 1.0, 5)
colors = pl.cm.Blues(np.linspace(0.3, 1.0, len(eps_range)))

# range of number of samples (observation) to embed
n_samples_range = np.logspace(1, 9, 9)

for eps, color in zip(eps_range, colors):
    min_n_components = johnson_lindenstrauss_bound(n_samples_range, eps=eps)
    pl.loglog(n_samples_range, min_n_components, color=color)

pl.legend(["eps = %0.1f" % eps for eps in eps_range], loc="lower right")
pl.xlabel("Number of observations to eps-embed")
pl.ylabel("Minimum number of dimensions")
pl.title("Johnson-Lindenstrauss bounds")
pl.show()


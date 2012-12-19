"""
=====================================================================
The Johnson-Lindenstrauss bound for embedding with random projections
=====================================================================


The `Johnson-Lindenstrauss lemma`_ states that any high dimensional
dataset can be randomly projected into a lower dimensional Euclidean
space while controlling the distortion in the pairwise distances.

.. _`Johnson-Lindenstrauss lemma`: http://en.wikipedia.org/wiki/Johnson%E2%80%93Lindenstrauss_lemma


Theoretical bounds
==================

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


The first two plots gives a visualization of the minimum number
dimensions ``n_components`` on the number of samples ``n_samples``
to embed through random projections various values of the admissible
distorion eps according to the lemma.


Empirical validation
====================

We validate those bounds on the 20 newsgroups text document (TF-IDF word
frequences) dataset: some 500 documents with 100k features in total are
projected using a sparse random matrix to smaller euclidean spaces with
various values for the target number of dimensions n_components.

For each value of n_components we plot:

- 2D distribution of sample pairs with pairwise distances in original
  and projected spaces as x and y axis respectively.

- 1D histogram of the ratio of those distances (projected / original).

We can see that for low values of n_components the distribution is wide
with many distorted paired and a skewed distribution (due to the hard
limit of zero ratio on the left as distances are always positives)
while for larger values of n_components the distortion is controlled
and the distances are well preserved by the random projection.

"""
from time import time
import numpy as np
import pylab as pl
from sklearn.random_projection import johnson_lindenstrauss_min_dim
from sklearn.random_projection import BernoulliRandomProjection
from sklearn.datasets import fetch_20newsgroups_vectorized
from sklearn.metrics.pairwise import euclidean_distances

# Part 1: plot the theoretical dependency between n_components_min and
# n_samples

# range of admissible distortions
eps_range = np.linspace(0.1, 1.0, 5)
colors = pl.cm.Blues(np.linspace(0.3, 1.0, len(eps_range)))

# range of number of samples (observation) to embed
n_samples_range = np.logspace(1, 9, 9)

pl.figure()
for eps, color in zip(eps_range, colors):
    min_n_components = johnson_lindenstrauss_min_dim(n_samples_range, eps=eps)
    pl.loglog(n_samples_range, min_n_components, color=color)

pl.legend(["eps = %0.1f" % eps for eps in eps_range], loc="lower right")
pl.xlabel("Number of observations to eps-embed")
pl.ylabel("Minimum number of dimensions")
pl.title("Johnson-Lindenstrauss bounds:\nn_samples vs n_components")


# range of admissible distortions
eps_range = np.linspace(0.01, 1.0, 100)

# range of number of samples (observation) to embed
n_samples_range = np.logspace(2, 6, 5)
colors = pl.cm.Blues(np.linspace(0.3, 1.0, len(n_samples_range)))

pl.figure()
for n_samples, color in zip(n_samples_range, colors):
    min_n_components = johnson_lindenstrauss_min_dim(n_samples, eps=eps_range)
    pl.semilogy(eps_range, min_n_components, color=color)

pl.legend(["n_samples = %d" % n for n in n_samples_range], loc="upper right")
pl.xlabel("Distortion eps")
pl.ylabel("Minimum number of dimensions")
pl.title("Johnson-Lindenstrauss bounds:\nn_components vs eps")


# Part 2: perform sparse random projection of some documents of the 20
# newsgroups data which is both high dimensional and sparse

data = fetch_20newsgroups_vectorized().data[:500]
n_samples, n_features = data.shape
print "Embedding %d faces with dim %d using various random projections" % (
    n_samples, n_features)

materialize = True  # switch to False to never allocate the random matrix
n_components_range = np.array([300, 1000, 10000])
dists = euclidean_distances(data, squared=True).ravel()

# select only non-identical samples pairs
nonzero = dists != 0
dists = dists[nonzero]

for n_components in n_components_range:
    t0 = time()
    rp = BernoulliRandomProjection(n_components=n_components)
    projected_data = rp.fit_transform(data)
    print "Projected %d samples from %d to %d in %0.3fs" % (
        n_samples, n_features, n_components, time() - t0)
    if hasattr(rp, 'components_'):
        n_bytes = rp.components_.data.nbytes
        n_bytes += rp.components_.indices.nbytes
        print "Random matrix with size: %0.3fMB" % (
            n_bytes / 1e6)

    projected_dists = euclidean_distances(
        projected_data, squared=True).ravel()[nonzero]

    pl.figure()
    pl.hexbin(dists, projected_dists, gridsize=100)
    pl.xlabel("Pairwise squared distances in original space")
    pl.ylabel("Pairwise squared distances in projected space")
    pl.title("Pairwise distances distribution for n_components=%d" %
             n_components)
    cb = pl.colorbar()
    cb.set_label('Sample pairs counts')

    rates = projected_dists / dists
    print "Mean distances rate: %0.2f (%0.2f)" % (
        np.mean(rates), np.std(rates))

    pl.figure()
    pl.hist(rates, bins=50, normed=True, range=(0., 2.))
    pl.xlabel("Squared distances rate: projected / original")
    pl.ylabel("Distribution of samples pairs")
    pl.title("Histogram of pairwise distance rates for n_components=%d" %
             n_components)

    # TODO: compute the expected value of eps and add them to the previous plot
    # as vertical lines / region

pl.show()

"""Power Iteration Clustering

Scalable alternative to Spectral Clustering for small number of centers.
"""

# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD

import os
import numpy as np

from .k_means_ import k_means
from ..utils.extmath import safe_sparse_dot
from ..utils import inplace_row_normalize
from ..utils import check_random_state


def make_plot(title):
    """Build a plot instance suitable for saving on the filesystem"""
    from pylab import Figure
    plot = Figure(figsize=(7, 5)).add_subplot(111)
    plot.grid(True)
    plot.set_title(title)
    return plot


def save_plot(plot, filename):
    """Save plot as a png file"""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    dirname = os.path.dirname(filename)
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)
    FigureCanvasAgg(plot.get_figure()).print_figure(filename, dpi=80)


def power_iteration_clustering(affinity, k=8, n_vectors=1, tol=1e-5,
                               random_state=0, max_iter=1000, verbose=False,
                               plot_vector=False):
    """Power Iteration Clustering: simple variant of spectral clustering

    One or more random vectors are multiplied several times to the
    row normalized affinity matrix so as to reach a local convergence
    (early stopping before reaching the convergence to the first
    eigen-vector).

    This process imprints the features of the major eigen-vectors into
    the vectors to make them suitable as clustering features given as
    input to the K-Means clustering algorithm.

    This method is supposed to be able to scale to large numbers of
    samples.

    For small to medium problems is it recommended to test the Spectral
    Clustering method first.

    Parameters
    ----------
    k: int
        Number of clusters to find

    n_vectors: int, optional, default to 1
        Number of random vectors to use during power iterations: increase
        value if the number of clusters is large w.r.t the number of
        samples.

    tol: float, optional, default is 1e-5 as recommended in reference paper
        The 'acceleration' convergence criterion (see Reference for details).

    max_iter: int, optional, default is 1000
        Stops after max_iter even if the convergence criterion is not met.

    random_state: a RandomState instance or an int seed (default is 0)
        Pseudo Random Number Generator used to initialize the random vectors
        and the K-Means algorithm.

    verbose: boolean, optional, default is False
        Print convergence info to stdout if True

    plot_vector: boolean, optional, false by default
        Plot the first random vector to files in a 'debug' folder (for
        debugging the convergence only: samples from the same ground
        truth clusters need to be contiguous in the affinity matrix to
        make sense of this).

    Returns
    --------
    labels: array of integer, shape: (n_samples, k)
        The cluster label assignement for each sample.

    Reference
    ---------

    W. Cohen, F. Lin, Power Iteration Clustering, ICML 2010
    http://www.cs.cmu.edu/~wcohen/postscript/icml2010-pic-final.pdf

    Complexity
    ----------

    TODO: this method is supposed to scale better to large n_samples than
    spectral clustering: this remains to be checked in practice

    """
    random_state = check_random_state(random_state)

    if not hasattr(affinity, 'todense'):
        # this is not a sparse matrix: check that this is an array like
        affinity = np.asanyarray(affinity)

    normalized = inplace_row_normalize(affinity.copy())

    n_samples = affinity.shape[0]

    if n_vectors == 1:
        # initialize a single vector deterministically
        sums = affinity.sum(axis=1)
        if hasattr(sums, 'A'):
            sums = sums.A.flatten()
        volume = sums.sum()
        vectors = (sums / volume).reshape((n_vectors, n_samples))
    else:
        # random init
        vectors = random_state.normal(size=(n_vectors, n_samples))

    previous_vectors = vectors.copy()
    delta = np.ones(vectors.size).reshape(vectors.shape)

    for i in range(max_iter):

        previous_vectors[:] = vectors
        previous_delta = delta

        vectors[:] = safe_sparse_dot(vectors, normalized)
        vectors /= np.abs(vectors).sum(axis=1)[:, np.newaxis]

        delta = np.abs(previous_vectors - vectors)
        stopping_gap = np.abs(previous_delta - delta).max() * n_samples

        if verbose and i % 10 == 0:
            print "Power Iteration %04d/%04d: gap=%f" % (
                i + 1, max_iter, stopping_gap)

        if plot_vector and i % 10 ==0:
            p = make_plot("First vector %04d" % (i + 1))
            p.plot(vectors[0])
            save_plot(p, "debug/power_iteration_%04d.png" % (i + 1))

        if stopping_gap < tol:
            break

    if verbose:
        print "Converged at iteration: %04d/%04d with delta=%f" % (
            i + 1, max_iter, delta.max())

    return k_means(vectors.T, k, verbose=verbose,
                   random_state=random_state)[1]

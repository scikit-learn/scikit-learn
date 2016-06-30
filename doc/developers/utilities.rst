.. _developers-utils:

========================
Utilities for Developers
========================

Scikit-learn contains a number of utilities to help with development.  These are
located in :mod:`sklearn.utils`, and include tools in a number of categories.
All the following functions and classes are in the module :mod:`sklearn.utils`.

.. warning ::

   These utilities are meant to be used internally within the scikit-learn
   package.  They are not guaranteed to be stable between versions of
   scikit-learn.  Backports, in particular, will be removed as the scikit-learn
   dependencies evolve.


.. currentmodule:: sklearn.utils

Validation Tools
================

These are tools used to check and validate input.  When you write a function
which accepts arrays, matrices, or sparse matrices as arguments, the following
should be used when applicable.

- :func:`assert_all_finite`: Throw an error if array contains NaNs or Infs.

- :func:`as_float_array`: convert input to an array of floats.  If a sparse
  matrix is passed, a sparse matrix will be returned.

- :func:`check_array`: convert input to 2d array, raise error on sparse
  matrices.  Allowed sparse matrix formats can be given optionally, as well as
  allowing 1d or nd arrays.  Calls :func:`assert_all_finite` by default.

- :func:`check_X_y`: check that X and y have consistent length, calls
  check_array on X, and column_or_1d on y. For multilabel classification or
  multitarget regression, specify multi_output=True, in which case check_array
  will be called on y.

- :func:`indexable`: check that all input arrays have consistent length and can
  be sliced or indexed using safe_index.  This is used to validate input for
  cross-validation.

If your code relies on a random number generator, it should never use
functions like ``numpy.random.random`` or ``numpy.random.normal``.  This
approach can lead to repeatability issues in unit tests.  Instead, a
``numpy.random.RandomState`` object should be used, which is built from
a ``random_state`` argument passed to the class or function.  The function
:func:`check_random_state`, below, can then be used to create a random
number generator object.

- :func:`check_random_state`: create a ``np.random.RandomState`` object from
  a parameter ``random_state``.

  - If ``random_state`` is ``None`` or ``np.random``, then a
    randomly-initialized ``RandomState`` object is returned.
  - If ``random_state`` is an integer, then it is used to seed a new
    ``RandomState`` object.
  - If ``random_state`` is a ``RandomState`` object, then it is passed through.

For example::

    >>> from sklearn.utils import check_random_state
    >>> random_state = 0
    >>> random_state = check_random_state(random_state)
    >>> random_state.rand(4)
    array([ 0.5488135 ,  0.71518937,  0.60276338,  0.54488318])


Efficient Linear Algebra & Array Operations
===========================================

- :func:`extmath.randomized_range_finder`: construct an orthonormal matrix
  whose range approximates the range of the input.  This is used in
  :func:`extmath.randomized_svd`, below.

- :func:`extmath.randomized_svd`: compute the k-truncated randomized SVD.
  This algorithm finds the exact truncated singular values decomposition
  using randomization to speed up the computations. It is particularly
  fast on large matrices on which you wish to extract only a small
  number of components.

- :func:`arrayfuncs.cholesky_delete`:
  (used in :func:`sklearn.linear_model.least_angle.lars_path`)  Remove an
  item from a cholesky factorization.

- :func:`arrayfuncs.min_pos`: (used in ``sklearn.linear_model.least_angle``)
  Find the minimum of the positive values within an array.

- :func:`extmath.norm`: computes Euclidean (L2) vector norm
  by directly calling the BLAS
  ``nrm2`` function.  This is more stable than ``scipy.linalg.norm``.  See
  `Fabian's blog post
  <http://fa.bianp.net/blog/2011/computing-the-vector-norm>`_ for a discussion.

- :func:`extmath.fast_logdet`: efficiently compute the log of the determinant
  of a matrix.

- :func:`extmath.density`: efficiently compute the density of a sparse vector

- :func:`extmath.safe_sparse_dot`: dot product which will correctly handle
  ``scipy.sparse`` inputs.  If the inputs are dense, it is equivalent to
  ``numpy.dot``.

- :func:`extmath.logsumexp`: compute the sum of X assuming X is in the log
  domain. This is equivalent to calling ``np.log(np.sum(np.exp(X)))``, but is
  robust to overflow/underflow errors.  Note that there is similar
  functionality in ``np.logaddexp.reduce``, but because of the pairwise nature
  of this routine, it is slower for large arrays.
  Scipy has a similar routine in ``scipy.misc.logsumexp`` (In scipy versions
  < 0.10, this is found in ``scipy.maxentropy.logsumexp``),
  but the scipy version does not accept an ``axis`` keyword.

- :func:`extmath.weighted_mode`: an extension of ``scipy.stats.mode`` which
  allows each item to have a real-valued weight.

- :func:`resample`: Resample arrays or sparse matrices in a consistent way.
  used in :func:`shuffle`, below.

- :func:`shuffle`: Shuffle arrays or sparse matrices in a consistent way.
  Used in ``sklearn.cluster.k_means``.


Efficient Random Sampling
=========================

- :func:`random.sample_without_replacement`: implements efficient algorithms
  for sampling ``n_samples`` integers from a population of size ``n_population``
  without replacement.


Efficient Routines for Sparse Matrices
======================================

The ``sklearn.utils.sparsefuncs`` cython module hosts compiled extensions to
efficiently process ``scipy.sparse`` data.

- :func:`sparsefuncs.mean_variance_axis`: compute the means and
  variances along a specified axis of a CSR matrix.
  Used for normalizing the tolerance stopping criterion in
  :class:`sklearn.cluster.k_means_.KMeans`.

- :func:`sparsefuncs.inplace_csr_row_normalize_l1` and
  :func:`sparsefuncs.inplace_csr_row_normalize_l2`: can be used to normalize
  individual sparse samples to unit L1 or L2 norm as done in
  :class:`sklearn.preprocessing.Normalizer`.

- :func:`sparsefuncs.inplace_csr_column_scale`: can be used to multiply the
  columns of a CSR matrix by a constant scale (one scale per column).
  Used for scaling features to unit standard deviation in
  :class:`sklearn.preprocessing.StandardScaler`.


Graph Routines
==============

- :func:`graph.single_source_shortest_path_length`:
  (not currently used in scikit-learn)
  Return the shortest path from a single source
  to all connected nodes on a graph.  Code is adapted from `networkx
  <https://networkx.github.io/>`_.
  If this is ever needed again, it would be far faster to use a single
  iteration of Dijkstra's algorithm from ``graph_shortest_path``.

- :func:`graph.graph_laplacian`:
  (used in :func:`sklearn.cluster.spectral.spectral_embedding`)
  Return the Laplacian of a given graph.  There is specialized code for
  both dense and sparse connectivity matrices.

- :func:`graph_shortest_path.graph_shortest_path`:
  (used in :class:`sklearn.manifold.Isomap`)
  Return the shortest path between all pairs of connected points on a directed
  or undirected graph.  Both the Floyd-Warshall algorithm and Dijkstra's
  algorithm are available.  The algorithm is most efficient when the
  connectivity matrix is a ``scipy.sparse.csr_matrix``.


Backports
=========

- :func:`fixes.expit`: Logistic sigmoid function. Replacement for SciPy 0.10's
  ``scipy.special.expit``.

- :func:`sparsetools.connected_components`
  (backported from ``scipy.sparse.connected_components`` in scipy 0.12).
  Used in ``sklearn.cluster.hierarchical``, as well as in tests for
  :mod:`sklearn.feature_extraction`.

- :func:`fixes.isclose`
  (backported from ``numpy.isclose`` in numpy 1.8.1).
  In versions before 1.7, this function was not available in
  numpy. Used in ``sklearn.metrics``.


ARPACK
------

- :func:`arpack.eigs`
  (backported from ``scipy.sparse.linalg.eigs`` in scipy 0.10)
  Sparse non-symmetric eigenvalue decomposition using the Arnoldi
  method.  A limited version of ``eigs`` is available in earlier
  scipy versions.

- :func:`arpack.eigsh`
  (backported from ``scipy.sparse.linalg.eigsh`` in scipy 0.10)
  Sparse non-symmetric eigenvalue decomposition using the Arnoldi
  method.  A limited version of ``eigsh`` is available in earlier
  scipy versions.

- :func:`arpack.svds`
  (backported from ``scipy.sparse.linalg.svds`` in scipy 0.10)
  Sparse non-symmetric eigenvalue decomposition using the Arnoldi
  method.  A limited version of ``svds`` is available in earlier
  scipy versions.


Benchmarking
------------

- :func:`bench.total_seconds` (back-ported from ``timedelta.total_seconds``
  in Python 2.7).  Used in ``benchmarks/bench_glm.py``.


Testing Functions
=================

- :func:`testing.assert_in`, :func:`testing.assert_not_in`: Assertions for
  container membership. Designed for forward compatibility with Nose 1.0.

- :func:`testing.assert_raise_message`: Assertions for checking the
  error raise message.

- :func:`testing.mock_mldata_urlopen`: Mocks the urlopen function to fake
  requests to mldata.org. Used in tests of :mod:`sklearn.datasets`.

- :func:`testing.all_estimators` : returns a list of all estimators in
  scikit-learn to test for consistent behavior and interfaces.

Multiclass and multilabel utility function
==========================================

- :func:`multiclass.is_multilabel`: Helper function to check if the task
  is a multi-label classification one.

- :func:`multiclass.is_label_indicator_matrix`: Helper function to check if
  a classification output is in label indicator matrix format.

- :func:`multiclass.unique_labels`: Helper function to extract an ordered
  array of unique labels from different formats of target.


Helper Functions
================

- :class:`gen_even_slices`: generator to create ``n``-packs of slices going up
  to ``n``.  Used in ``sklearn.decomposition.dict_learning`` and
  ``sklearn.cluster.k_means``.

- :func:`safe_mask`: Helper function to convert a mask to the format expected
  by the numpy array or scipy sparse matrix on which to use it (sparse
  matrices support integer indices only while numpy arrays support both
  boolean masks and integer indices).

- :func:`safe_sqr`: Helper function for unified squaring (``**2``) of
  array-likes, matrices and sparse matrices.


Hash Functions
==============

- :func:`murmurhash3_32` provides a python wrapper for the
  ``MurmurHash3_x86_32`` C++ non cryptographic hash function. This hash
  function is suitable for implementing lookup tables, Bloom filters,
  Count Min Sketch, feature hashing and implicitly defined sparse
  random projections::

    >>> from sklearn.utils import murmurhash3_32
    >>> murmurhash3_32("some feature", seed=0) == -384616559
    True

    >>> murmurhash3_32("some feature", seed=0, positive=True) == 3910350737
    True

  The ``sklearn.utils.murmurhash`` module can also be "cimported" from
  other cython modules so as to benefit from the high performance of
  MurmurHash while skipping the overhead of the Python interpreter.


Warnings and Exceptions
=======================

- :class:`deprecated`: Decorator to mark a function or class as deprecated.

- :class:`sklearn.exceptions.ConvergenceWarning`: Custom warning to catch
  convergence problems. Used in ``sklearn.covariance.graph_lasso``.

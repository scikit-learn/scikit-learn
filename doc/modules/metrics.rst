.. _metrics:

Pairwise metrics, Affinities and Kernels
========================================

The :mod:`sklearn.metrics.pairwise` submodule implements utilities to evaluate
pairwise distances or affinity of sets of samples.

This module contains both distance metrics and kernels. A brief summary is
given on the two here.

Distance metrics are a function d(a, b) such that d(a, b) < d(a, c) if objects
a and b are considered "more similar" to objects a and c. Two objects exactly
alike would have a distance of zero.
One of the most popular examples is Euclidean distance.
To be a 'true' metric, it must obey the following four conditions::

    1. d(a, b) >= 0, for all a and b
    2. d(a, b) == 0, if and only if a = b, positive definiteness
    3. d(a, b) == d(b, a), symmetry
    4. d(a, c) <= d(a, b) + d(b, c), the triangle inequality

Kernels are measures of similarity, i.e. ``s(a, b) > s(a, c)``
if objects ``a`` and ``b`` are considered "more similar" to objects
``a`` and ``c``. A kernel must also be positive semi-definite.

There are a number of ways to convert between a distance metric and a
similarity measure, such as a kernel. Let D be the distance, and S be the
kernel:

    1. ``S = np.exp(-D * gamma)``, where one heuristic for choosing
       ``gamma`` is ``1 / num_features``
    2. ``S = 1. / (D / np.max(D))``


.. currentmodule:: sklearn.metrics.pairwise

Cosine similarity
-----------------
:func:`cosine_similarity` computes the L2-normalized dot product of vectors.
That is, if :math:`x` and :math:`y` are row vectors,
their cosine similarity :math:`k` is defined as:

.. math::

    k(x, y) = \frac{x \dot y^\top}{\|x\| \|y\|}

This is called cosine similarity, because Euclidean (L2) normalization
projects the vectors onto the unit sphere,
and their dot product is then the cosine of the angle between the points
denoted by the vectors.

This kernel is a popular choice for computing the similarity of documents
represented as tf-idf vectors.
:func:`cosine_similarity` accepts ``scipy.sparse`` matrices.
(Note that the tf-idf functionality in ``sklearn.feature_extraction.text``
can produce normalized vectors, in which case :func:`cosine_similarity`
is equivalent to :func:`linear_kernel`, only slower.)

.. topic:: References:

    * C.D. Manning, P. Raghavan and H. Schütze (2008). Introduction to
      Information Retrieval. Cambridge University Press.
      http://nlp.stanford.edu/IR-book/html/htmledition/the-vector-space-model-for-scoring-1.html

Chi Squared Kernel
------------------
The chi squared kernel is a very popular choice for training non-linear SVMs in
Computer Vision applications.
It can be computed using :func:`chi2_kernel` and then passed to an
:class:`sklearn.svm.SVC` with ``kernel="precomputed"``::

    >>> from sklearn.svm import SVC
    >>> from sklearn.metrics.pairwise import chi2_kernel
    >>> X = [[0, 1], [1, 0], [.2, .8], [.7, .3]]
    >>> y = [0, 1, 0, 1]
    >>> K = chi2_kernel(X, gamma=.5)
    >>> K                        # doctest: +ELLIPSIS
    array([[ 1.        ,  0.36...,  0.89...,  0.58...],
           [ 0.36...,  1.        ,  0.51...,  0.83...],
           [ 0.89...,  0.51...,  1.        ,  0.77... ],
           [ 0.58...,  0.83...,  0.77... ,  1.        ]])

    >>> svm = SVC(kernel='precomputed').fit(K, y)
    >>> svm.predict(K)
    array([0, 1, 0, 1])

It can also be directly used as the ``kernel`` argument::

    >>> svm = SVC(kernel=chi2_kernel).fit(X, y)
    >>> svm.predict(X)
    array([0, 1, 0, 1])


The chi squared kernel is given by

.. math::

        k(x, y) = \exp \left (-\gamma \sum_i \frac{(x[i] - y[i]) ^ 2}{x[i] + y[i]} \right )

The data is assumed to be non-negative, and is often normalized to have an L1-norm of one.
The normalization is rationalized with the connection to the chi squared distance,
which is a distance between discrete probability distributions.

The chi squared kernel is most commonly used on histograms (bags) of visual words.

.. topic:: References:

    * Zhang, J. and Marszalek, M. and Lazebnik, S. and Schmid, C.
      Local features and kernels for classification of texture and object
      categories: A comprehensive study
      International Journal of Computer Vision 2007
      http://eprints.pascal-network.org/archive/00002309/01/Zhang06-IJCV.pdf


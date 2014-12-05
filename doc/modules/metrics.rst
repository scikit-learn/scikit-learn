.. _metrics:

Pairwise metrics, Affinities and Kernels
========================================

The :mod:`sklearn.metrics.pairwise` submodule implements utilities to evaluate
pairwise distances or affinity of sets of samples.

This module contains both distance metrics and kernels. A brief summary is
given on the two here.

Distance metrics are functions ``d(a, b)`` such that ``d(a, b) < d(a, c)``
if objects ``a`` and ``b`` are considered "more similar" than objects ``a``
and ``c``. Two objects exactly alike would have a distance of zero.
One of the most popular examples is Euclidean distance.
To be a 'true' metric, it must obey the following four conditions::

    1. d(a, b) >= 0, for all a and b
    2. d(a, b) == 0, if and only if a = b, positive definiteness
    3. d(a, b) == d(b, a), symmetry
    4. d(a, c) <= d(a, b) + d(b, c), the triangle inequality

Kernels are measures of similarity, i.e. ``s(a, b) > s(a, c)``
if objects ``a`` and ``b`` are considered "more similar" than objects
``a`` and ``c``. A kernel must also be positive semi-definite.

There are a number of ways to convert between a distance metric and a
similarity measure, such as a kernel. Let ``D`` be the distance, and ``S`` be
the kernel:

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

    k(x, y) = \frac{x y^\top}{\|x\| \|y\|}

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

Linear kernel
-------------
The function :func:`linear_kernel` computes the linear kernel, that is, a
special case of :func:`polynomial_kernel` with ``degree=1`` and ``coef0=0`` (homogeneous).
If ``x`` and ``y`` are column vectors, their linear kernel is:

.. math::

    k(x, y) = x^\top y

Polynomial kernel
-----------------
The function :func:`polynomial_kernel` computes the degree-d polynomial kernel
between two vectors. The polynomial kernel represents the similarity between two
vectors. Conceptually, the polynomial kernels considers not only the similarity
between vectors under the same dimension, but also across dimensions. When used
in machine learning algorithms, this allows to account for feature interaction.

The polynomial kernel is defined as:

.. math::

    k(x, y) = (\gamma x^\top y +c_0)^d

where:

    * ``x``, ``y`` are the input vectors
    * ``d`` is the kernel degree

If :math:`c_0 = 0` the kernel is said to be homogeneous.

Sigmoid kernel
--------------
The function :func:`sigmoid_kernel` computes the sigmoid kernel between two
vectors. The sigmoid kernel is also known as hyperbolic tangent, or Multilayer
Perceptron (because, in the neural network field, it is often used as neuron
activation function). It is defined as:

.. math::

    k(x, y) = \tanh( \gamma x^\top y + c_0)

where:

    * ``x``, ``y`` are the input vectors
    * :math:`\gamma` is known as slope
    * :math:`c_0` is known as intercept

RBF kernel
----------
The function :func:`rbf_kernel` computes the radial basis function (RBF) kernel
between two vectors. This kernel is defined as:

.. math::

    k(x, y) = \exp( -\gamma \| x-y \|^2)

where ``x`` and ``y`` are the input vectors. If :math:`\gamma = \sigma^{-2}`
the kernel is known as the Gaussian kernel of variance :math:`\sigma^2`.

Matérn kernel
-------------
The function :func:`matern_kernel` is a generalization of the RBF kernel. It has
an additional parameter :math:`\nu` (set via the keyword coef0) which controls
the smoothness of the resulting function. The general functional form of a
Matérn is given by

.. math::

    k(d) = \sigma^2\frac{1}{\Gamma(\nu)2^{\nu-1}}\Bigg(\gamma\sqrt{2\nu} d\Bigg)^\nu K_\nu\Bigg(\gamma\sqrt{2\nu} d\Bigg),

where :math:`d=\| x-y \|` and ``x`` and ``y`` are the input vectors. 

As :math:`\nu\rightarrow\infty`, the Matérn kernel converges to the RBF kernel.
When :math:`\nu = 1/2`, the Matérn kernel becomes identical to the absolute
exponential kernel, i.e.,

.. math::
    k(d) = \sigma^2 \exp \Bigg(-\gamma d \Bigg) \quad \quad \nu= \tfrac{1}{2}

In particular, :math:`\nu = 3/2`:

.. math::
    k(d) = \sigma^2 \Bigg(1 + \gamma \sqrt{3} d \Bigg) \exp \Bigg(-\gamma \sqrt{3}d \Bigg) \quad \quad \nu= \tfrac{3}{2}

and :math:`\nu = 5/2`:

.. math::
    k(d) = \sigma^2 \Bigg(1 + \gamma \sqrt{5}d +\frac{5}{3} \gamma^2d^2 \Bigg) \exp \Bigg(-\gamma \sqrt{5}d \Bigg) \quad \quad \nu= \tfrac{5}{2}

are popular choices for learning functions that are not infinitely
differentiable (as assumed by the RBF kernel) but at least once (:math:`\nu =
3/2`) or twice differentiable (:math:`\nu = 5/2`).

The following example illustrates how the Matérn kernel's covariance decreases
with increasing dissimilarity of the two inputs for different values of coef0
(the parameter :math:`\nu` of the Matérn kernel):

.. figure:: ../auto_examples/metrics/images/plot_matern_kernel_001.png
    :target: ../auto_examples/metrics/plot_matern_kernel.html
    :align: center

The flexibility of controlling the smoothness of the learned function via coef0
allows adapting to the properties of the true underlying functional relation.
The following example shows that support vector regression with Matérn kernel
with smaller values of coef0 can better approximate a discontinuous 
step-function:

.. figure:: ../auto_examples/svm/images/plot_svm_matern_kernel_001.png
    :target: ../auto_examples/svm/plot_svm_matern_kernel.html
    :align: center

See Rasmussen and Williams 2006, pp84 for further details regarding the
different variants of the Matérn kernel.


Chi-squared kernel
------------------
The chi-squared kernel is a very popular choice for training non-linear SVMs in
computer vision applications.
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

    * Rasmussen, C. E. and Williams, C.
      Gaussian Processes for Machine Learning
      The MIT Press, 2006
      http://www.gaussianprocess.org/gpml/chapters/


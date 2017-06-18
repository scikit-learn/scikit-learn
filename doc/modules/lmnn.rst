.. _lmnn:

=============================
Large Margin Nearest Neighbor
=============================

.. sectionauthor:: John Chiotellis <johnyc.code@gmail.com>

.. currentmodule:: sklearn.neighbors

Large Margin Nearest Neighbor (:class:`lmnn.LargeMarginNearestNeighbor`) is
a distance metric learning algorithm.

This transformer aims to learn a Mahalanobis metric that improves the
accuracy of k-nearest neighbor classification compared to the standard
Euclidean metric.

As in support vector machines (SVMs), this method attempts to maximize the
margin between data samples that belong to different classes, therefore the
name of the method. Instead of a set of hyper-planes in a high or infinite
dimensional space, LMNN constructs a set of hyper-spheres centered at each
training point that are supposed to contain only samples that share the same
label as the central point.

In the beginning of the algorithm each training sample fixes :math:`k`
target neighbors, namely :math:`k` nearest training samples according to some
initial - most commonly the Euclidean - metric that share the same label.
The first goal is to find a metric such that the distances of each sample to
its target neighbors are minimized.

The second goal of the algorithm is to create a large margin between the
farthest target neighbor of each sample and data samples that belong to
different classes. Data samples from different classes that violate this
margin are called impostors.

Combined with a nearest neighbors classifier (:class:`classification.KNeighborsClassifier`)
this method is attractive for classification because it requires no
modification or extension for multi-class (as opposed to binary) problems
and only a single parameter (``n_neighbors``) has to be selected by the user
before training.

Large Margin Nearest Neighbor classification has been proven to work well in
practice for data sets of varying size and difficulty.


.. Enter ../auto_examples/classification/images/sphx_glr_plot_lda_lmnn_001.png
        :target: ../auto_examples/classification/plot_lda_lmnn.html
        :scale: 80

.. Enter centered:: |lmnn|

.. The plot shows decision boundaries for nearest neighbor classification and
    large margin nearest neighbor classification. The bottom row shows the
    same comparison after first doing PCA on the inputs.

.. Enter topic:: Examples:
    .. Comment
    :ref:`sphx_glr_auto_examples_classification_plot_lda_lmnn.py`: Comparison
    of LDA and LMNN on synthetic data.

Dimensionality reduction
========================

:class:`lmnn.LargeMarginNearestNeighbor` can be used to
perform supervised dimensionality reduction, by projecting the input data to a
linear subspace consisting of the directions which maximize the nearest
neighbor classification accuracy and therefore the separation between
clusters of samples of the same class.

This is implemented in :func:`lmnn.LargeMarginNearestNeighbor.transform`.
The desired dimensionality can be set using the ``n_features_out`` parameter
in the constructor.

.. Comment topic:: Examples:
    ..
    :ref:`sphx_glr_auto_examples_decomposition_plot_pca_vs_lda_vs_lmnn.py`:
    Comparison of PCA, LDA and LMNN
    for dimensionality reduction of the Olivetti faces dataset

Mathematical formulation
========================

The LMNN objective function consists of two competing terms, the pull loss
that pulls target neighbors together and the push loss that pushes impostors
 away:

.. math::

    \varepsilon_{\text{pull}} (L) = \sum_{i, j \rightsquigarrow i} ||L
    (x_i - x_j)||^2

    \varepsilon_{\text{push}} (L) = \sum_{i, j \rightsquigarrow i}
    \sum_{l} (1 - y_{il}) [1 + || L(x_i - x_j)||^2 - || L
    (x_i - x_l)||^2]_+
    \text{ where } y_{il} = 1 \text{ if } $y_i = y_l$ \text{ and 0 otherwise}.

    \varepsilon(L) = (1 - \mu) \varepsilon_{\text{pull}} (L) +
    \mu \varepsilon_{\text{push}} (L) \text{, } \quad \mu \in [0,1]


The term inside the sum of the push loss is called the hingle loss and
the notation amounts to :math:`[x]_+ = \max(0, x)`.
The parameter :math:`\mu` is a trade-off between penalizing large distances
between target neighbors and penalizing margin violations. In practice,
the two terms are weighted equally.

To use this model for classification, we just need to use a nearest neighbors
classifier (:class:`classification.KNeighborsClassifier`) on the data
transformed by the linear transformation :math:`L` found by LMNN.

The loss function is not convex in :math:`L`, but there is an alternative
formulation:

.. math::

    \varepsilon_{\text{pull}} (M) = \sum_{i, j \rightsquigarrow i} (x_i-x_j)
    ^T M (x_i-x_j)

    \varepsilon_{\text{push}} (M) = \sum_{i, j \rightsquigarrow i}
    \sum_{l} (1 - y_{il}) [1 + (x_i-x_j)^T M (x_i-x_j) - (x_i-x_l)^T M (x_i-x_l)]_+
    \text{ where } y_{il} = 1 \text{ if } $y_i = y_l$ \text{ and 0 otherwise}.

    \varepsilon(M) = (1 - \mu) \varepsilon_{\text{pull}} (M) +
    \mu \varepsilon_{\text{push}} (M) \text{, } \quad \mu \in [0,1]


This objective is now convex in :math:`M` but the matrix :math:`M = L^T L`
must be constrained to be positive semidefinite as opposed to :math:`L`
which could be any real matrix.


In contrast to related methods such as LDA, LMNN does not make any
assumptions about the class distributions. The nearest neighbor
classification can naturally produce highly irregular decision boundaries.
See [#1]_ for more details.



Implementation
==============

This implementation follows closely the MATLAB implementation found at
https://bitbucket.org/mlcircus/lmnn which solves the unconstrained problem.
It finds a linear transformation :math:`L` by optimisation with
L-BFGS instead of solving the constrained problem that finds the globally
optimal metric. Different from the paper, the problem solved by this
implementation is with the squared hinge loss (to make the problem
differentiable).


.. topic:: References:

   .. [#1] "Distance Metric Learning for Large Margin Nearest Neighbor
            Classification.", Weinberger, Kilian Q., and Lawrence K. Saul, Journal
            of Machine Learning Research, Vol. 10, Feb. 2009, pp. 207-244.

    http://jmlr.csail.mit.edu/papers/volume10/weinberger09a/weinberger09a.pdf

   .. [#2] Wikipedia entry on Large Margin Nearest Neighbor

    https://en.wikipedia.org/wiki/Large_margin_nearest_neighbor

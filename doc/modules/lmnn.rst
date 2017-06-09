.. _lda_qda:

=============================
Large Margin Nearest Neighbor
=============================

.. sectionauthor:: John Chiotellis <johnyc.code@gmail.com>

.. currentmodule:: sklearn.neighbors

Large Margin Nearest Neighbor (:class:`lmnn.LargeMarginNearestNeighbor`) is
a method for distance metric learning.

This transformer aims to learn a Mahalanobis metric that improves the
accuracy of k-nearest neighbor classification which by default uses the
Euclidean metric.
As in support vector machines (SVMs), this method attempts to maximize the
margin between samples that belong to different classes, thus its name.
Combined with a nearest neighbors classifier (:class:`lmnn.KNeighborsClassifier`)
this method is attractive for classification because it requires no
modification or extension for multi-class (as opposed to binary) problems
and only a single parameter (``n_neighbors``) has to be selected by the user
before training.

Large Margin Nearest Neighbor classification has been proven to work well in
practice for data sets of varying size and difficulty.

.. |ldaqda| image:: ../auto_examples/classification/images/sphx_glr_plot_lda_qda_001.png
        :target: ../auto_examples/classification/plot_lda_qda.html
        :scale: 80

.. centered:: |ldaqda|

The plot shows decision boundaries for simple nearest neighbor
classification and large margin nearest neighbor classification.
The bottom row shows the same comparison after first doing PCA on the inputs.

.. topic:: Examples:

    :ref:`sphx_glr_auto_examples_classification_plot_lda_qda.py`: Comparison of LDA and QDA
    on synthetic data.

Dimensionality reduction using Large Margin Nearest Neighbor
============================================================

:class:`lmnn.LargeMarginNearestNeighbor` can be used to
perform supervised dimensionality reduction, by projecting the input data to a
linear subspace consisting of the directions which maximize the separation
between classes (in a precise sense discussed in the mathematics section
below).

This is implemented in
:func:`lmnn.LargeMarginNearestNeighbor.transform`. The desired
dimensionality can be set using the ``n_features_out`` constructor parameter.
This parameter has no influence on :func:`lmnn.LargeMarginNearestNeighbor.fit`.

.. topic:: Examples:

    :ref:`sphx_glr_auto_examples_decomposition_plot_pca_vs_lda.py`:
    Comparison of PCA, LDA and LMNN
    for dimensionality reduction of the Iris dataset

Mathematical formulation of the Large Margin Nearest Neighbor
=============================================================

The LMNN objective function consists of two competing terms, the pull loss
that pulls target neighbors together and the push loss that pushes impostors
 away:

.. math::
    \newcommand{\bL}{\textbf{L}}

    \varepsilon_{\text{pull}} (\bL) = \sum_{i, j \rightsquigarrow i} ||\bL(\vec{x}_i - \vec{x}_j)||^2

    \varepsilon_{\text{push}} (\bL) = \sum_{i, j \rightsquigarrow i}
\sum_{l} (1 - y_{il}) [1 + ||\bL(\vec{x}_i - \vec{x}_j)||^2 - ||\bL(\vec{x}_i - \vec{x}_l)||^2]_+
\text{ where } y_{il} = 1 \text{ if } $y_i = y_l$ \text{ and 0 otherwise}.

    \varepsilon(\bL) = (1 - \mu) \varepsilon_{\text{pull}} (\bL) +
\mu \varepsilon_{\text{push}} (\bL) \text{, } \quad \mu \in [0,1]



and we minimize this loss with L-BFGS.

The second term of the pull loss is the hingle loss, namely :math:`[x]_+ =
\max(0, x)`.


To use this model for classification, we just need to use a nearest neighbors
classifier (:func:`KNeighborsClassifier`) in the input space, transformed by
the linear transformation found by LMNN.


In contrast to related methods such as LDA, LMNN does not need to make any
assumptions about the class distributions since the nearest neighbor
classification naturally produces non-linear decision boundaries. See [#1]_
for more details.


Mathematical formulation of LMNN dimensionality reduction
=========================================================

LMNN can be used for dimensionality reduction by setting the
``n_features_out`` parameter in the constructor. The algorithm will find an
optimal transformation for a subspace of dimension ``n_features_out``.



Implementation
==============

This implementation follows closely the original MATLAB implementation found
at https://bitbucket.org/mlcircus/lmnn which solves the unconstrained
problem. It finds a linear transformation by optimisation with L-BFGS instead
of solving the constrained problem that finds the globally optimal metric.
Different from the paper, the problem solved by this implementation is
with the squared hinge loss (to make the problem differentiable).


.. topic:: Examples:

    :ref:`sphx_glr_auto_examples_classification_plot_lda.py`: Comparison of LDA classifiers
    with and without shrinkage.

.. topic:: References:

   .. [#1] "Distance Metric Learning for Large Margin Nearest Neighbor
      Classification.", Weinberger, Kilian Q., and Lawrence K. Saul, Journal
      of Machine Learning Research, Vol. 10, Feb. 2009, pp. 207-244.
    (http://jmlr.csail.mit.edu/papers/volume10/weinberger09a/weinberger09a.pdf)

   .. [#2] Wikipedia entry on Large Margin Nearest Neighbor
    (https://en.wikipedia.org/wiki/Large_margin_nearest_neighbor)

.. _lda_qda:

==========================================
Linear and quadratic discriminant analysis
==========================================

.. currentmodule:: sklearn

Linear discriminant analysis (:class:`lda.LDA`) and
quadratic discriminant analysis (:class:`qda.QDA`)
are two classic classifiers, with, as their names suggest, a linear and a
quadratic decision surface, respectively.

These classifiers are attractive because they have closed-form solutions that
can be easily computed, are inherently multiclass,
and have proven to work well in practice.
Also there are no parameters to tune for these algorithms.

.. |ldaqda| image:: ../auto_examples/images/plot_lda_qda_1.png
        :target: ../auto_examples/plot_lda_qda.html
        :scale: 80

.. centered:: |ldaqda|

The plot shows decision boundaries for LDA and QDA. The bottom row
demonstrates that LDA can only learn linear boundaries, while QDA can learn
quadratic boundaries and is therefore more flexible.

.. topic:: Examples:

    :ref:`example_plot_lda_qda.py`: Comparison of LDA and QDA on synthetic data.

.. topic:: References:

     .. [3] "The Elements of Statistical Learning", Hastie T., Tibshirani R.,
        Friedman J., 2008.


Dimensionality reduction using LDA
==================================

:class:`lda.LDA` can be used to perform supervised dimensionality reduction by
projecting the input data to a subspace consisting of the most
discriminant directions.
This is implemented in :func:`lda.LDA.transform`. The desired
dimensionality can be set using the `n_components` constructor
parameter. This parameter has no influence on :func:`lda.LDA.fit` or :func:`lda.LDA.predict`.


Mathematical Idea
=================

Both methods work by modeling the class conditional distribution of the data :math:`P(X|y=k)`
for each class `k`. Predictions can be obtained by using Bayes' rule:

.. math::
    P(y | X) = P(X | y) \cdot P(y) / P(X) = P(X | y) \cdot P(Y) / ( \sum_{y'} P(X | y') \cdot p(y'))

In linear and quadratic discriminant analysis, `P(X|y)` is modelled as a Gaussian distribution.
In the case of LDA, the Gaussians for each class are assumed to share the same covariance matrix.
This leads to a linear decision surface, as can be seen by comparing the the log-probability rations
:math:`log[P(y=k | X) / P(y=l | X)]`.

In the case of QDA, there are no assumptions on the covariance matrices of the Gaussians,
leading to a quadratic decision surface.

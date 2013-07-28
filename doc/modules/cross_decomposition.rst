.. _cross_decomposition:

===================
Cross decomposition
===================

.. currentmodule:: sklearn.cross_decomposition

The cross decomposition module contains two main families of algorithms: the
partial least squares (PLS) and the canonical correlation analysis (CCA).

These families of algorithms are useful to find linear relations between two
multivariate datasets: the `X` and `Y` arguments of the `fit` method
are 2D arrays.

.. figure:: ../auto_examples/cross_decomposition/images/plot_compare_cross_decomposition_1.png
   :target: ../auto_examples/cross_decomposition/plot_compare_cross_decomposition.html
   :scale: 75%
   :align: center


Cross decomposition algorithms find the fundamental relations between two
matrices (X and Y). They are latent variable approaches to modeling the
covariance structures in these two spaces. They will try to find the
multidimensional direction in the X space that explains the maximum
multidimensional variance direction in the Y space. PLS-regression is
particularly suited when the matrix of predictors has more variables than
observations, and when there is multicollinearity among X values. By contrast,
standard regression will fail in these cases.

Classes included in this module are :class:`PLSRegression`
:class:`PLSCanonical`, :class:`CCA` and :class:`PLSSVD`


.. topic:: Reference:

   * JA Wegelin
     `A survey of Partial Least Squares (PLS) methods, with emphasis on the two-block case <https://www.stat.washington.edu/www/research/reports/2000/tr371.pdf>`_

.. topic:: Examples:

    * :ref:`example_cross_decomposition_plot_compare_cross_decomposition.py`

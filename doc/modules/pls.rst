.. _pls:

======================
Partial Least Squares
======================

.. currentmodule:: sklearn.pls

Partial least squares (PLS) models are useful to find linear relations
between two multivariate datasets: in PLS the `X` and `Y` arguments of
the `fit` method are 2D arrays.

.. figure:: ../auto_examples/images/plot_pls_1.png
   :target: ../auto_examples/plot_pls.html
   :scale: 75%
   :align: center


PLS finds the fundamental relations between two matrices
(X and Y): it is a latent variable approach to modeling the covariance
structures in these two spaces. A PLS model will try to find the
multidimensional direction in the X space that explains the maximum
multidimensional variance direction in the Y space. PLS-regression is
particularly suited when the matrix of predictors has more variables
than observations, and when there is multicollinearity among X
values. By contrast, standard regression will fail in these cases.

Classes included in this module are :class:`PLSRegression`
:class:`PLSCanonical`, :class:`CCA` and :class:`PLSSVD`


.. topic:: Reference:

   * JA Wegelin
     `A survey of Partial Least Squares (PLS) methods, with emphasis on the two-block case <https://www.stat.washington.edu/www/research/reports/2000/tr371.pdf>`_

.. topic:: Examples:

    * :ref:`example_plot_pls.py`


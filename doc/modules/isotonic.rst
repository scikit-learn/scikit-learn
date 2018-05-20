.. _isotonic:

===================
Isotonic regression
===================

.. currentmodule:: sklearn.isotonic

The class :class:`IsotonicRegression` fits a non-decreasing function to data.
It solves the following problem:

  minimize :math:`\sum_i w_i (y_i - \hat{y}_i)^2`

  subject to :math:`\hat{y}_{min} = \hat{y}_1 \le \hat{y}_2 ... \le \hat{y}_n = \hat{y}_{max}`

where each :math:`w_i` is strictly positive and each :math:`y_i` is an
arbitrary real number. It yields the vector which is composed of non-decreasing
elements the closest in terms of mean squared error. In practice this list
of elements forms a function that is piecewise linear.

.. figure:: ../auto_examples/images/sphx_glr_plot_isotonic_regression_001.png
   :target: ../auto_examples/plot_isotonic_regression.html
   :align: center

.. _isotonic:

===================
Isotonic regression
===================

.. currentmodule:: sklearn.isotonic

The class :class:`IsotonicRegression` fits a non-decreasing real function to
1-dimensional data. It solves the following problem:

  minimize :math:`\sum_i w_i (y_i - \hat{y}_i)^2`

  subject to :math:`\hat{y}_i \le \hat{y}_j` whenever :math:`X_i \le X_j`,

where the weights :math:`w_i` are strictly positive, and both `X` and `y` are
arbitrary real quantities.

The `increasing` parameter changes the constraint to
:math:`\hat{y}_i \ge \hat{y}_j` whenever :math:`X_i \le X_j`. Setting it to
'auto' will automatically choose the constraint based on `Spearman's rank
correlation coefficient
<https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient>`_.

:class:`IsotonicRegression` produces a series of predictions
:math:`\hat{y}_i` for the training data which are the closest to the targets
:math:`y` in terms of mean squared error. These predictions are interpolated
for predicting to unseen data. The predictions of :class:`IsotonicRegression`
thus form a function that is piecewise linear:

.. figure:: ../auto_examples/miscellaneous/images/sphx_glr_plot_isotonic_regression_001.png
   :target: ../auto_examples/miscellaneous/plot_isotonic_regression.html
   :align: center

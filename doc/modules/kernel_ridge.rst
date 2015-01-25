.. _kernel_ridge:

===========================
Kernel ridge regression
===========================

.. currentmodule:: sklearn.kernel_ridge

Kernel ridge regression (KRR) [M2012]_ combines :ref:`ridge_regression`
(linear least squares with l2-norm regularization) with the kernel trick. It
thus learns a linear function in the space induced by the respective kernel and
the data. For non-linear kernels, this corresponds to a non-linear
function in the original space.

The form of the model learned by :class:`KernelRidge` is identical to support
vector regression (:class:`SVR`). However, different loss functions are used:
KRR uses squared error loss while support vector regression uses
:math:`\epsilon`-insensitive loss, both combined with l2 regularization.  In
contrast to :class:`SVR`, fitting :class:`KernelRidge` can be done in
closed-form and is typically faster for medium-sized datasets. On the other
hand, the learned model is non-sparse and thus slower than SVR, which learns
a sparse model for :math:`\epsilon > 0`, at prediction-time.

The following figure compares :class:`KernelRidge` and :class:`SVR` on
an artificial dataset, which consists of a sinusoidal target function and
strong noise added to every fifth datapoint. The learned model of
:class:`KernelRidge` and :class:`SVR` is plotted, where both
complexity/regularization and bandwidth of the RBF kernel have been optimized
using grid-search. The learned functions are very similar; however, fitting
:class:`KernelRidge` is approx. seven times faster than fitting :class:`SVR`
(both with grid-search). However, prediction of 100000 target values is more
than three times faster with SVR since it has learned a sparse model using only
approx. 1/3 of the 100 training datapoints as support vectors.

.. figure:: ../auto_examples/images/plot_kernel_ridge_regression_001.png
   :target: ../auto_examples/plot_kernel_ridge_regression.html
   :align: center

The next figure compares the time for fitting and prediction of
:class:`KernelRidge` and :class:`SVR` for different sizes of the training set.
Fitting :class:`KernelRidge` is faster than :class:`SVR` for medium-sized
training sets (less than 1000 samples); however, for larger training sets
:class:`SVR` scales better. With regard to prediction time, :class:`SVR` is
faster than :class:`KernelRidge` for all sizes of the training set because of
the learned sparse solution. Note that the degree of sparsity and thus the
prediction time depends on the parameters :math:`\epsilon` and :math:`C` of the
:class:`SVR`; :math:`\epsilon = 0` would correspond to a dense model.

.. figure:: ../auto_examples/images/plot_kernel_ridge_regression_002.png
   :target: ../auto_examples/plot_kernel_ridge_regression.html
   :align: center


.. topic:: References:

    .. [M2012] "Machine Learning: A Probabilistic Perspective"
      Murphy, K. P. - chapter 14.4.3, pp. 492-493, The MIT Press, 2012

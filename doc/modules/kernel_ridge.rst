.. _kernel_ridge:

===========================
Kernelized ridge regression
===========================

.. currentmodule:: sklearn.kernel_ridge

Kernelized ridge regression implemented in the class :class:`KernelRidge`
combines ridge regression (linear least squares plus l2-norm  regularization)
with the kernel trick. It thus learns a linear function in the Reproducing
kernel Hilbert space induced by the respective kernel. This corresponds to a
non-linear function in the original space.

The model learned by :class:`KernelRidge` is identical to support vector regression
(:class:`SVR`). However, different loss functions are used (ridge versus
epsilon-insensitive loss). In contrast to :class:`SVR`, fitting
:class:`KernelRidge` can be done in closed-form and is typically faster for
medium-sized datasets. On the other hand, the learned model is non-sparse and
thus slower than SVR at prediction-time.

Using :class:`KernelRidge` with a linear kernel can be advisable over 
non-kernelized Ridge in situations when the number of feature dimensions D is
considerably larger than the number of training datapoints N since the
computational cost of solving the dual (kernelized) problem is :math:`O(N^3)`
while the standard Ridge requires :math:`O(D^3)`.

The following figure compares :class:`KernelRidge` and :class:`SVR` on
an artificial dataset, which consists of a sinusoidal target function and
strong noise added to every fifth datapoint. The learned model of
:class:`KernelRidge` and :class:`SVR` is plotted, where both
complexity/regularization and bandwidth of the RBF kernel have been optimized
using grid-search. The learned functions are very similar; however, fitting
:class:`KernelRidge` is approx. seven times faster than fitting :class:`SVR`
(both with grid-search). However, prediction of 100000 target values is more
than tree times faster with SVR since it has learned a sparse model using only
approx. 1/3 of the 100 training datapoints as support vectors.

.. figure:: ../auto_examples/images/plot_kernel_ridge_regression_001.png
   :target: ../auto_examples/images/plot_kernel_ridge_regression.html
   :align: center

The next figure compares the time for fitting and prediction of
:class:`KernelRidge` and :class:`SVR` for different sizes of the training set.
Fitting :class:`KernelRidge` is faster than :class:`SVR` for medium-sized
training sets (less than 1000 samples); however, for larger training sets
:class:`SVR` scales better. With regard to prediction time, :class:`SVR` is
faster than :class:`KernelRidge` for all sizes of the training set because of
the learned sparse solution. Note that the degree of sparsity and thus the
prediction time depends on the parameters :math:`epsilon` and :math:`C` of the
:class:`SVR`.

.. figure:: ../auto_examples/images/plot_kernel_ridge_regression_002.png
   :target: ../auto_examples/images/plot_kernel_ridge_regression.html
   :align: center
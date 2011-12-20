.. _kernel_approximation:

Kernel Approximation
====================

This submodule contains functions that approximate the feature mappings that
correspond to certain kernels, as they are used for example in support vector
machines (see :ref:`svm`).
The following feature functions perform non-linear transformations of the
input, which can serve as a basis for linear classification or other
algorithms.

.. currentmodule:: sklearn.linear_model

The advantage of using approximate explicit feature maps compared to the 
`kernel trick <http://en.wikipedia.org/wiki/Kernel_trick>`_, 
which makes use of feature maps implicitly, is that explicit mappings
can be better suited for online learning and can significantly reduce the cost
of learning with very large datasets.
Standard kernelized SVMs do not scale well to large datasets, but using an
approximate kernel map it is possible to use much more efficient linear SVMs.
In particularly the combination of kernel map approximations with
:class:`SGDClassifier` can make nonlinear learning on large datasets possible.

Since there has not been much empirical work using approximate embeddings, it
is advisable to compare results against exact kernel methods when possible.


Radial Basis Function Kernel
----------------------------

.. currentmodule:: sklearn.kernel_approximation

The :class:`RBFSampler` constructs an approximate mapping for the radial basis
function kernel. 

The mapping relies on a Monte Carlo approximation to the
kernel values. The ``fit`` function performs the Monte Carlo sampling, whereas
the ``transform`` method performs the mapping of the data.  Because of the
inherent randomness of the process, results may vary between different calls to
the ``fit`` function.

The ``fit`` function takes two arguments:
`n_components`, which is the target dimensionality of the feature transform,
and `gamma`, the parameter of the RBF-kernel.  A higher `n_components` will
result in a better approximation of the kernel and will yield results more
similar to those produced by a kernel SVM. Note that "fitting" the feature
function does not actually depend on the data given to the ``fit`` function.
Only the dimensionality of the data is used.
Details on the method can be found in [RR2007]_.

.. figure:: ../auto_examples/images/plot_kernel_approximation_2.png
    :target: ../auto_examples/plot_kernel_approximation.html
    :scale: 50%
    :align: center

    Comparing an exact RBF kernel (left) with the approximation (right)

.. topic:: Examples:

    * :ref:`example_plot_kernel_approximation.py`


Additive Chi Squared Kernel
---------------------------

The chi squared kernel is a kernel on histograms, often used in computer vision.

The chi squared kernel is given by

.. math::

        k(x, y) = \sum_i \frac{2x_iy_i}{x_i+y_i}

Since the kernel is additive, it is possible to treat all components
:math:`x_i` separately for embedding. This makes it possible to sample
the Fourier transform in regular intervals, instead of approximating
using Monte Carlo sampling.

The class :class:`AdditiveChi2Kernel` implements this component wise
deterministic sampling. Each component is sampled `n` times, yielding
`2n+1` dimensions per input dimension (the multiple of two stems
from the real and complex part of the Fourier transform).
In the literature, `n` is usually choosen to be `1` or `2`, transforming
the dataset to size `n_samples x 5 * n_features` (in the case of `n=2`).

The approximate feature map provided by :class:`AdditiveChi2Sampler` can be combined
with the approximate feature map provided by :class:`RBFSampler` to yield an approximate
feature map for the exponentiated chi squared kernel.
See the [VZ2010]_ for details and [VVZ2010]_ for combination with the :class:`RBFSampler`.


Skewed Chi Squared Kernel
-------------------------

The skewed chi squared kernel is given by:

.. math::

        k(x,y) = \prod_i \frac{2\sqrt{x_i+c}\sqrt{y_i+c}}{x_i + y_i + 2c}


It has properties that are similar to the exponentiated chi squared kernel
often used in computer vision, but allows for a simple Monte Carlo
approximation of the feature map. 

The usage of the :class:`SkewedChi2Sampler` is the same as the usage described
above for the :class:`RBFSampler`. The only difference is in the free
parameter, that is called `c`.
For a motivation for this mapping and the mathematical details see [LS2010]_.


Mathematical Details
--------------------

Kernel methods like support vector machines or kernelized
PCA rely on a property of reproducing kernel Hilbert spaces.
For any positive definite kernel function `k` (a so called Mercer kernel),
it is guaranteed that there exists a mapping :math:`\phi` into a Hilber space :math:`\mathcal{H}`,
such that

.. math::

        k(x,y) = < \phi(x), \phi(y)>

Where :math:`< \cdot, \cdot >` denotes the inner product in the
Hilbert space.

If an algorithm, such as a linear support vector machine or PCA,
relies only on the scalar product of data points :math:`x_i`, one may use
the value of :math:`k(x_i, x_j)`, which corresponds to applying the algorithm
to the mapped data points :math:`\phi(x_i)`.
The advantage of using `k` is that the mapping :math:`\phi` never has
to be calculated explicitly, allowing for arbitrary large
features (even infinite).

One drawback of kernel methods is, that it might be necessary
to store many kernel values :math:`k(x_i, x_j)` during optimization.
If a kernelized classifier is applied to new data :math:`y_j`,
:math:`k(x_i, y_j)` needs to be computed to make predictions,
possibly for many different :math:`x_i` in the training set.

The classes in this submodule allow to approximate the embedding
:math:`\phi`, thereby working explicitly with the representations
:math:`\phi(x_i)`, which obviates the need to apply the kernel
or store training examples.


.. topic:: References:

    .. [RR2007] `"Random features for large-scale kernel machines"
      <http://webmail.robots.ox.ac.uk/~vgg/rg/papers/randomfeatures.pdf>`_
      Rahimi, A. and Recht, B. - Advances in neural information processing 2007,
    .. [LS2010] `"Random Fourier approximations for skewed multiplicative histogram kernels"
      <http://sminchisescu.ins.uni-bonn.de/papers/lis_dagm10.pdf>`_
      Random Fourier approximations for skewed multiplicative histogram kernels
      - Lecture Notes for Computer Sciencd (DAGM)
    .. [VZ2010] `"Efficient additive kernels via explicit feature maps"
      <http://eprints.pascal-network.org/archive/00006964/01/vedaldi10.pdf>`_
      Vedaldi, A. and Zisserman, A. - Computer Vision and Pattern Recognition 2010
    .. [VVZ2010] `"Generalized RBF feature maps for Efficient Detection"
      <http://eprints.pascal-network.org/archive/00007024/01/inproceedings.pdf.8a865c2a5421e40d.537265656b616e7468313047656e6572616c697a65642e706466.pdf>`_
      Vempati, S. and Vedaldi, A. and Zisserman, A. and Jawahar, CV - 2010

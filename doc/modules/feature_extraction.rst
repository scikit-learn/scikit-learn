
.. _feature_extraction:

==================
Feature extraction
==================

.. currentmodule:: sklearn.feature_extraction

The :mod:`sklearn.feature_extraction` module can be used to extract
features in a format supported by machine learning algorithms from datasets
consisting of formats such as text and image.

Kernel Approximation
====================

.. currentmodule:: sklearn.feature_extraction.kernel_approximation

This submodule contains functions that approximate the feature mappings that
correspond to certain kernels, as they are used for example in support vector
machines (see :ref:`svm`).

The following feature functions perform non-linear transformations of the
input, which can serve as a basis for linear classification or other
algorithms.

The advantage of using appoximate explicit feature maps compared to the kernel
trick, which makes use of feature maps implicitly, is that explicit mappings
can be better suited for online learning and can significantly reduce the cost
of learning with large datasets.

Since there has not been much empirical work using approximate embeddings, it
is advised to compare results against exact kernel methods.

Radial Basis Function Kernel
----------------------------
The :class:`RBFSampler` constructes an approximate mapping for the radial basis
function kernel. 

The mapping relies on a Monte Carlo approximation to the
kernel values. The ``fit`` function performs the Monte Carlo sampling, whereas
the ``transform`` method performs the mapping of the data.  Because of the
inherent randomness of the process, results may vary between different calls to
the func:``fit`` function.

The ``fit`` function takes two arguments:
`n_components`, which is the target dimensionality of the feature transform,
and `gamma`, the parameter of the RBF-kernel.  A higher `n_components` will
result in a better approximation of the kernel and will yield results more
similar to those produced by a kernel SVM. Note that "fitting" the feature
function does not actually depend on the data given to the ``fit`` function.
Only the dimensionality of the data is used.


Skewed Chi Squared Kernel
-------------------------
The skewed chi squared kernel is given by:

It has properties that are similar to the exponentiated chi squared kernel
often used in computer vision, but allows for a simple Monte Carlo
approximation of the feature map. 

The usage of the :class:`SkewedChi2Sampler` is the same as the usage described
above for the :class:`RBFSampler`. The only difference is in the free
parameter, that is called `c`.

Examples
--------
* :ref:`example_kernel_approximation.py`

Mathematical Details
--------------------
Kernel methods like support vector machines or kernelized
PCA rely on a property of reproducing kernel Hilbert spaces.
For any positive definite kernel function `k` (a so called Mercer kernel),
it is guaranteed that there exists a mapping :math:`\phi` into a Hilber space :math:`\mathcal{H}`,
such that
.. math::
        k(x,y) = < \phi(x), \phi(y)>
Where :math:`< \dotc, \dotc, >` denotes the inner product in the
Hilbert space.

If an algorithm, such as a linear support vector machine or PCA,
relies only on the scalar product of data points `x_i`, one may use
the value of :math:`k(x_i, x_j)`, which corresponds to applying the algorithm
to the mapped data points :math:`\phi(x_i)`.
The advantage of using `k` is that the mapping :math:`\phi` never has
to be calculated explicitly, allowing for arbitrary large
features (even infinite).

One drawback of kernel methods is, that it might be nesseccary
to store many kernel values :math:`k(x_i, x_j)` during optimization.
If a kernelized classifier is applied to new data :math:`y_j`,
:math:`k(x_i, y_j)` needs to be computed to make predictions,
possibly for many different `x_i` in the training set.

The classes in this submodule allow to approximate the embedding
:math:`\phi`, thereby working explicitly with the representations
:math:`\phi(x_i)`, which obliverates the need to apply the kernel
or store training examples.


References
----------


Text feature extraction
=======================

.. currentmodule:: sklearn.feature_extraction.text

XXX: a lot to do here


Image feature extraction
========================

.. currentmodule:: sklearn.feature_extraction.image

Patch extraction
----------------

The :func:`extract_patches_2d` function extracts patches from an image stored
as a two-dimensional array, or three-dimensional with color information along
the third axis. For rebuilding an image from all its patches, use
:func:`reconstruct_from_patches_2d`. For example let use generate a 4x4 pixel
picture with 3 color channels (e.g. in RGB format)::

    >>> import numpy as np
    >>> from sklearn.feature_extraction import image

    >>> one_image = np.arange(4 * 4 * 3).reshape((4, 4, 3))
    >>> one_image[:, :, 0]  # R channel of a fake RGB picture
    array([[ 0,  3,  6,  9],
           [12, 15, 18, 21],
           [24, 27, 30, 33],
           [36, 39, 42, 45]])

    >>> patches = image.extract_patches_2d(one_image, (2, 2), max_patches=2,
    ...     random_state=0)
    >>> patches.shape
    (2, 2, 2, 3)
    >>> patches[:, :, :, 0]
    array([[[ 0,  3],
            [12, 15]],
    <BLANKLINE>
           [[15, 18],
            [27, 30]]])
    >>> patches = image.extract_patches_2d(one_image, (2, 2))
    >>> patches.shape
    (9, 2, 2, 3)
    >>> patches[4, :, :, 0]
    array([[15, 18],
           [27, 30]])

Let us now try to reconstruct the original image from the patches by averaging
on overlapping areas::

    >>> reconstructed = image.reconstruct_from_patches_2d(patches, (4, 4, 3))
    >>> np.testing.assert_array_equal(one_image, reconstructed)

The :class:`PatchExtractor` class works in the same way as
:func:`extract_patches_2d`, only it supports multiple images as input. It is
implemented as an estimator, so it can be used in pipelines. See::

    >>> five_images = np.arange(5 * 4 * 4 * 3).reshape(5, 4, 4, 3)
    >>> patches = image.PatchExtractor((2, 2)).transform(five_images)
    >>> patches.shape
    (45, 2, 2, 3)

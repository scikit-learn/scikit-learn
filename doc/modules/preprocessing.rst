.. _preprocessing:

==================
Preprocessing data
==================

.. currentmodule:: sklearn.preprocessing

The ``sklearn.preprocessing`` package provides several common
utility functions and transformer classes to change raw feature vectors
into a representation that is more suitable for the downstream estimators.

.. _preprocessing_scaler:

Standardization or Mean Removal and Variance Scaling
====================================================

**Standardization** of datasets is a **common requirement for many
machine learning estimators** implemented in the scikit: they might behave
badly if the individual feature do not more or less look like standard
normally distributed data: Gaussian with **zero mean and unit variance**.

In practice we often ignore the shape of the distribution and just
transform the data to center it by removing the mean value of each
feature, then scale it by dividing non-constant features by their
standard deviation.

For instance, many elements used in the objective function of
a learning algorithm (such as the RBF kernel of Support Vector
Machines or the l1 and l2 regularizers of linear models) assume that
all features are centered around zero and have variance in the same
order. If a feature has a variance that is orders of magnitude larger
that others, it might dominate the objective function and make the
estimator unable to learn from other features correctly as expected.


The function :func:`scale` provides a quick and easy way to perform this
operation on a single array-like dataset::

  >>> from sklearn import preprocessing
  >>> X = [[ 1., -1.,  2.],
  ...      [ 2.,  0.,  0.],
  ...      [ 0.,  1., -1.]]
  >>> X_scaled = preprocessing.scale(X)

  >>> X_scaled                                          # doctest: +ELLIPSIS
  array([[ 0.  ..., -1.22...,  1.33...],
         [ 1.22...,  0.  ..., -0.26...],
         [-1.22...,  1.22..., -1.06...]])

..
        >>> import numpy as np
        >>> print_options = np.get_printoptions()
        >>> np.set_printoptions(suppress=True)

Scaled data has zero mean and unit variance::

  >>> X_scaled.mean(axis=0)
  array([ 0.,  0.,  0.])

  >>> X_scaled.std(axis=0)
  array([ 1.,  1.,  1.])

..    >>> print_options = np.set_printoptions(print_options)

The ``preprocessing`` module further provides a utility class
:class:`Scaler` that implements the ``Transformer`` API to compute
the mean and standard deviation on a training set so as to be
able to later reapply the same transformation on the testing set.
This class is hence suitable for use in the early steps of a
:class:`sklearn.pipeline.Pipeline`::

  >>> scaler = preprocessing.Scaler().fit(X)
  >>> scaler
  Scaler(copy=True, with_mean=True, with_std=True)

  >>> scaler.mean_                                      # doctest: +ELLIPSIS
  array([ 1. ...,  0. ...,  0.33...])

  >>> scaler.std_                                       # doctest: +ELLIPSIS
  array([ 0.81...,  0.81...,  1.24...])

  >>> scaler.transform(X)                               # doctest: +ELLIPSIS
  array([[ 0.  ..., -1.22...,  1.33...],
         [ 1.22...,  0.  ..., -0.26...],
         [-1.22...,  1.22..., -1.06...]])


The scaler instance can then be used on new data to transform it the
same way it did on the training set::

  >>> scaler.transform([[-1.,  1., 0.]])                # doctest: +ELLIPSIS
  array([[-2.44...,  1.22..., -0.26...]])

It is possible to disable either centering or scaling by either
passing ``with_mean=False`` or ``with_std=False`` to the constructor
of :class:`Scaler`.


.. topic:: References:

  Further discussion on the importance of centering and scaling data is
  available on this FAQ: `Should I normalize/standardize/rescale the data?
  <http://www.faqs.org/faqs/ai-faq/neural-nets/part2/section-16.html>`_

.. topic:: Scaling vs Whitening

  It is sometimes not enough to center and scale the features
  independently, since a downstream model can further make some assumption
  on the linear independence of the features.

  To address this issue you can use :class:`sklearn.decomposition.PCA`
  or :class:`sklearn.decomposition.RandomizedPCA` with ``whiten=True``
  to further remove the linear correlation across features.

.. topic:: Sparse input

  :func:`scale` and :class:`Scaler` accept ``scipy.sparse`` matrices
  as input **only when with_mean=False is explicitly passed to the
  constructor**. Otherwise a ``ValueError`` will be raised as
  silently centering would break the sparsity and would often crash the
  execution by allocating excessive amounts of memory unintentionally.

  If the centered data is expected to be small enough, explicitly convert
  the input to an array using the ``toarray`` method of sparse matrices
  instead.

  For sparse input the data is **converted to the Compressed Sparse Rows
  representation** (see ``scipy.sparse.csr_matrix``).
  To avoid unnecessary memory copies, it is recommended to choose the CSR
  representation upstream.


Normalization
=============

**Normalization** is the process of **scaling individual samples to have
unit norm**. This process can be useful if you plan to use a quadratic form
such as the dot-product or any other kernel to quantify the similarity
of any pair of samples.

This assumption is the base of the `Vector Space Model
<http://en.wikipedia.org/wiki/Vector_Space_Model>`_ often used in text
classification and clustering contexts.

The function :func:`normalize` provides a quick and easy way to perform this
operation on a single array-like dataset, either using the ``l1`` or ``l2``
norms::

  >>> X = [[ 1., -1.,  2.],
  ...      [ 2.,  0.,  0.],
  ...      [ 0.,  1., -1.]]
  >>> X_normalized = preprocessing.normalize(X, norm='l2')

  >>> X_normalized                                      # doctest: +ELLIPSIS
  array([[ 0.40..., -0.40...,  0.81...],
         [ 1.  ...,  0.  ...,  0.  ...],
         [ 0.  ...,  0.70..., -0.70...]])

The ``preprocessing`` module further provides a utility class
:class:`Normalizer` that implements the same operation using the
``Transformer`` API (even though the ``fit`` method is useless in this case:
the class is stateless as this operation treats samples independently).

This class is hence suitable for use in the early steps of a
:class:`sklearn.pipeline.Pipeline`::

  >>> normalizer = preprocessing.Normalizer().fit(X)  # fit does nothing
  >>> normalizer
  Normalizer(copy=True, norm='l2')


The normalizer instance can then be used on sample vectors as any transformer::

  >>> normalizer.transform(X)                            # doctest: +ELLIPSIS
  array([[ 0.40..., -0.40...,  0.81...],
         [ 1.  ...,  0.  ...,  0.  ...],
         [ 0.  ...,  0.70..., -0.70...]])

  >>> normalizer.transform([[-1.,  1., 0.]])             # doctest: +ELLIPSIS
  array([[-0.70...,  0.70...,  0.  ...]])


.. topic:: Sparse input

  :func:`normalize` and :class:`Normalizer` accept **both dense array-like
  and sparse matrices from scipy.sparse as input**.

  For sparse input the data is **converted to the Compressed Sparse Rows
  representation** (see ``scipy.sparse.csr_matrix``) before being fed to
  efficient Cython routines. To avoid unnecessary memory copies, it is
  recommended to choose the CSR representation upstream.


Binarization
============

Feature binarization
--------------------

**Feature binarization** is the process of **thresholding numerical
features to get boolean values**. This can be useful for downsteam
probabilistic estimators that make assumption that the input data
is distributed according to a multi-variate `Bernoulli distribution
<http://en.wikipedia.org/wiki/Bernoulli_distribution>`_. For instance,
this is the case for the most common class of `(Restricted) Boltzmann
Machines <http://en.wikipedia.org/wiki/Boltzmann_machine>`_
(not yet implemented in the scikit).

It is also commmon among the text processing community to use binary
feature values (probably to simplify the probabilistic reasoning) even
if normalized counts (a.k.a. term frequencies) or TF-IDF valued features
often perform slightly better in practice.

As for the :class:`Normalizer`, the utility class
:class:`Binarizer` is meant to be used in the early stages of
:class:`sklearn.pipeline.Pipeline`. The ``fit`` method does nothing
as each sample is treated independently of others::

  >>> X = [[ 1., -1.,  2.],
  ...      [ 2.,  0.,  0.],
  ...      [ 0.,  1., -1.]]

  >>> binarizer = preprocessing.Binarizer().fit(X)  # fit does nothing
  >>> binarizer
  Binarizer(copy=True, threshold=0.0)

  >>> binarizer.transform(X)
  array([[ 1.,  0.,  1.],
         [ 1.,  0.,  0.],
         [ 0.,  1.,  0.]])

It is possible to adjust the threshold of the binarizer::

  >>> binarizer = preprocessing.Binarizer(threshold=1.1)
  >>> binarizer.transform(X)
  array([[ 0.,  0.,  1.],
         [ 1.,  0.,  0.],
         [ 0.,  0.,  0.]])

As for the :class:`Scaler` and :class:`Normalizer` classes, the
preprocessing module provides a companion function :func:`binarize`
to be used when the transformer API is not necessary.

.. topic:: Sparse input

  :func:`binarize` and :class:`Binarizer` accept **both dense array-like
  and sparse matrices from scipy.sparse as input**.

  For sparse input the data is **converted to the Compressed Sparse Rows
  representation** (see ``scipy.sparse.csr_matrix``).
  To avoid unnecessary memory copies, it is recommended to choose the CSR
  representation upstream.

.. TODO

  Label binarization
  ------------------

  Please @mblondel or someone else write me!


  Kernel centering
  ================

  Please @mblondel or someone else write me!

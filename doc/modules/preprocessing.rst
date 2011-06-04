.. _preprocessing:

==================
Preprocessing data
==================

.. currentmodule:: scikits.learn.preprocessing

The ``scikits.learn.preprocessing`` package provides several common utility
functions and transformer classes to change raw feature vectors into a
representation that is more suitable for the downstream estimators.


Mean removal and variance scaling
=================================

The **Standardazition** of a dataset is a **common requirement for many
machine learning estimators** implemented in the scikit: they might behave
badly if the individual feature do not more or less look like standard
normally distributed data: Gaussian with **zero mean and unit variance**.

In practice we often ignore the shape of the distribution and just
transform the data to center it by removing the mean value of each
feature and then scale it by dividing non constant features by their
standard deviation.

For instance many elements used in the objective function of
a learning algorithm (such as the RBF kernel of Support Vector
Machines or the l1 and l2 regularizers of linear models) assume that
all features are centered around zero and have variance in the same
order. If a feature has a variance that is orders of magnitude larger
that others, it might dominate the objective function and make the
estimator unable to learn from other features correctly as expected.

The function :func:`scale` provides a quick and easy way to perform this
operation on a single array-like dataset::

  >>> from scikits.learn import preprocessing
  >>> X = [[ 1., -1.,  2.],
  ...      [ 2.,  0.,  0.],
  ...      [ 0.,  1., -1.]]
  >>> X_scaled = preprocessing.scale(X)

  >>> X_scaled                                          # doctest: +ELLIPSIS
  array([[ 0.  ..., -1.22...,  1.33...],
         [ 1.22...,  0.  ..., -0.26...],
         [-1.22...,  1.22..., -1.06...]])

Scaled data has zero mean and unit variance::

  >>> X_scaled.mean(axis=0)
  array([ 0.,  0.,  0.])

  >>> X_scaled.std(axis=0)
  array([ 1.,  1.,  1.])

The ``preprocessing`` module further provides a utility class that
implements the ``Transformer`` API to compute the mean and standard
deviation on a training set so as to be able to later reapply the same
transformation on the testing set.  This class is hence suitable for
use in the eary steps of a :class:`scikits.learn.pipeline.Pipeline`::

  >>> scaler = preprocessing.Scaler().fit(X)
  >>> scaler.mean_                                      # doctest: +ELLIPSIS
  array([ 1. ...,  0. ...,  0.33...])

  >>> scaler.std_                                       # doctest: +ELLIPSIS
  array([ 0.81...,  0.81...,  1.24...])

The scaler instance can then be used on new data to transform it the
same way it did on the training set::

  >>> scaler.transform([[-1.,  1., 0.]])                # doctest: +ELLIPSIS
  array([[-2.44...,  1.22..., -0.26...]])

.. topic:: References:

  Further discussion on the importance of centering and scaling data is
  available on this FAQ: `Should I normalize/standardize/rescale the data?
  <http://www.faqs.org/faqs/ai-faq/neural-nets/part2/section-16.html>`_

.. topic:: Notes

  It is sometimes not enough to center and scale the features independently
  since downstream model can further make assumption on the linear independence
  of the features.

  To address this issue you can use :class:`scikits.learn.decomposition.PCA`
  or :class:`scikits.learn.decomposition.RandomizedPCA` with ``whiten=True``
  to further remove the linear correlation across features.

  Also note that the current implementation of :func:`scale` and
  :class:`Scaler` **do not yet work with ``scipy.sparse`` matrices**.


Normalization
=============

TODO

.. topic:: Notes

  :func:`normalize` and :class:`Normalizer` **accept both dense array-like
  sparse matrices from scipy.sparse as input**.

  For sparse input the data is **converted to the Compressed Sparse Rows
  representation** (see ``scipy.sparse.csr_matrix``) before being fed to
  efficient cython routines. To avoid un-necessary memory copy it is
  therefore recommended to choose the CSR representation upsteam.


Binarization
============

Feature binarization
--------------------

TODO

.. topic:: Notes

  :class:`Binarizer` **accepts both dense array-like sparse matrices
  from scipy.sparse as input**.

  For sparse input the data is **converted to the Compressed Sparse Rows
  representation** (see ``scipy.sparse.csr_matrix``).
  To avoid un-necessary memory copy it is therefore recommended to choose
  the CSR representation upsteam.


.. TODO

  Label binarization
  ------------------

  write me!

..
    For doctests:

    >>> import numpy as np
    >>> import os
    >>> from sklearn import datasets
    >>> datasets.mldata.urllib2 = mock_urllib2

.. _datasets:

=========================
Dataset loading utilities
=========================

.. currentmodule:: sklearn.datasets

The ``sklearn.datasets`` package embeds some small toy datasets
as introduced in the "Getting Started" section.

To evaluate the impact of the scale of the dataset (``n_samples`` and
``n_features``) while controlling the statistical properties of the data
(typically the correlation and informativeness of the features), it is
also possible to generate synthetic data.

This package also features helpers to fetch larger datasets commonly
used by the machine learning community to benchmark algorithm on data
that comes from the 'real world'.

Toy datasets
============

scikit-learn comes with a few small standard datasets that do not
require to download any file from some external website.

.. autosummary::

   :toctree: generated/
   :template: function.rst

   load_boston
   load_iris
   load_diabetes
   load_digits
   load_linnerud

These datasets are useful to quickly illustrate the behavior of the
various algorithms implemented in the scikit. They are however often too
small to be representative of real world machine learning tasks.

Sample images
=============

The scikit also embed a couple of sample JPEG images published under Creative
Commons license by their authors. Those image can be useful to test algorithms
and pipeline on 2D data.

.. autosummary::

   load_sample_images
   load_sample_image

.. note::

  The default coding of images is based on the ``uint8`` dtype to
  spare memory.  Often machine learning algorithms work best if the
  input is converted to a floating point representation first.  Also,
  if you plan to use ``pylab.imshow`` don't forget to scale to the range
  0 - 1 as done in the following example.

.. topic:: Examples:

    * :ref:`example_cluster_plot_color_quantization.py`


.. _sample_generators:

Sample generators
=================

In addition, scikit-learn includes various random sample generators that
can be used to build artifical datasets of controled size and complexity.

.. autosummary::

   :toctree: generated/
   :template: function.rst

   make_classification
   make_regression
   make_blobs
   make_friedman1
   make_friedman2
   make_friedman3
   make_low_rank_matrix
   make_sparse_coded_signal
   make_sparse_uncorrelated
   make_spd_matrix
   make_swiss_roll
   make_s_curve


.. _libsvm_loader:

Datasets in svmlight / libsvm format
====================================

scikit-learn includes a fast utility function, ``load_svmlight_format``,  to load
datasets in the svmlight / libsvm format. In this format, each line
takes the form ``<label> <feature-id>:<feature-value>
<feature-id>:<feature-value> ...``. This format is especially suitable for sparse datasets.
Scipy sparse CSR matrices are used for ``X`` and numpy arrays are used for ``y``.

You may load a dataset like this::

  >>> from sklearn.datasets import load_svmlight_file
  >>> X_train, y_train = load_svmlight_file("/path/to/train_dataset.txt")
  ...                                                         # doctest: +SKIP

You may also load two datasets at once::

  >>> X_train, y_train, X_test, y_test = load_svmlight_file(
  ...     "/path/to/train_dataset.txt",
  ...     "/path/to/test_dataset.txt")                        # doctest: +SKIP

In this case, ``X_train`` and ``X_test`` are guaranteed to have the same number
of features. Another way to achieve the same result is to fix the number of
features::

  >>> X_test, y_test = load_svmlight_file(
  ...     "/path/to/test_dataset.txt", n_features=X_train.shape[1])
  ...                                                         # doctest: +SKIP

.. topic:: Public datasets:

 _`Public datasets in svmlight / libsvm format`: http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/


.. include:: olivetti_faces.rst 

.. include:: twenty_newsgroups.rst

.. include:: mldata.rst

.. include:: labeled_faces.rst


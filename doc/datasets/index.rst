.. _datasets:

=========================
Dataset loading utilities
=========================

.. currentmodule:: sklearn.datasets

The ``sklearn.datasets`` package embeds some small toy datasets
as introduced in the :ref:`Getting Started <loading_example_dataset>` section.

To evaluate the impact of the scale of the dataset (``n_samples`` and
``n_features``) while controlling the statistical properties of the data
(typically the correlation and informativeness of the features), it is
also possible to generate synthetic data.

This package also features helpers to fetch larger datasets commonly
used by the machine learning community to benchmark algorithm on data
that comes from the 'real world'.

General dataset API
===================

There are three distinct kinds of dataset interfaces for different types
of datasets.
The simplest one is the interface for sample images, which is described
below in the :ref:`sample_images` section.

The dataset generation functions and the svmlight loader share a simplistic
interface, returning a tuple ``(X, y)`` consisting of a n_samples x n_features
numpy array X and an array of length n_samples containing the targets y.

The toy datasets as well as the 'real world' datasets and the datasets
fetched from mldata.org have more sophisticated structure.
These functions return a dictionary-like object holding at least two items:
an array of shape ``n_samples`` * `` n_features`` with key ``data``
(except for 20newsgroups)
and a NumPy array of length ``n_features``, containing the target values,
with key ``target``.

The datasets also contain a description in ``DESCR`` and some contain
``feature_names`` and ``target_names``.
See the dataset descriptions below for details.


Toy datasets
============

scikit-learn comes with a few small standard datasets that do not
require to download any file from some external website.

.. autosummary::

   :toctree: ../modules/generated/
   :template: function.rst

   load_boston
   load_iris
   load_diabetes
   load_digits
   load_linnerud

These datasets are useful to quickly illustrate the behavior of the
various algorithms implemented in the scikit. They are however often too
small to be representative of real world machine learning tasks.

.. _sample_images:

Sample images
=============

The scikit also embed a couple of sample JPEG images published under Creative
Commons license by their authors. Those image can be useful to test algorithms
and pipeline on 2D data.

.. autosummary::

   :toctree: ../modules/generated/
   :template: function.rst

   load_sample_images
   load_sample_image

.. image:: ../auto_examples/cluster/images/plot_color_quantization_1.png
   :target: ../auto_examples/cluster/plot_color_quantization.html
   :scale: 30
   :align: right


.. warning::

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
can be used to build artificial datasets of controlled size and complexity.

.. image:: ../auto_examples/datasets/images/plot_random_dataset_1.png
   :target: ../auto_examples/datasets/plot_random_dataset.html
   :scale: 50
   :align: center

.. autosummary::

   :toctree: ../modules/generated/
   :template: function.rst

   make_classification
   make_multilabel_classification
   make_regression
   make_blobs
   make_friedman1
   make_friedman2
   make_friedman3
   make_hastie_10_2
   make_low_rank_matrix
   make_sparse_coded_signal
   make_sparse_uncorrelated
   make_spd_matrix
   make_swiss_roll
   make_s_curve
   make_sparse_spd_matrix
   make_biclusters
   make_checkerboard

.. _libsvm_loader:

Datasets in svmlight / libsvm format
====================================

scikit-learn includes utility functions for loading
datasets in the svmlight / libsvm format. In this format, each line
takes the form ``<label> <feature-id>:<feature-value>
<feature-id>:<feature-value> ...``. This format is especially suitable for sparse datasets.
In this module, scipy sparse CSR matrices are used for ``X`` and numpy arrays are used for ``y``.

You may load a dataset like as follows::

  >>> from sklearn.datasets import load_svmlight_file
  >>> X_train, y_train = load_svmlight_file("/path/to/train_dataset.txt")
  ...                                                         # doctest: +SKIP

You may also load two (or more) datasets at once::

  >>> X_train, y_train, X_test, y_test = load_svmlight_files(
  ...     ("/path/to/train_dataset.txt", "/path/to/test_dataset.txt"))
  ...                                                         # doctest: +SKIP

In this case, ``X_train`` and ``X_test`` are guaranteed to have the same number
of features. Another way to achieve the same result is to fix the number of
features::

  >>> X_test, y_test = load_svmlight_file(
  ...     "/path/to/test_dataset.txt", n_features=X_train.shape[1])
  ...                                                         # doctest: +SKIP

.. topic:: Related links:

 _`Public datasets in svmlight / libsvm format`: http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/

 _`Faster API-compatible implementation`: https://github.com/mblondel/svmlight-loader


.. include:: olivetti_faces.rst

.. include:: twenty_newsgroups.rst

.. include:: mldata.rst

.. include:: labeled_faces.rst

.. include:: covtype.rst

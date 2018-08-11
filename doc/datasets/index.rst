.. _datasets:

=========================
Dataset loading utilities
=========================

.. currentmodule:: sklearn.datasets

The ``sklearn.datasets`` package embeds some small toy datasets
as introduced in the :ref:`Getting Started <loading_example_dataset>` section.

This package also features helpers to fetch larger datasets commonly
used by the machine learning community to benchmark algorithms on data
that comes from the 'real world'.

To evaluate the impact of the scale of the dataset (``n_samples`` and
``n_features``) while controlling the statistical properties of the data
(typically the correlation and informativeness of the features), it is
also possible to generate synthetic data.

General dataset API
===================

There are three main kinds of dataset interfaces that can be used to get 
datasets depending on the desired type of dataset.
  
**The dataset loaders.** They can be used to load small standard datasets, 
described in the :ref:`toy_datasets` section.  

**The dataset fetchers.** They can be used to download and load larger datasets,
described in the :ref:`real_world_datasets` section.

Both loaders and fetchers functions return a dictionary-like object holding 
at least two items: an array of shape ``n_samples`` * ``n_features`` with 
key ``data`` (except for 20newsgroups) and a numpy array of 
length ``n_samples``, containing the target values, with key ``target``.

It's also possible for almost all of these function to constrain the output
to be a tuple containing only the data and the target, by setting the 
``return_X_y`` parameter to ``True``.

The datasets also contain a full description in their ``DESCR`` attribute and 
some contain ``feature_names`` and ``target_names``. See the dataset 
descriptions below for details.  

**The dataset generation functions.** They can be used to generate controlled 
synthetic datasets, described in the :ref:`sample_generators` section.

These functions return a tuple ``(X, y)`` consisting of a ``n_samples`` *
``n_features`` numpy array ``X`` and an array of length ``n_samples``
containing the targets ``y``.

In addition, there are also miscellanous tools to load datasets of other 
formats or from other locations, described in the :ref:`loading_other_datasets`
section. 

.. _toy_datasets:

Toy datasets
============

scikit-learn comes with a few small standard datasets that do not require to 
download any file from some external website. 

They can be loaded using the following functions:

.. autosummary::

   :toctree: ../modules/generated/
   :template: function.rst

   load_boston
   load_iris
   load_diabetes
   load_digits
   load_linnerud
   load_wine
   load_breast_cancer

These datasets are useful to quickly illustrate the behavior of the
various algorithms implemented in scikit-learn. They are however often too
small to be representative of real world machine learning tasks.

.. toctree::
    :maxdepth: 2
    :hidden:

    boston_house_prices
    iris
    diabetes
    digits
    linnerud
    wine_data
    breast_cancer

.. include:: ../../sklearn/datasets/descr/boston_house_prices.rst

.. include:: ../../sklearn/datasets/descr/iris.rst

.. include:: ../../sklearn/datasets/descr/diabetes.rst

.. include:: ../../sklearn/datasets/descr/digits.rst

.. include:: ../../sklearn/datasets/descr/linnerud.rst

.. include:: ../../sklearn/datasets/descr/wine_data.rst

.. include:: ../../sklearn/datasets/descr/breast_cancer.rst

.. _real_world_datasets:

Real world datasets
===================

scikit-learn provides tools to load larger datasets, downloading them if
necessary.

They can be loaded using the following functions:

.. autosummary::

   :toctree: ../modules/generated/
   :template: function.rst

   fetch_olivetti_faces
   fetch_20newsgroups
   fetch_20newsgroups_vectorized
   fetch_lfw_people
   fetch_lfw_pairs
   fetch_covtype
   fetch_rcv1
   fetch_kddcup99
   fetch_california_housing

.. toctree::
    :maxdepth: 2
    :hidden:

    olivetti_faces
    twenty_newsgroups
    labeled_faces
    covtype
    rcv1
    kddcup99
    california_housing

.. include:: ../../sklearn/datasets/descr/olivetti_faces.rst

.. include:: ../../sklearn/datasets/descr/twenty_newsgroups.rst

.. include:: ../../sklearn/datasets/descr/lfw.rst

.. include:: ../../sklearn/datasets/descr/covtype.rst

.. include:: ../../sklearn/datasets/descr/rcv1.rst

.. include:: ../../sklearn/datasets/descr/kddcup99.rst

.. include:: ../../sklearn/datasets/descr/california_housing.rst

.. _sample_generators:

Generated datasets
==================

In addition, scikit-learn includes various random sample generators that
can be used to build artificial datasets of controlled size and complexity.

Generators for classification and clustering
--------------------------------------------

These generators produce a matrix of features and corresponding discrete
targets.

Single label
~~~~~~~~~~~~

Both :func:`make_blobs` and :func:`make_classification` create multiclass
datasets by allocating each class one or more normally-distributed clusters of
points.  :func:`make_blobs` provides greater control regarding the centers and
standard deviations of each cluster, and is used to demonstrate clustering.
:func:`make_classification` specialises in introducing noise by way of:
correlated, redundant and uninformative features; multiple Gaussian clusters
per class; and linear transformations of the feature space.

:func:`make_gaussian_quantiles` divides a single Gaussian cluster into
near-equal-size classes separated by concentric hyperspheres.
:func:`make_hastie_10_2` generates a similar binary, 10-dimensional problem.

.. image:: ../auto_examples/datasets/images/sphx_glr_plot_random_dataset_001.png
   :target: ../auto_examples/datasets/plot_random_dataset.html
   :scale: 50
   :align: center

:func:`make_circles` and :func:`make_moons` generate 2d binary classification
datasets that are challenging to certain algorithms (e.g. centroid-based
clustering or linear classification), including optional Gaussian noise.
They are useful for visualisation. produces Gaussian
data with a spherical decision boundary for binary classification.

Multilabel
~~~~~~~~~~

:func:`make_multilabel_classification` generates random samples with multiple
labels, reflecting a bag of words drawn from a mixture of topics. The number of
topics for each document is drawn from a Poisson distribution, and the topics
themselves are drawn from a fixed random distribution. Similarly, the number of
words is drawn from Poisson, with words drawn from a multinomial, where each
topic defines a probability distribution over words. Simplifications with
respect to true bag-of-words mixtures include:

* Per-topic word distributions are independently drawn, where in reality all
  would be affected by a sparse base distribution, and would be correlated.
* For a document generated from multiple topics, all topics are weighted
  equally in generating its bag of words.
* Documents without labels words at random, rather than from a base
  distribution.

.. image:: ../auto_examples/datasets/images/sphx_glr_plot_random_multilabel_dataset_001.png
   :target: ../auto_examples/datasets/plot_random_multilabel_dataset.html
   :scale: 50
   :align: center

Biclustering
~~~~~~~~~~~~

.. autosummary::

   :toctree: ../modules/generated/
   :template: function.rst

   make_biclusters
   make_checkerboard


Generators for regression
-------------------------

:func:`make_regression` produces regression targets as an optionally-sparse
random linear combination of random features, with noise. Its informative
features may be uncorrelated, or low rank (few features account for most of the
variance).

Other regression generators generate functions deterministically from
randomized features.  :func:`make_sparse_uncorrelated` produces a target as a
linear combination of four features with fixed coefficients.
Others encode explicitly non-linear relations:
:func:`make_friedman1` is related by polynomial and sine transforms;
:func:`make_friedman2` includes feature multiplication and reciprocation; and
:func:`make_friedman3` is similar with an arctan transformation on the target.

Generators for manifold learning
--------------------------------

.. autosummary::

   :toctree: ../modules/generated/
   :template: function.rst

   make_s_curve
   make_swiss_roll

Generators for decomposition
----------------------------

.. autosummary::

   :toctree: ../modules/generated/
   :template: function.rst

   make_low_rank_matrix
   make_sparse_coded_signal
   make_spd_matrix
   make_sparse_spd_matrix


.. _loading_other_datasets:

Loading other datasets
======================

.. _sample_images:

Sample images
-------------

Scikit-learn also embed a couple of sample JPEG images published under Creative
Commons license by their authors. Those image can be useful to test algorithms
and pipeline on 2D data.

.. autosummary::

   :toctree: ../modules/generated/
   :template: function.rst

   load_sample_images
   load_sample_image

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_color_quantization_001.png
   :target: ../auto_examples/cluster/plot_color_quantization.html
   :scale: 30
   :align: right


.. warning::

  The default coding of images is based on the ``uint8`` dtype to
  spare memory.  Often machine learning algorithms work best if the
  input is converted to a floating point representation first.  Also,
  if you plan to use ``matplotlib.pyplpt.imshow`` don't forget to scale to the range
  0 - 1 as done in the following example.

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_cluster_plot_color_quantization.py`

.. _libsvm_loader:

Datasets in svmlight / libsvm format
------------------------------------

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

 _`Public datasets in svmlight / libsvm format`: https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets

 _`Faster API-compatible implementation`: https://github.com/mblondel/svmlight-loader

..
    For doctests:

    >>> import numpy as np
    >>> import os
    >>> import tempfile
    >>> # Create a temporary folder for the data fetcher
    >>> custom_data_home = tempfile.mkdtemp()
    >>> os.makedirs(os.path.join(custom_data_home, 'mldata'))


.. _mldata:

Downloading datasets from the mldata.org repository
---------------------------------------------------

`mldata.org <http://mldata.org>`_ is a public repository for machine learning
data, supported by the `PASCAL network <http://www.pascal-network.org>`_ .

The ``sklearn.datasets`` package is able to directly download data
sets from the repository using the function
:func:`sklearn.datasets.fetch_mldata`.

For example, to download the MNIST digit recognition database::

  >>> from sklearn.datasets import fetch_mldata
  >>> mnist = fetch_mldata('MNIST original', data_home=custom_data_home)

The MNIST database contains a total of 70000 examples of handwritten digits
of size 28x28 pixels, labeled from 0 to 9::

  >>> mnist.data.shape
  (70000, 784)
  >>> mnist.target.shape
  (70000,)
  >>> np.unique(mnist.target)
  array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9.])

After the first download, the dataset is cached locally in the path
specified by the ``data_home`` keyword argument, which defaults to
``~/scikit_learn_data/``::

  >>> os.listdir(os.path.join(custom_data_home, 'mldata'))
  ['mnist-original.mat']

Data sets in `mldata.org <http://mldata.org>`_ do not adhere to a strict
naming or formatting convention. :func:`sklearn.datasets.fetch_mldata` is
able to make sense of the most common cases, but allows to tailor the
defaults to individual datasets:

* The data arrays in `mldata.org <http://mldata.org>`_ are most often
  shaped as ``(n_features, n_samples)``. This is the opposite of the
  ``scikit-learn`` convention, so :func:`sklearn.datasets.fetch_mldata`
  transposes the matrix by default. The ``transpose_data`` keyword controls
  this behavior::

    >>> iris = fetch_mldata('iris', data_home=custom_data_home)
    >>> iris.data.shape
    (150, 4)
    >>> iris = fetch_mldata('iris', transpose_data=False,
    ...                     data_home=custom_data_home)
    >>> iris.data.shape
    (4, 150)

* For datasets with multiple columns, :func:`sklearn.datasets.fetch_mldata`
  tries to identify the target and data columns and rename them to ``target``
  and ``data``. This is done by looking for arrays named ``label`` and
  ``data`` in the dataset, and failing that by choosing the first array to be
  ``target`` and the second to be ``data``. This behavior can be changed with
  the ``target_name`` and ``data_name`` keywords, setting them to a specific
  name or index number (the name and order of the columns in the datasets
  can be found at its `mldata.org <http://mldata.org>`_ under the tab "Data"::

    >>> iris2 = fetch_mldata('datasets-UCI iris', target_name=1, data_name=0,
    ...                      data_home=custom_data_home)
    >>> iris3 = fetch_mldata('datasets-UCI iris', target_name='class',
    ...                      data_name='double0', data_home=custom_data_home)


..
    >>> import shutil
    >>> shutil.rmtree(custom_data_home)

.. _external_datasets:

Loading from external datasets
------------------------------

scikit-learn works on any numeric data stored as numpy arrays or scipy sparse
matrices. Other types that are convertible to numeric arrays such as pandas
DataFrame are also acceptable.
 
Here are some recommended ways to load standard columnar data into a 
format usable by scikit-learn: 

* `pandas.io <https://pandas.pydata.org/pandas-docs/stable/io.html>`_ 
  provides tools to read data from common formats including CSV, Excel, JSON
  and SQL. DataFrames may also be constructed from lists of tuples or dicts.
  Pandas handles heterogeneous data smoothly and provides tools for
  manipulation and conversion into a numeric array suitable for scikit-learn.
* `scipy.io <https://docs.scipy.org/doc/scipy/reference/io.html>`_ 
  specializes in binary formats often used in scientific computing 
  context such as .mat and .arff
* `numpy/routines.io <https://docs.scipy.org/doc/numpy/reference/routines.io.html>`_
  for standard loading of columnar data into numpy arrays
* scikit-learn's :func:`datasets.load_svmlight_file` for the svmlight or libSVM
  sparse format
* scikit-learn's :func:`datasets.load_files` for directories of text files where
  the name of each directory is the name of each category and each file inside
  of each directory corresponds to one sample from that category

For some miscellaneous data such as images, videos, and audio, you may wish to
refer to:

* `skimage.io <http://scikit-image.org/docs/dev/api/skimage.io.html>`_ or
  `Imageio <https://imageio.readthedocs.io/en/latest/userapi.html>`_ 
  for loading images and videos into numpy arrays
* `scipy.io.wavfile.read 
  <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.io.wavfile.read.html>`_ 
  for reading WAV files into a numpy array

Categorical (or nominal) features stored as strings (common in pandas DataFrames) 
will need converting to numerical features using :class:`sklearn.preprocessing.OneHotEncoder`
or :class:`sklearn.preprocessing.OrdinalEncoder` or similar.
See :ref:`preprocessing`.

Note: if you manage your own numerical data it is recommended to use an 
optimized file format such as HDF5 to reduce data load times. Various libraries
such as H5Py, PyTables and pandas provides a Python interface for reading and 
writing data in that format.

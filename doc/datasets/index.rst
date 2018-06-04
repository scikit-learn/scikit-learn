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
interface, returning a tuple ``(X, y)`` consisting of a ``n_samples`` *
``n_features`` numpy array ``X`` and an array of length ``n_samples``
containing the targets ``y``.

The toy datasets as well as the 'real world' datasets and the datasets
fetched from mldata.org have more sophisticated structure.
These functions return a dictionary-like object holding at least two items:
an array of shape ``n_samples`` * ``n_features`` with key ``data``
(except for 20newsgroups)
and a numpy array of length ``n_samples``, containing the target values,
with key ``target``.

The datasets also contain a description in ``DESCR`` and some contain
``feature_names`` and ``target_names``.
See the dataset descriptions below for details.

.. _loading_datasets:

Loading datasets
================

*quick desc*

.. _toy_datasets:

Toy datasets
------------

scikit-learn comes with a few small standard datasets that do not
require to download any file from some external website.

These datasets are useful to quickly illustrate the behavior of the
various algorithms implemented in scikit-learn. They are however often too
small to be representative of real world machine learning tasks.

.. _boston_dataset:

Boston House Prices dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Notes
+++++

Data Set Characteristics:  

    :Number of Instances: 506 

    :Number of Attributes: 13 numeric/categorical predictive
    
    :Median Value (attribute 14) is usually the target

    :Attribute Information (in order):
        - CRIM     per capita crime rate by town
        - ZN       proportion of residential land zoned for lots over 25,000 sq.ft.
        - INDUS    proportion of non-retail business acres per town
        - CHAS     Charles River dummy variable (= 1 if tract bounds river; 0 otherwise)
        - NOX      nitric oxides concentration (parts per 10 million)
        - RM       average number of rooms per dwelling
        - AGE      proportion of owner-occupied units built prior to 1940
        - DIS      weighted distances to five Boston employment centres
        - RAD      index of accessibility to radial highways
        - TAX      full-value property-tax rate per $10,000
        - PTRATIO  pupil-teacher ratio by town
        - B        1000(Bk - 0.63)^2 where Bk is the proportion of blacks by town
        - LSTAT    % lower status of the population
        - MEDV     Median value of owner-occupied homes in $1000's

    :Missing Attribute Values: None

    :Creator: Harrison, D. and Rubinfeld, D.L.

This is a copy of UCI ML housing dataset.
https://archive.ics.uci.edu/ml/machine-learning-databases/housing/


This dataset was taken from the StatLib library which is maintained at Carnegie Mellon University.

The Boston house-price data of Harrison, D. and Rubinfeld, D.L. 'Hedonic
prices and the demand for clean air', J. Environ. Economics & Management,
vol.5, 81-102, 1978.   Used in Belsley, Kuh & Welsch, 'Regression diagnostics
...', Wiley, 1980.   N.B. Various transformations are used in the table on
pages 244-261 of the latter.

The Boston house-price data has been used in many machine learning papers that address regression
problems.   
     
References
++++++++++

   - Belsley, Kuh & Welsch, 'Regression diagnostics: Identifying Influential Data and Sources of Collinearity', Wiley, 1980. 244-261.
   - Quinlan,R. (1993). Combining Instance-Based and Model-Based Learning. In Proceedings on the Tenth International Conference of Machine Learning, 236-243, University of Massachusetts, Amherst. Morgan Kaufmann.

.. _iris_dataset:

Iris Plants Database
~~~~~~~~~~~~~~~~~~~~

Notes
+++++

Data Set Characteristics:
    :Number of Instances: 150 (50 in each of three classes)
    :Number of Attributes: 4 numeric, predictive attributes and the class
    :Attribute Information:
        - sepal length in cm
        - sepal width in cm
        - petal length in cm
        - petal width in cm
        - class:
                - Iris-Setosa
                - Iris-Versicolour
                - Iris-Virginica
    :Summary Statistics:

    ============== ==== ==== ======= ===== ====================
                    Min  Max   Mean    SD   Class Correlation
    ============== ==== ==== ======= ===== ====================
    sepal length:   4.3  7.9   5.84   0.83    0.7826
    sepal width:    2.0  4.4   3.05   0.43   -0.4194
    petal length:   1.0  6.9   3.76   1.76    0.9490  (high!)
    petal width:    0.1  2.5   1.20  0.76     0.9565  (high!)
    ============== ==== ==== ======= ===== ====================

    :Missing Attribute Values: None
    :Class Distribution: 33.3% for each of 3 classes.
    :Creator: R.A. Fisher
    :Donor: Michael Marshall (MARSHALL%PLU@io.arc.nasa.gov)
    :Date: July, 1988

The famous Iris database, first used by Sir R.A. Fisher. The dataset is taken
from Fisher's paper. Note that it's the same as in R, but not as in the UCI
Machine Learning Repository, which has two wrong data points.

This is perhaps the best known database to be found in the
pattern recognition literature.  Fisher's paper is a classic in the field and
is referenced frequently to this day.  (See Duda & Hart, for example.)  The
data set contains 3 classes of 50 instances each, where each class refers to a
type of iris plant.  One class is linearly separable from the other 2; the
latter are NOT linearly separable from each other.

References
++++++++++

   - Fisher, R.A. "The use of multiple measurements in taxonomic problems"
     Annual Eugenics, 7, Part II, 179-188 (1936); also in "Contributions to
     Mathematical Statistics" (John Wiley, NY, 1950).
   - Duda, R.O., & Hart, P.E. (1973) Pattern Classification and Scene Analysis.
     (Q327.D83) John Wiley & Sons.  ISBN 0-471-22361-1.  See page 218.
   - Dasarathy, B.V. (1980) "Nosing Around the Neighborhood: A New System
     Structure and Classification Rule for Recognition in Partially Exposed
     Environments".  IEEE Transactions on Pattern Analysis and Machine
     Intelligence, Vol. PAMI-2, No. 1, 67-71.
   - Gates, G.W. (1972) "The Reduced Nearest Neighbor Rule".  IEEE Transactions
     on Information Theory, May 1972, 431-433.
   - See also: 1988 MLC Proceedings, 54-64.  Cheeseman et al"s AUTOCLASS II
     conceptual clustering system finds 3 classes in the data.
   - Many, many more ...

.. _diabetes_dataset:

Diabetes dataset
~~~~~~~~~~~~~~~~

*to fill*

.. _digits_dataset:

Optical Recognition of Handwritten Digits Data Set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*to fill*

.. _sample_images:

Sample images
-------------

.. _real_world_datasets:

Real world datasets
-------------------

.. _olivetti_faces:

The Olivetti faces dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~

`This dataset contains a set of face images`_ taken between April 1992 and April
1994 at AT&T Laboratories Cambridge. The
:func:`sklearn.datasets.fetch_olivetti_faces` function is the data
fetching / caching function that downloads the data
archive from AT&T.

.. _This dataset contains a set of face images: http://www.cl.cam.ac.uk/research/dtg/attarchive/facedatabase.html

As described on the original website:

    There are ten different images of each of 40 distinct subjects. For some
    subjects, the images were taken at different times, varying the lighting,
    facial expressions (open / closed eyes, smiling / not smiling) and facial
    details (glasses / no glasses). All the images were taken against a dark
    homogeneous background with the subjects in an upright, frontal position (with
    tolerance for some side movement).

The image is quantized to 256 grey levels and stored as unsigned 8-bit integers;
the loader will convert these to floating point values on the interval [0, 1],
which are easier to work with for many algorithms.

The "target" for this database is an integer from 0 to 39 indicating the
identity of the person pictured; however, with only 10 examples per class, this
relatively small dataset is more interesting from an unsupervised or
semi-supervised perspective.

The original dataset consisted of 92 x 112, while the version available here
consists of 64x64 images.

When using these images, please give credit to AT&T Laboratories Cambridge.

.. _20newsgroups:

The 20 newsgroups text dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 20 newsgroups dataset comprises around 18000 newsgroups posts on
20 topics split in two subsets: one for training (or development)
and the other one for testing (or for performance evaluation). The split
between the train and test set is based upon a messages posted before
and after a specific date.

This module contains two loaders. The first one,
:func:`sklearn.datasets.fetch_20newsgroups`,
returns a list of the raw texts that can be fed to text feature
extractors such as :class:`sklearn.feature_extraction.text.CountVectorizer`
with custom parameters so as to extract feature vectors.
The second one, :func:`sklearn.datasets.fetch_20newsgroups_vectorized`,
returns ready-to-use features, i.e., it is not necessary to use a feature
extractor.

Usage
+++++

The :func:`sklearn.datasets.fetch_20newsgroups` function is a data
fetching / caching functions that downloads the data archive from
the original `20 newsgroups website`_, extracts the archive contents
in the ``~/scikit_learn_data/20news_home`` folder and calls the
:func:`sklearn.datasets.load_files` on either the training or
testing set folder, or both of them::

  >>> from sklearn.datasets import fetch_20newsgroups
  >>> newsgroups_train = fetch_20newsgroups(subset='train')

  >>> from pprint import pprint
  >>> pprint(list(newsgroups_train.target_names))
  ['alt.atheism',
   'comp.graphics',
   'comp.os.ms-windows.misc',
   'comp.sys.ibm.pc.hardware',
   'comp.sys.mac.hardware',
   'comp.windows.x',
   'misc.forsale',
   'rec.autos',
   'rec.motorcycles',
   'rec.sport.baseball',
   'rec.sport.hockey',
   'sci.crypt',
   'sci.electronics',
   'sci.med',
   'sci.space',
   'soc.religion.christian',
   'talk.politics.guns',
   'talk.politics.mideast',
   'talk.politics.misc',
   'talk.religion.misc']

The real data lies in the ``filenames`` and ``target`` attributes. The target
attribute is the integer index of the category::

  >>> newsgroups_train.filenames.shape
  (11314,)
  >>> newsgroups_train.target.shape
  (11314,)
  >>> newsgroups_train.target[:10]
  array([12,  6,  9,  8,  6,  7,  9,  2, 13, 19])

It is possible to load only a sub-selection of the categories by passing the
list of the categories to load to the
:func:`sklearn.datasets.fetch_20newsgroups` function::

  >>> cats = ['alt.atheism', 'sci.space']
  >>> newsgroups_train = fetch_20newsgroups(subset='train', categories=cats)

  >>> list(newsgroups_train.target_names)
  ['alt.atheism', 'sci.space']
  >>> newsgroups_train.filenames.shape
  (1073,)
  >>> newsgroups_train.target.shape
  (1073,)
  >>> newsgroups_train.target[:10]
  array([1, 1, 1, 0, 1, 0, 0, 1, 1, 1])

Converting text to vectors
++++++++++++++++++++++++++

In order to feed predictive or clustering models with the text data,
one first need to turn the text into vectors of numerical values suitable
for statistical analysis. This can be achieved with the utilities of the
``sklearn.feature_extraction.text`` as demonstrated in the following
example that extract `TF-IDF`_ vectors of unigram tokens
from a subset of 20news::

  >>> from sklearn.feature_extraction.text import TfidfVectorizer
  >>> categories = ['alt.atheism', 'talk.religion.misc',
  ...               'comp.graphics', 'sci.space']
  >>> newsgroups_train = fetch_20newsgroups(subset='train',
  ...                                       categories=categories)
  >>> vectorizer = TfidfVectorizer()
  >>> vectors = vectorizer.fit_transform(newsgroups_train.data)
  >>> vectors.shape
  (2034, 34118)

The extracted TF-IDF vectors are very sparse, with an average of 159 non-zero
components by sample in a more than 30000-dimensional space
(less than .5% non-zero features)::

  >>> vectors.nnz / float(vectors.shape[0])
  159.01327433628319

:func:`sklearn.datasets.fetch_20newsgroups_vectorized` is a function which returns
ready-to-use tfidf features instead of file names.

.. _`20 newsgroups website`: http://people.csail.mit.edu/jrennie/20Newsgroups/
.. _`TF-IDF`: https://en.wikipedia.org/wiki/Tf-idf


Filtering text for more realistic training
++++++++++++++++++++++++++++++++++++++++++

It is easy for a classifier to overfit on particular things that appear in the
20 Newsgroups data, such as newsgroup headers. Many classifiers achieve very
high F-scores, but their results would not generalize to other documents that
aren't from this window of time.

For example, let's look at the results of a multinomial Naive Bayes classifier,
which is fast to train and achieves a decent F-score::

  >>> from sklearn.naive_bayes import MultinomialNB
  >>> from sklearn import metrics
  >>> newsgroups_test = fetch_20newsgroups(subset='test',
  ...                                      categories=categories)
  >>> vectors_test = vectorizer.transform(newsgroups_test.data)
  >>> clf = MultinomialNB(alpha=.01)
  >>> clf.fit(vectors, newsgroups_train.target)
  >>> pred = clf.predict(vectors_test)
  >>> metrics.f1_score(newsgroups_test.target, pred, average='macro')
  0.88213592402729568

(The example :ref:`sphx_glr_auto_examples_text_document_classification_20newsgroups.py` shuffles
the training and test data, instead of segmenting by time, and in that case
multinomial Naive Bayes gets a much higher F-score of 0.88. Are you suspicious
yet of what's going on inside this classifier?)

Let's take a look at what the most informative features are:

  >>> import numpy as np
  >>> def show_top10(classifier, vectorizer, categories):
  ...     feature_names = np.asarray(vectorizer.get_feature_names())
  ...     for i, category in enumerate(categories):
  ...         top10 = np.argsort(classifier.coef_[i])[-10:]
  ...         print("%s: %s" % (category, " ".join(feature_names[top10])))
  ...
  >>> show_top10(clf, vectorizer, newsgroups_train.target_names)
  alt.atheism: sgi livesey atheists writes people caltech com god keith edu
  comp.graphics: organization thanks files subject com image lines university edu graphics
  sci.space: toronto moon gov com alaska access henry nasa edu space
  talk.religion.misc: article writes kent people christian jesus sandvik edu com god

You can now see many things that these features have overfit to:

- Almost every group is distinguished by whether headers such as
  ``NNTP-Posting-Host:`` and ``Distribution:`` appear more or less often.
- Another significant feature involves whether the sender is affiliated with
  a university, as indicated either by their headers or their signature.
- The word "article" is a significant feature, based on how often people quote
  previous posts like this: "In article [article ID], [name] <[e-mail address]>
  wrote:"
- Other features match the names and e-mail addresses of particular people who
  were posting at the time.

With such an abundance of clues that distinguish newsgroups, the classifiers
barely have to identify topics from text at all, and they all perform at the
same high level.

For this reason, the functions that load 20 Newsgroups data provide a
parameter called **remove**, telling it what kinds of information to strip out
of each file. **remove** should be a tuple containing any subset of
``('headers', 'footers', 'quotes')``, telling it to remove headers, signature
blocks, and quotation blocks respectively.

  >>> newsgroups_test = fetch_20newsgroups(subset='test',
  ...                                      remove=('headers', 'footers', 'quotes'),
  ...                                      categories=categories)
  >>> vectors_test = vectorizer.transform(newsgroups_test.data)
  >>> pred = clf.predict(vectors_test)
  >>> metrics.f1_score(pred, newsgroups_test.target, average='macro')
  0.77310350681274775

This classifier lost over a lot of its F-score, just because we removed
metadata that has little to do with topic classification.
It loses even more if we also strip this metadata from the training data:

  >>> newsgroups_train = fetch_20newsgroups(subset='train',
  ...                                       remove=('headers', 'footers', 'quotes'),
  ...                                       categories=categories)
  >>> vectors = vectorizer.fit_transform(newsgroups_train.data)
  >>> clf = MultinomialNB(alpha=.01)
  >>> clf.fit(vectors, newsgroups_train.target)
  >>> vectors_test = vectorizer.transform(newsgroups_test.data)
  >>> pred = clf.predict(vectors_test)
  >>> metrics.f1_score(newsgroups_test.target, pred, average='macro')
  0.76995175184521725

Some other classifiers cope better with this harder version of the task. Try
running :ref:`sphx_glr_auto_examples_model_selection_grid_search_text_feature_extraction.py` with and without
the ``--filter`` option to compare the results.

.. topic:: Recommendation

  When evaluating text classifiers on the 20 Newsgroups data, you
  should strip newsgroup-related metadata. In scikit-learn, you can do this by
  setting ``remove=('headers', 'footers', 'quotes')``. The F-score will be
  lower because it is more realistic.

.. topic:: Examples

   * :ref:`sphx_glr_auto_examples_model_selection_grid_search_text_feature_extraction.py`

   * :ref:`sphx_glr_auto_examples_text_document_classification_20newsgroups.py`

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
will need converting to integers, and integer categorical variables may be best 
exploited when encoded as one-hot variables 
(:class:`sklearn.preprocessing.OneHotEncoder`) or similar. 
See :ref:`preprocessing`.

Note: if you manage your own numerical data it is recommended to use an 
optimized file format such as HDF5 to reduce data load times. Various libraries
such as H5Py, PyTables and pandas provides a Python interface for reading and 
writing data in that format.
.. _generating_datasets:

Generating datasets
===================

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


.. _feature_extraction:

==================
Feature extraction
==================

.. currentmodule:: sklearn.feature_extraction

The :mod:`sklearn.feature_extraction` module can be used to extract
features in a format supported by machine learning algorithms from datasets
consisting of formats such as text and image.

.. note::

   Feature extraction is very different from :ref:`feature_selection`:
   the former consists in transforming arbitrary data, such as text or
   images, into numerical features usable for machine learning. The later
   is a machine learning technique applied on these features.

.. _dict_feature_extraction:

Loading features from dicts
===========================

The class :class:`DictVectorizer` can be used to convert feature
arrays represented as lists of standard Python ``dict`` objects to the
NumPy/SciPy representation used by scikit-learn estimators.

While not particularly fast to process, Python's ``dict`` has the
advantages of being convenient to use, being sparse (absent features
need not be stored) and storing feature names in addition to values.

:class:`DictVectorizer` implements what is called one-of-K or "one-hot"
coding for categorical (aka nominal, discrete) features. Categorical
features are "attribute-value" pairs where the value is restricted
to a list of discrete of possibilities without ordering (e.g. topic
identifiers, types of objects, tags, names...).

In the following, "city" is a categorical attribute while "temperature"
is a traditional numerical feature::

  >>> measurements = [
  ...     {'city': 'Dubai', 'temperature': 33.},
  ...     {'city': 'London', 'temperature': 12.},
  ...     {'city': 'San Fransisco', 'temperature': 18.},
  ... ]

  >>> from sklearn.feature_extraction import DictVectorizer
  >>> vec = DictVectorizer()

  >>> vec.fit_transform(measurements).toarray()
  array([[  1.,   0.,   0.,  33.],
         [  0.,   1.,   0.,  12.],
         [  0.,   0.,   1.,  18.]])

  >>> vec.get_feature_names()
  ['city=Dubai', 'city=London', 'city=San Fransisco', 'temperature']

:class:`DictVectorizer` is also a useful representation transformation
for training sequence classifiers in Natural Language Processing models
that typically work by extracting feature windows around a particular
word of interest.

For example, suppose that we have a first algorithm that extracts Part of
Speech (PoS) tags that we want to use as complementary tags for training
a sequence classifier (e.g. a chunker). The following dict could be
such a window of feature extracted around the word 'sat' in the sentence
'The cat sat on the mat.'::

  >>> pos_window = [
  ...     {
  ...         'word-2': 'the',
  ...         'pos-2': 'DT',
  ...         'word-1': 'cat',
  ...         'pos-1': 'NN',
  ...         'word+1': 'on',
  ...         'pos+1': 'PP',
  ...     },
  ...     # in a real application one would extract many such dictionaries
  ... ]

This description can be vectorized into a sparse two-dimensional matrix
suitable for feeding into a classifier (maybe after being piped into a
:class:`text.TfidfTransformer` for normalization)::

  >>> vec = DictVectorizer()
  >>> pos_vectorized = vec.fit_transform(pos_window)
  >>> pos_vectorized                     # doctest: +NORMALIZE_WHITESPACE
  <1x6 sparse matrix of type '<type 'numpy.float64'>'
      with 6 stored elements in COOrdinate format>
  >>> pos_vectorized.toarray()
  array([[ 1.,  1.,  1.,  1.,  1.,  1.]])
  >>> vec.get_feature_names()
  ['pos+1=PP', 'pos-1=NN', 'pos-2=DT', 'word+1=on', 'word-1=cat', 'word-2=the']

As you can imagine, if one extracts such a context around each individual
word of a corpus of documents the resulting matrix will be very wide
(many one-hot-features) with most of them being valued to zero most
of the time. So as to make the resulting data structure able to fit in
memory the ``DictVectorizer`` class uses a ``scipy.sparse`` matrix by
default instead of a ``numpy.ndarray``.


.. _text_feature_extraction:

Text feature extraction
=======================

.. currentmodule:: sklearn.feature_extraction.text


The Bag of Words representation
-------------------------------

Text Analysis is a major application field for machine learning
algorithms. However the raw data, a sequence of symbols cannot be fed
directly to the algorithms themselves as most of them expect numerical
feature vectors with a fixed size rather than the raw text documents
with variable length.

In order to address this, scikit-learn provides utilities for the most
common ways to extract numerical features from text content, namely:

- **tokenizing** strings and giving an integer id for each possible token,
  for instance by using whitespaces and punctuation as token separators.

- **counting** the occurrences of tokens in each document.

- **normalizing** and weighting with diminishing importance tokens that
  occur in the majority of samples / documents.

In this scheme, features and samples are defined as follows:

- each **individual token occurrence frequency** (normalized or not)
  is treated as a **feature**.

- the vector of all the token frequencies for a given **document** is
  considered a multivariate **sample**.

A corpus of documents can thus be represented by a matrix with one row
per document and one column per token (e.g. word) occurring in the corpus.

We call **vectorization** the general process of turning a collection
of text documents into numerical feature vectors. This specific stragegy
(tokenization, counting and normalization) is called the **Bag of Words**
or "Bag of n-grams" representation. Documents are described by word
occurrences while completely ignoring the relative position information
of the words in the document.

When combined with :ref:`tfidf`, the bag of words encoding is also known
as the `Vector Space Model
<https://en.wikipedia.org/wiki/Vector_space_model>`_.


Sparsity
--------

As most documents will typically use a very subset of a the words used in
the corpus, the resulting matrix will have many feature values that are
zeros (typically more than 99% of them).

For instance a collection of 10,000 short text documents (such as emails)
will use a vocabulary with a size in the order of 100,000 unique words in
total while each document will use 100 to 1000 unique words individually.

In order to be able to store such a matrix in memory but also to speed
up algebraic operations matrix / vector, implementations will typically
use a sparse representation such as the implementations available in the
``scipy.sparse`` package.


Common Vectorizer usage
-----------------------

:class:`CountVectorizer` implements both tokenization and occurrence
counting in a single class::

  >>> from sklearn.feature_extraction.text import CountVectorizer

This model has many parameters, however the default values are quite
reasonable (please see  the :ref:`reference documentation
<text_feature_extraction_ref>` for the details)::

  >>> vectorizer = CountVectorizer()
  >>> vectorizer
  CountVectorizer(analyzer='word', binary=False, charset='utf-8',
          charset_error='strict', dtype=<type 'long'>, input='content',
          lowercase=True, max_df=1.0, max_features=None, max_n=1, min_n=1,
          preprocessor=None, stop_words=None, strip_accents=None,
          token_pattern=u'\\b\\w\\w+\\b', tokenizer=None, vocabulary=None)

Let's use it to tokenize and count the word occurrences of a minimalistic
corpus of text documents::

  >>> corpus = [
  ...     'This is the first document.',
  ...     'This is the second second document.',
  ...     'And the third one.',
  ...     'Is this the first document?',
  ... ]
  >>> X = vectorizer.fit_transform(corpus)
  >>> X                                       # doctest: +NORMALIZE_WHITESPACE
  <4x9 sparse matrix of type '<type 'numpy.int64'>'
      with 19 stored elements in COOrdinate format>

The default configuration tokenizes the string by extracting words of
at least 2 letters. The specific function that does this step can be
requested explicitly::

  >>> analyze = vectorizer.build_analyzer()
  >>> analyze("This is a text document to analyze.")
  [u'this', u'is', u'text', u'document', u'to', u'analyze']

Each term found by the analyzer during the fit is assigned a unique
integer index corresponding to a column in the resulting matrix. This
interpretation of the columns can be retrieved as follows::

  >>> vectorizer.get_feature_names()
  [u'and', u'document', u'first', u'is', u'one', u'second', u'the', u'third', u'this']

  >>> X.toarray()           # doctest: +ELLIPSIS
  array([[0, 1, 1, 1, 0, 0, 1, 0, 1],
         [0, 1, 0, 1, 0, 2, 1, 0, 1],
         [1, 0, 0, 0, 1, 0, 1, 1, 0],
         [0, 1, 1, 1, 0, 0, 1, 0, 1]]...)

The converse mapping from feature name to column index is stored in the
``vocabulary_`` attribute of the vectorizer::

  >>> vectorizer.vocabulary_.get('document')
  1

Hence words that were not seen in the training corpus will be completely
ignored in future calls to the transform method::

  >>> vectorizer.transform(['Something completely new.']).toarray()
  ...                           # doctest: +ELLIPSIS
  array([[0, 0, 0, 0, 0, 0, 0, 0, 0]]...)

Note that in the previous corpus, the first and the last documents have
exactly the same words hence are encoded in equal vectors. In particular
we lose the information that the last document is an interogative form. To
preserve some of the local ordering information we can extract 2-grams
of words in addition to the 1-grams (the word themselvs)::

  >>> bigram_vectorizer = CountVectorizer(min_n=1, max_n=2,
  ...                                     token_pattern=ur'\b\w+\b')
  >>> analyze = bigram_vectorizer.build_analyzer()
  >>> analyze('Bi-grams are cool!')
  [u'bi', u'grams', u'are', u'cool', u'bi grams', u'grams are', u'are cool']

The vocabulary extracted by this vectorizer is hence much bigger and
can now resolve ambiguities encoded in local positioning patterns::

  >>> X_2 = bigram_vectorizer.fit_transform(corpus).toarray()
  >>> X_2
  ...                           # doctest: +ELLIPSIS
  array([[0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0],
         [0, 0, 1, 0, 0, 1, 1, 0, 0, 2, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0],
         [1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0],
         [0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1]]...)


In particular the interogative form "Is this" is only present in the
last document::

  >>> feature_index = bigram_vectorizer.vocabulary_.get(u'is this')
  >>> X_2[:, feature_index]     # doctest: +ELLIPSIS
  array([0, 0, 0, 1]...)


.. _tfidf:

TF-IDF normalization
--------------------

In a large text corpus, some words will be very present (e.g. "the", "a",
"is" in English) hence carrying very little meaningul information about
the actual contents of the document. If we were to feed the direct count
data directly to a classifier those very frequent terms would shadow
the frequencies of rarer yet more interesting terms.

In order to re-weight the count features into floating point values
suitable for usage by a classifier it is very common to use the tf–idf
transform.

Tf means **term-frequency** while tf–idf means term-frequency times
**inverse document-frequency**. This is a orginally a term weighting
scheme developed for information retrieval (as a ranking function
for search engines results), that has also found good use in document
classification and clustering.

This normalization is implemented by the :class:`text.TfidfTransformer`
class::

  >>> from sklearn.feature_extraction.text import TfidfTransformer
  >>> transformer = TfidfTransformer()
  >>> transformer
  TfidfTransformer(norm='l2', smooth_idf=True, sublinear_tf=False, use_idf=True)

Again please see the :ref:`reference documentation
<text_feature_extraction_ref>` for the details on all the parameters.

Let's take an example with the following counts. The first term is present
100% of the time hence not very interesting. The two other features only
in less than 50% of the time hence probably more representative of the
content of the documents::

  >>> counts = [[3, 0, 1],
  ...           [2, 0, 0],
  ...           [3, 0, 0],
  ...           [4, 0, 0],
  ...           [3, 2, 0],
  ...           [3, 0, 2]]
  ...
  >>> tfidf = transformer.fit_transform(counts)
  >>> tfidf                                  # doctest: +NORMALIZE_WHITESPACE
  <6x3 sparse matrix of type '<type 'numpy.float64'>'
      with 9 stored elements in Compressed Sparse Row format>

  >>> tfidf.toarray()                        # doctest: +ELLIPSIS
  array([[ 0.85...,  0.  ...,  0.52...],
         [ 1.  ...,  0.  ...,  0.  ...],
         [ 1.  ...,  0.  ...,  0.  ...],
         [ 1.  ...,  0.  ...,  0.  ...],
         [ 0.55...,  0.83...,  0.  ...],
         [ 0.63...,  0.  ...,  0.77...]])

Each row is normalized to have unit euclidean norm. The weights of each
feature computed by the ``fit`` method call are stored in a model
attribute::

  >>> transformer.idf_                       # doctest: +ELLIPSIS
  array([ 1. ...,  2.25...,  1.84...])


As tf–idf is a very often used for text features, there is also another
class called :class:`TfidfVectorizer` that combines all the option of
:class:`CountVectorizer` and :class:`TfidfTransformer` in a single model::

  >>> from sklearn.feature_extraction.text import TfidfVectorizer
  >>> vectorizer = TfidfVectorizer()
  >>> vectorizer.fit_transform(corpus)
  ...                                       # doctest: +NORMALIZE_WHITESPACE
  <4x9 sparse matrix of type '<type 'numpy.float64'>'
      with 19 stored elements in Compressed Sparse Row format>

While the tf–idf normalization is often very useful, there might
be cases where the binary occurrence markers might offer better
features. This can be achieved by using the ``binary`` parameter
of :class:`CountVectorizer`. In particular, some estimators such as
:ref:`bernoulli_naive_bayes` explicitly model discrete boolean random
variables. Also very short text are likely to have noisy tf–idf values
while the binary occurrence info is more stable.

As usual the only way how to best adjust the feature extraction parameters
is to use a cross-validated grid search, for instance by pipelining the
feature extractor with a classifier:

 * :ref:`example_grid_search_text_feature_extraction.py`


Applications and examples
-------------------------

The bag of words representation is quite simplistic but surprisingly
useful in practice.

In particular in a **supervised setting** it can be successfully combined
with fast and scalable linear models to train **document classificers**,
for instance:

 * :ref:`example_document_classification_20newsgroups.py`

In an **unsupervised setting** it can be used to group similar documents
together by applying clustering algorithms such as :ref:`k_means`:

  * :ref:`example_document_clustering.py`

Finally it is possible to discover the main topics of a corpus by
relaxing the hard assignement constraint of clustering, for instance by
using :ref:`NMF`:

  * :ref:`example_applications_topics_extraction_with_nmf.py`


Limitations of the Bag of Words representation
----------------------------------------------

While some local positioning information can be preserved by extracting
n-grams instead of individual words, Bag of Words and Bag of n-grams
destroy most of the inner structure of the document and hence most of
the meaning carried by that internal structure.

In order to address the wider task of Natural Language Understanding,
the local structure of sentences and paragraphs should thus be taken
into account. Many such models will thus be casted as "Structured output"
problems which are currently outside of the scope of scikit-learn.


Customizing the vectorizer classes
-----------------------------------

It is possible to customize the behavior by passing some callable as
parameters of the vectorizer::

  >>> def my_tokenizer(s):
  ...     return s.split()
  ...
  >>> vectorizer = CountVectorizer(tokenizer=my_tokenizer)
  >>> vectorizer.build_analyzer()(u"Some... punctuation!")
  [u'some...', u'punctuation!']

In particular we name:

  * ``preprocessor`` a callable that takes a string as input and return
    another string (removing HTML tags or converting to lower case for
    instance)

  * ``tokenizer`` a callable that takes a string as input and output a
    sequence of feature occurrences (a.k.a. the tokens).

  * ``analyzer`` a callable that wraps calls to the preprocessor and
    tokenizer and further perform some filtering or n-grams extractions
    on the tokens.

To make the preprocessor, tokenizer and analyzers aware of the model
parameters it is possible to derive from the class and override the
``build_preprocessor``, ``build_tokenizer``` and ``build_analyzer``
factory method instead.

Customizing the vectorizer can be very useful to handle Asian languages
that do not use an explicit word separator such as the whitespace for
instance.


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

Connectivity graph of an image
-------------------------------

Several estimators in the scikit-learn can use connectivity information between
features or samples. For instance Ward clustering
(:ref:`hierarchical_clustering`) can cluster together only neighboring pixels
of an image, thus forming contiguous patches:

.. figure:: ../auto_examples/cluster/images/plot_lena_ward_segmentation_1.png
   :target: ../auto_examples/cluster/plot_lena_ward_segmentation.html
   :align: center
   :scale: 40

For this purpose, the estimators use a 'connectivity' matrix, giving
which samples are connected.

The function :func:`img_to_graph` returns such a matrix from a 2D or 3D
image. Similarly, :func:`grid_to_graph` build a connectivity matrix for
images given the shape of these image.

These matrices can be used to impose connectivity in estimators that use
connectivity information, such as Ward clustering
(:ref:`hierarchical_clustering`), but also to build precomputed kernels,
or similarity matrices.

.. note:: **Examples**

   * :ref:`example_cluster_plot_lena_ward_segmentation.py`

   * :ref:`example_cluster_plot_segmentation_toy.py`

   * :ref:`example_cluster_plot_feature_agglomeration_vs_univariate_selection.py`


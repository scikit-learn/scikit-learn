
.. _feature_extraction:

==================
Feature extraction
==================

.. currentmodule:: sklearn.feature_extraction

The :mod:`sklearn.feature_extraction` module can be used to extract
features in a format supported by machine learning algorithms from datasets
consisting of formats such as text and image.


Text feature extraction
=======================

.. currentmodule:: sklearn.feature_extraction.text


The Bags of Words representation
--------------------------------

Natural Language Processing is a major application field for machine
learning algorithms. However the raw data, a sequence of symbols cannot
be fed directly to the algorithms themselves as most of them expect
numerical feature vectors with a fixed size rather than the raw text
documents with variable length.

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

We call **vectorization** the general process of turning a collection of
text documents into numerical feature vectors. This specific stragegy
(tokenization, counting and normalization) is called the **Bags of
Words** representation as documents are descriped by word occurrences
while completely ignoring the relative position information of the words
in the document.


Sparsity
--------

As most documents will typically use a very subset of a the words used in
the corpus, the resulting matrix will have many feature values that are
zeros (typically more than 99% of them).

For instance a collection of 10,000 short text documents (such as emails)
will use a vocabulary with a size in the order of 100,000 unique words in
total while each document will use 100 to 1000 unique words individually.

In order to be able to store this such a matrix in memory but also to
speed up algebraic operations matrix / vector, implementations will
typically use a sparse representation such as the implementations
available in the ``scipy.sparse`` package.


Common Vectorizer usage
-----------------------

:class:`CountVectorizer` implements both tokenization and occurrence
counting in a single class::

  >>> from sklearn.feature_extraction.text import CountVectorizer

This model has many parameters, however the default values are quite
reasonable (please refer to the :ref:`reference documentation
<text_feature_extraction_ref>` for the details)::

  >>> vectorizer = CountVectorizer()
  >>> vectorizer
  CountVectorizer(analyzer='word', binary=False, charset='utf-8',
          charset_error='strict', dtype=<type 'long'>, input='content',
          lowercase=True, max_df=1.0, max_features=None, max_n=1, min_n=1,
          preprocessor=None, stop_words=None, strip_accents='ascii',
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

The default configuration tokenize the string by extracting words of
at least 2 letters. The specific function that does this step can be
requested explicitly::

  >>> analyze = vectorizer.build_analyzer()
  >>> analyze("This is a text document to analyze.")
  ['this', 'is', 'text', 'document', 'to', 'analyze']

Each term found by the analyzer during the fit is assigned a unique
integer index to assign it a column in the resulting matrix.  This
interpretation of the columns can be retrieved as follows::

  >>> list(vectorizer.get_feature_names())
  ['and', 'document', 'first', 'is', 'one', 'second', 'the', 'third', 'this']

  >>> X.toarray()
  array([[0, 1, 1, 1, 0, 0, 1, 0, 1],
         [0, 1, 0, 1, 0, 2, 1, 0, 1],
         [1, 0, 0, 0, 1, 0, 1, 1, 0],
         [0, 1, 1, 1, 0, 0, 1, 0, 1]])

The converse mapping from feature name to column index is stored in the
``vocabulary_`` attribute of the vectorizer::

  >>> vectorizer.vocabulary_.get('document')
  1

Hence words that were not seen in the training corpus will be completely
ignored in future calls to the transform method::

  >>> vectorizer.transform(['Something completely new.']).toarray()
  array([[0, 0, 0, 0, 0, 0, 0, 0, 0]])

Note that in the previous corpus, the first and the last documents have
exaclty the same words hence are encoded in equal vectors. In particular
we lose the info that the last document is an interogative form. To
preserve some of the local ordering information we can extract 2-grams
of words in addition to the 1-grams (the word themselvs)::

  >>> bigram_vectorizer = CountVectorizer(min_n=1, max_n=2,
  ...                                     token_pattern=ur'\b\w+\b')
  >>> analyze = bigram_vectorizer.build_analyzer()
  >>> analyze('Bi-grams are cool!')
  [u'bi', u'grams', u'are', u'cool', u'bi grams', u'grams are', u'are cool']

The vocabulary extracted by this vectorizer is hence much bigger and
can now resolve ambiguities encoded in local positioning patterns::

  >>> X_2 = bigram_vectorizer.fit_transform(corpus)
  >>> X_2.toarray()
  array([[0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0],
         [0, 0, 1, 0, 0, 1, 1, 0, 0, 2, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0],
         [1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0],
         [0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1]])


In particular the interogative form "Is this" is only present in the
last document::

  >>> bigram_vectorizer.vocabulary_.get(u'is this')
  7


TF-IDF normalization
--------------------

TODO: Explain

TODO: Limitations of TF-IDF: short text, binary occurrences better

TODO: Point to grid search example.


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

Finally is possible to discover the main topics of a corpus by relaxing
the hard assignement constraint of clustering, for instance by using
:ref:`NMF`:

  * :ref:`example_applications_topics_extraction_with_nmf.py`


Limitations of the Bag of Words representation
----------------------------------------------

While some local positioning information can be preserved by extracting
n-grams instead of individual words, Bag of Words and Bag of n-grams
destroy most of the inner structure of the document and hence most of
the meaning carried by that internal structure.

In order to address the wider task of Natural Language Understanding,
the local structure of sentences and paragraphs should thus be taken
into account. Many such models will thus be casted as "Structured Ouput"
problems which are currently outside of the scope of scikit-learn.


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

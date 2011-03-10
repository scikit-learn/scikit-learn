Working with text data
======================

The goal of this section is to explore some of the main ``scikit-learn``
tools on a single practical task: analysing a collection of text
documents (newsgroups posts) on twenty different topics.

In this section we will see how to:

  - load the file contents and the categories

  - extract feature vectors suitable for machine learning

  - train a linear model to perform categorization

  - use a grid search strategy to find a good configuration of both
    the feature extraction components and the classifier


Downloading the data and loading it from Python
-----------------------------------------------

The dataset is called "Twenty Newsgroups". Here is the official
description, quoted from the `website
<http://people.csail.mit.edu/jrennie/20Newsgroups/>`_:

  The 20 Newsgroups data set is a collection of approximately 20,000
  newsgroup documents, partitioned (nearly) evenly across 20 different
  newsgroups. To the best of my knowledge, it was originally collected
  by Ken Lang, probably for his Newsweeder: Learning to filter
  netnews paper, though he does not explicitly mention this collection.
  The 20 newsgroups collection has become a popular data set for
  experiments in text applications of machine learning techniques,
  such as text classification and text clustering.

To download the dataset, go to ``$TUTORIAL_HOME/twenty_newsgroups``
run the ``fetch_data.py`` script.

Once the data is downloaded, fire an ipython shell from and define
a variable to hold the list of categories to load. In order to get
fast execution times for this first example we will work on a partial
dataset with only 4 categories out of the 20 available in the
dataset::

  >>> categories = ['alt.atheism', 'soc.religion.christian',
  ...               'comp.graphics', 'sci.med']

We can now load the list of files matching those categories as follows::

  >>> from scikits.learn.datasets import load_files
  >>> twenty_train = load_files('data/twenty_newsgroups/20news-bydate-train',
  ...                           categories=categories)


The returned dataset is a ``scikit-learn`` "bunch": a simple holder
object with fields that can be both accessed as python ``dict``
keys or ``object`` attributes for convenience, for instance the
``target_names`` holds the list of the requested category names::

  >>> twenty_train.target_names
  ['alt.atheism', 'comp.graphics', 'sci.med', 'soc.religion.christian']

The files them-selves are not loaded in memory yet::

  >>> twenty_train.filenames.shape
  (2257,)
  >>> twenty_train.filenames[0]
  'data/twenty_newsgroups/20news-bydate-train/comp.graphics/38244'

Let us print the first 2 lines of the first file::

  >>> print "".join(open(twenty_train.filenames[0]).readlines()[:2]).strip()
  From: clipper@mccarthy.csd.uwo.ca (Khun Yee Fung)
  Subject: Re: looking for circle algorithm faster than Bresenhams

Supervised learning algorithms will require the category to predict
for each document. In this case the category is the name of the
newsgroup which also happens to be the name of folder holding the
individual documents.

For speed and space efficiency reasons ``scikit-learn`` loads the
target attribute as an array of integers that corresponds to the
index of the category name in the ``target_names`` list. The category
integer id of each sample is stored in the ``target`` attribute::

  >>> twenty_train.target[:10]
  array([1, 0, 2, 2, 0, 1, 1, 3, 3, 2])

It is possible to get back the category names as follows::

  >>> for t in twenty_train.target[:10]:
  ...     print twenty_train.target_names[t]
  ...
  comp.graphics
  alt.atheism
  sci.med
  sci.med
  alt.atheism
  comp.graphics
  comp.graphics
  soc.religion.christian
  soc.religion.christian
  sci.med

You can notices that the samples have been shuffled randomly (with
a fixed RNG seed): this is useful if you select only the first
samples to quickly train a model and get a first idea of the results
before re-training on the complete dataset later.


Extracting features from text files
-----------------------------------

In order to perform machine learning on text documents, one first
need to turn the text content into numerical feature vectors.


Bags of words
~~~~~~~~~~~~~

The most intuitive way to do so is the bags of words representation:

  1. assign a fixed integer id to each word occurring in any document
     of the training set (for instance by building a dictionary
     from words to integer indices).

  2. for each document #i, count the number of occurrences of each
     word w and store it in ``X[i, j]`` as the the value of feature
     #j where j is the index of word w in the dictionary

The bags of words representation implies that ``n_features`` is
the number of distinct words in the corpus: this number is typically
larger that 100,000.

If ``n_samples == 10000``, storing ``X`` as a numpy array of type
float32 would require 10000 x 100000 x 4 bytes = **4GB in RAM** which
is barely manageable on today's computers.

Furtunately, **most values in X will be zeros** since for a given
document less than a couple thousands of distinct words will be
used. For this reason we say that bags of words are typically
**high-dimensional sparse datasets**.

Hence it is highly recommended to use ``scipy.sparse`` matrices
instead of numpy arrays to store the extracted features of a text
corpus.


Tokenizing text with ``scikit-learn``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``scikit-learn`` offers a couple of basic yet useful utilities to
work with text data. The first one is a preprocessor that removes
accents and convert to lowercase on roman languages::

  >>> from scikits.learn.feature_extraction.text import RomanPreprocessor
  >>> text = u"J'ai bien mang\xe9."
  >>> print RomanPreprocessor().preprocess(text)
  j'ai bien mange.

The second one is a utility that extract that splits the text into words after
having applied the preprocessor::

  >>> from scikits.learn.feature_extraction.text import WordNGramAnalyzer
  >>> WordNGramAnalyzer().analyze(text)
  ['ai', 'bien', 'mange']

Note that punctuation and single letter words have automatically
been removed.

It is further possible to configure ``WordNGramAnalyzer`` to extract n-grams
instead of single words::

  >>> WordNGramAnalyzer(min_n=1, max_n=2).analyze(text)
  [u'ai', u'bien', u'mange', u'ai bien', u'bien mange']

These tools are wrapped into a higher level component that is able to build a
dictionary of features::

  >>> from scikits.learn.feature_extraction.text.sparse import CountVectorizer
  >>> count_vect = CountVectorizer()
  >>> _ = count_vect.fit(open(f).read() for f in twenty_train.filenames)

Once fitted, the vectorizer has build a dictionary of feature indices::

  >>> count_vect.vocabulary.get(u'algorithm')
  1513

The index value of a word in the vocabulary is linked to its frequency
in the whole training corpus.

Once the vocabulary is built, it is possible to rescan the training
set so as to perform the actual feature extraction::

  >>> X_train_counts = count_vect.transform(
  ...     open(f).read() for f in twenty_train.filenames)
  >>> X_train_counts.shape
  (2257, 33881)


From occurrences to frequencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Occurrence count is a good start but there is an issue: longer
documents will have higher average count values than shorter document,
even though they might talk about the same topics.

To avoid these potential discrepancies is suffice to divide the
number of occurrences of each words in document by the total number
of words in the document: this new features are called TF for Term
Frequency.

Another refinement on top of TF is to downscale term weights for
words that occur in many documents in the corpus and are therefore
less informative than those that occur only in a smaller portion of the
corpus.

This downscaling is called `TF-IDF`_ for "Term Frequency times
Inverse Document Frequency".

.. _`TF-IDF`: http://en.wikipedia.org/wiki/Tf-idf


Both TF and TF-IDF can be computed as follows::

  >>> from scikits.learn.feature_extraction.text.sparse import TfidfTransformer
  >>> tf_transformer = TfidfTransformer(use_idf=False).fit(X_train_counts)
  >>> X_train_tf = tf_transformer.transform(X_train_counts)
  >>> X_train_tf.shape
  (2257, 33881)

  >>> tfidf_transformer = TfidfTransformer().fit(X_train_counts)
  >>> X_train_tfidf = tfidf_transformer.transform(X_train_counts)
  >>> X_train_tfidf.shape
  (2257, 33881)






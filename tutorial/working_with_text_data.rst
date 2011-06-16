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
  newsgroups. To the best of our knowledge, it was originally collected
  by Ken Lang, probably for his paper "Newsweeder: Learning to filter
  netnews," though he does not explicitly mention this collection.
  The 20 newsgroups collection has become a popular data set for
  experiments in text applications of machine learning techniques,
  such as text classification and text clustering.

To download the dataset, go to ``$TUTORIAL_HOME/data/twenty_newsgroups``
and run the ``fetch_data.py`` script.

Once the data is downloaded, start a Python interpreter (or IPython shell)
in the ``$TUTORIAL_HOME`` folder and define a variable to hold the list
of categories to load. In order to get fast execution times for
this first example we will work on a partial dataset with only 4
categories out of the 20 available in the dataset::

  >>> categories = ['alt.atheism', 'soc.religion.christian',
  ...               'comp.graphics', 'sci.med']

We can now load the list of files matching those categories as follows::

  >>> from scikits.learn.datasets import load_filenames
  >>> twenty_train = load_filenames('data/twenty_newsgroups/20news-bydate-train',
  ...                           categories=categories)


The returned dataset is a ``scikit-learn`` "bunch": a simple holder
object with fields that can be both accessed as python ``dict``
keys or ``object`` attributes for convenience, for instance the
``target_names`` holds the list of the requested category names::

  >>> twenty_train.target_names
  ['alt.atheism', 'comp.graphics', 'sci.med', 'soc.religion.christian']

The files themselves are not loaded in memory yet::

  >>> twenty_train.filenames.shape
  (2257,)
  >>> twenty_train.filenames[0]
  'data/twenty_newsgroups/20news-bydate-train/comp.graphics/38244'

Let's print the first 2 lines of the first file::

  >>> print "".join(open(twenty_train.filenames[0]).readlines()[:2]).strip()
  From: clipper@mccarthy.csd.uwo.ca (Khun Yee Fung)
  Subject: Re: looking for circle algorithm faster than Bresenhams

Supervised learning algorithms will require a category label for each
document in the training set. In this case the category is the name of the
newsgroup which also happens to be the name of the folder holding the
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

You can notice that the samples have been shuffled randomly (with
a fixed RNG seed): this is useful if you select only the first
samples to quickly train a model and get a first idea of the results
before re-training on the complete dataset later.


Extracting features from text files
-----------------------------------

In order to perform machine learning on text documents, we first need to
turn the text content into numerical feature vectors.


Bags of words
~~~~~~~~~~~~~

The most intuitive way to do so is the bags of words representation:

  1. assign a fixed integer id to each word occurring in any document
     of the training set (for instance by building a dictionary
     from words to integer indices).

  2. for each document #i, count the number of occurrences of each
     word w and store it in ``X[i, j]`` as the value of feature
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
**high-dimensional sparse datasets**. We can save a lot of memory by
only storing the non-zero parts of the feature vectors in memory.

``scipy.sparse`` matrices are data structures that do exactly this,
and ``scikit-learn`` has built-in support for these structures.


Tokenizing text with ``scikit-learn``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``scikit-learn`` offers a couple of basic yet useful utilities to
work with text data. The first one is a preprocessor that removes
accents and converts to lowercase on roman languages::

  >>> from scikits.learn.feature_extraction.text import RomanPreprocessor
  >>> text = u"J'ai bien mang\xe9."
  >>> print RomanPreprocessor().preprocess(text)
  j'ai bien mange.

The second one is a utility that splits the text into words after
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
dictionary of features and transform documents to feature vectors::

  >>> from scikits.learn.feature_extraction.text import CountVectorizer
  >>> count_vect = CountVectorizer()
  >>> docs_train = [open(f).read() for f in twenty_train.filenames]
  >>> X_train_counts = count_vect.fit_transform(docs_train)
  >>> X_train_counts.shape
  (2257, 33881)

Once fitted, the vectorizer has built a dictionary of feature indices::

  >>> count_vect.vocabulary.get(u'algorithm')
  1513

The index value of a word in the vocabulary is linked to its frequency
in the whole training corpus.

.. note:

  The method ``count_vect.fit_transform`` performs two actions:
  it learns the vocabulary and transforms the documents into count vectors.
  It's possible to separate these steps by calling
  ``count_vect.fit(docs_train)`` followed by
  ``X_train_counts = count_vect.transform(docs_train)``,
  but doing so would read and tokenize each text file twice.


From occurrences to frequencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Occurrence count is a good start but there is an issue: longer
documents will have higher average count values than shorter documents,
even though they might talk about the same topics.

To avoid these potential discrepancies it suffices to divide the
number of occurrences of each word in a document by the total number
of words in the document: these new features are called "tf" for Term
Frequencies.

Another refinement on top of tf is to downscale weights for words
that occur in many documents in the corpus and are therefore less
informative than those that occur only in a smaller portion of the
corpus.

This downscaling is called `tf–idf`_ for "Term Frequency times
Inverse Document Frequency".

.. _`tf–idf`: http://en.wikipedia.org/wiki/Tf–idf


Both tf and tf–idf can be computed as follows::

  >>> from scikits.learn.feature_extraction.text import TfidfTransformer
  >>> tf_transformer = TfidfTransformer(use_idf=False).fit(X_train_counts)
  >>> X_train_tf = tf_transformer.transform(X_train_counts)
  >>> X_train_tf.shape
  (2257, 33881)

  >>> tfidf_transformer = TfidfTransformer()
  >>> X_train_tfidf = tfidf_transformer.fit_transform(X_train_counts)
  >>> X_train_tfidf.shape
  (2257, 33881)


Training a linear classifier
----------------------------

Now that we have our feature, we can train a classifier to try to predict
the category of a post. Let's start with a naïve Bayes classifier, which
provides a nice baseline for this task. ``scikit-learn`` includes several
variants of this classifier; the one most suitable for word counts is the
multinomial variant::

  >>> from scikits.learn.naive_bayes import MultinomialNB
  >>> clf = MultinomialNB().fit(X_train_tfidf, twenty_train.target)

To try to predict the outcome on a new document we need to extract
the features using almost the same feature extracting chain as before.
The difference is that we call ``transform`` instead of ``fit_transform``
on the transformers, since they have already been fit to the training set::


  >>> docs_new = ['God is love', 'OpenGL on the GPU is fast']
  >>> X_new_counts = count_vect.transform(docs_new)
  >>> X_new_tfidf = tfidf_transformer.transform(X_new_counts)

  >>> predicted = clf.predict(X_new_tfidf)

  >>> for doc, category in zip(docs_new, predicted):
  ...     print '%r => %s' % (doc, twenty_train.target_names[category])
  ...
  'God is love' => soc.religion.christian
  'OpenGL on the GPU is fast' => comp.graphics


Building a pipeline
-------------------

In order to make the vectorizer => transformer => classifier easier
to work with, ``scikit-learn`` provides a ``Pipeline`` class that behaves
like a compound classifier::

  >>> from scikits.learn.pipeline import Pipeline
  >>> text_clf = Pipeline([
  ...     ('vect', CountVectorizer()),
  ...     ('tfidf', TfidfTransformer()),
  ...     ('clf', MultinomialNB()),
  ... ])

The names ``vect``, ``tfidf`` and ``clf`` (classifier) are arbitrary.
We shall see their use in the section on grid search, below.
We can now train the model with a single command::

  >>> _ = text_clf.fit(docs_train, twenty_train.target)


Evaluation of the performance on the test set
---------------------------------------------

Evaluating the predictive accuracy of the model is equally easy::

  >>> import numpy as np
  >>> twenty_test = load_filenames('data/twenty_newsgroups/20news-bydate-test',
  ...                              categories=categories)
  ...
  >>> docs_test = [open(f).read() for f in twenty_test.filenames]
  >>> predicted = text_clf.predict(docs_test)
  >>> np.mean(predicted == twenty_test.target)
  0.86884154460719043

I.e., we achieved 86.9% accuracy. Let's see if we can do better with a
linear support vector machine (SVM), which is widely regarded as one of
the best text classification algorithms (although it's also a bit slower
than naïve Bayes). We can change the learner by just plugging a different
classifier object into our pipeline::

  >>> from scikits.learn.svm.sparse import LinearSVC
  >>> text_clf = Pipeline([
  ...     ('vect', CountVectorizer()),
  ...     ('tfidf', TfidfTransformer()),
  ...     ('clf', LinearSVC()),
  ... ])
  0.92410119840213045

``scikit-learn`` further provides utilities for more detailed performance
analysis of the results::

  >>> from scikits.learn import metrics
  >>> print metrics.classification_report(
  ...     twenty_test.target, predicted,
  ...     target_names=twenty_test.target_names)
  ...

                          precision    recall  f1-score   support
  
             alt.atheism       0.95      0.80      0.87       319
           comp.graphics       0.96      0.97      0.96       389
                 sci.med       0.95      0.95      0.95       396
  soc.religion.christian       0.86      0.96      0.90       398
  
             avg / total       0.93      0.92      0.92      1502
  
  >>> metrics.confusion_matrix(twenty_test.target, predicted)
  array([[254,   4,  11,  50],
         [  3, 376,   6,   4],
         [  1,   9, 377,   9],
         [  9,   4,   4, 381]])


.. note:

  SVC stands for support vector classifier. ``scikit-learn`` also
  includes support vector machine for regression tasks, which are
  called SVR.

Parameter tuning using grid search
----------------------------------

We've already encountered some parameters such as ``use_idf`` in the
``TfidfTransformer``. Classifiers tend to have many parameters as well;
e.g., ``MultinomialNB`` includes a smoothing parameter ``alpha``
and ``LinearSVC`` has a penalty parameter ``C``
(see the module documentation, or use the Python ``help`` function,
to get a description of these).

Instead of tweaking the parameters of the various components of the
chain, it is possible to run an exhaustive search of the best
parameters on a grid of possible values. We try out all classifiers
on either words or bigrams, with or without idf, and with a penalty
parameter of either 100 or 1000 for the linear SVM::

  >>> from scikits.learn.grid_search import GridSearchCV
  >>> parameters = {
  ...     'vect__analyzer__max_n': (1, 2),
  ...     'tfidf__use_idf': (True, False),
  ...     'clf__C': (100, 1000),
  ... }

Obviously, such an exhaustive search can be expensive. If we have multiple
CPU cores at our disposal, we can tell the grid searcher to try these eight
parameter combinations in parallel with the ``n_jobs`` parameter. If we give
this parameter a value of ``-1``, grid search will detect how many cores are installed and uses them all::

  >>> gs_clf = GridSearchCV(text_clf, parameters, n_jobs=-1)

The grid search instance behaves like a normal ``scikit-learn``
model. Let's perform the search on a smaller subset of the training data
to speed up the computation::

  >>> gs_clf = gs_clf.fit(docs_train[:400], twenty_train.target[:400])

The result of calling ``fit`` on a ``GridSearchCV`` object is a classifier
that we can use to ``predict``::

  >>> twenty_train.target_names[gs_clf.predict(['God is love'])]
  'alt.atheism'

but otherwise, it's a pretty large and clumsy object. We can, however, get the
optimal parameters out by inspecting the object's ``grid_scores_`` attribute,
which is a list of parameters/score pairs. To get the best scoring attributes,
we can do::

  >>> best_parameters, score = max(gs_clf.grid_scores_, key=lambda x: x[1])
  >>> for param_name in sorted(parameters.keys()):
  ...     print "%s: %r" % (param_name, best_parameters[param_name])
  ...
  clf__C: 100
  tfidf__use_idf: True
  vect__analyzer__max_n: 2
  >>> score
  0.92507387872666735

.. note:

  A ``GridSearchCV`` object also stores the best classifier that it trained
  as its ``best_estimator`` attribute. In this case, that isn't much use as
  we trained on a small, 400-document subset of our full training set.



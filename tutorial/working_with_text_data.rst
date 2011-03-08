Working with text data
======================

The goal of this section is to explore some of the main ``scikit-learn``
tools on a single practical task: analysing a collection of  text
documents (newsgroups posts) on twenty different topics.


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
a variable to hold the list of categories to load (we will work on
a partial dataset to work faster)::

  >>> categories = ['alt.atheism', 'soc.religion.christian',
  ...               'comp.graphics', 'sci.med']

.. note::

  To be able to copy and paste examples without taking care of the leading
  ``>>>`` and ``...`` prompt signs, enable the ipython doctest mode with:
  ``%doctest_mode``

We can now load the list of files matching those categories as follows::

  >>> from scikits.learn.datasets import load_files
  >>> training_set = load_files('data/twenty_newsgroups/20news-bydate-train',
  ...                           categories=categories)

.. note::

    Whenever you import a scikit-learn class or function of the first time,
    you are advised to read the docstring by using the ``?`` magic suffix
    of ipython, for instance type: ``load_files?``.


The returned dataset is a ``scikit-learn`` "bunch": a simple class
holder with fields that can be both accessed as python ``dict``
keys or ``object`` attributes for convenience, for instance the
``target_names`` holds the list of the requested category names::

  >>> training_set.target_names
  ['alt.atheism', 'comp.graphics', 'sci.med', 'soc.religion.christian']

The files them-selves are not loaded in memory yet::

  >>> training_set.filenames.shape
  (2257,)
  >>> training_set.filenames[0]
  'data/twenty_newsgroups/20news-bydate-train/comp.graphics/38244'

  >>> print open(training_set.filenames[0]).read()[:114]
  From: clipper@mccarthy.csd.uwo.ca (Khun Yee Fung)
  Subject: Re: looking for circle algorithm faster than Bresenhams

Supervised learning algorithms will require information on the right
category to predict for document. In this case the category is the
name of the newsgroups, which also happen to be the name of folder
holding the individual documents.

For speed and space efficiency ``scikit-learn`` loads the target
signal as integers that corresponds to the position of the category
name in the ``target_names`` list. The category integer id of each
sample is stored in the ``target`` attribute::

  >>> training_set.target[:10]
  array([1, 0, 2, 2, 0, 1, 1, 3, 3, 2])

You can notices that the samples have been shuffled randomly (with
a fixed RNG seed): this is useful if you select only the first
samples to quickly train a model and get a first idea of the results
before re-training on the complete dataset later.


Extracting features from text files
-----------------------------------




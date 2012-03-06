The 20 newsgroups text dataset
==============================

The 20 newsgroups dataset comprises around 18000 newsgroups posts on
20 topics splitted in two subsets: one for training (or development)
and the other one for testing (or for performance evaluation). The split
between the train and test set is based upon a messages posted before
and after a specific date.

This module contains two loaders. The first one, 
``sklearn.datasets.fetch_20newsgroups``,
returns a list of the raw text files that can be fed to text feature
extractors such as :class:`sklearn.feature_extraction.text.Vectorizer`
with custom parameters so as to extract feature vectors.
The second one, ``sklearn.datasets.fetch_20newsgroups_vectorized``,
returns ready-to-use features, i.e., it is not necessary to use a feature
extractor.

Usage
-----

The ``sklearn.datasets.fetch_20newsgroups`` function is a data
fetching / caching functions that downloads the data archive from
the original `20 newsgroups website`_, extracts the archive contents
in the ``~/scikit_learn_data/20news_home`` folder and calls the
``sklearn.datasets.load_file`` on either the training or
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
list of the categories to load to the ``fetch_20newsgroups`` function::

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

In order to feed predictive or clustering models with the text data,
one first need to turn the text into vectors of numerical values suitable
for statistical analysis. This can be achieved with the utilities of the
``sklearn.feature_extraction.text`` as demonstrated in the following
example that extract `TF-IDF`_ vectors of unigram tokens::


  >>> from sklearn.feature_extraction.text import Vectorizer
  >>> documents = [open(f).read() for f in newsgroups_train.filenames]
  >>> vectorizer = Vectorizer()
  >>> vectors = vectorizer.fit_transform(documents)
  >>> vectors.shape
  (1073, 21108)

The extracted TF-IDF vectors are very sparse with an average of 118 non zero
components by sample in a more than 20000 dimensional space (less than 1% non
zero features)::

  >>> vectors.nnz / vectors.shape[0]
  118

``sklearn.datasets.fetch_20newsgroups_vectorized`` is a function which returns 
ready-to-use tfidf features instead of file names.

.. _`20 newsgroups website`: http://people.csail.mit.edu/jrennie/20Newsgroups/
.. _`TF-IDF`: http://en.wikipedia.org/wiki/Tf-idf


.. topic:: Examples

   * :ref:`example_grid_search_text_feature_extraction.py`

   * :ref:`example_document_classification_20newsgroups.py`



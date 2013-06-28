.. _20newsgroups:

The 20 newsgroups text dataset
==============================

The 20 newsgroups dataset comprises around 18000 newsgroups posts on
20 topics split in two subsets: one for training (or development)
and the other one for testing (or for performance evaluation). The split
between the train and test set is based upon a messages posted before
and after a specific date.

This module contains two loaders. The first one,
``sklearn.datasets.fetch_20newsgroups``,
returns a list of the raw texts that can be fed to text feature
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

Converting text to vectors
--------------------------

In order to feed predictive or clustering models with the text data,
one first need to turn the text into vectors of numerical values suitable
for statistical analysis. This can be achieved with the utilities of the
``sklearn.feature_extraction.text`` as demonstrated in the following
example that extract `TF-IDF`_ vectors of unigram tokens::

  >>> from sklearn.feature_extraction.text import TfidfVectorizer
  >>> newsgroups_train = fetch_20newsgroups(subset='train')
  >>> vectorizer = TfidfVectorizer()
  >>> vectors = vectorizer.fit_transform(newsgroups_train.data)
  >>> vectors.shape
  (11314, 129792)

The extracted TF-IDF vectors are very sparse, with an average of 110 non zero
components by sample in a more than 20000 dimensional space (less than 1% non
zero features)::

  >>> vectors.nnz / vectors.shape[0]
  110

``sklearn.datasets.fetch_20newsgroups_vectorized`` is a function which returns
ready-to-use tfidf features instead of file names.

.. _`20 newsgroups website`: http://people.csail.mit.edu/jrennie/20Newsgroups/
.. _`TF-IDF`: http://en.wikipedia.org/wiki/Tf-idf


Filtering text for more realistic training
------------------------------------------
It is easy for a classifier to overfit on particular things that appear in the
20 Newsgroups data, such as newsgroup headers. Many classifiers achieve very
high F-scores, but their results would not generalize to other documents that
aren't from this window of time.

For example, let's look at the results of a Bernoulli Naive Bayes classifier,
which is fast to train and achieves a decent F-score::

  >>> from sklearn.naive_bayes import BernoulliNB
  >>> from sklearn import metrics
  >>> newsgroups_test = fetch_20newsgroups(subset='test')
  >>> vectors_test = vectorizer.transform(newsgroups_test.data)
  >>> clf = BernoulliNB(alpha=.01)
  >>> clf.fit(vectors, newsgroups_train.target)
  >>> pred = clf.predict(vectors_test)
  >>> metrics.f1_score(pred, newsgroups_test.target)
  0.78117467868044399

(The example :ref:`example_document_classification_20newsgroups.py` shuffles
the training and test data, instead of segmenting by time, and in that case
Bernoulli Naive Bayes gets a much higher F-score of 0.88. Are you suspicious
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
  alt.atheism: god say think people don com nntp host posting article
  comp.graphics: like com article know thanks graphics university nntp host posting
  comp.os.ms-windows.misc: know thanks use com article nntp host posting university windows
  comp.sys.ibm.pc.hardware: just does article know thanks com university nntp host posting
  comp.sys.mac.hardware: thanks does apple know article mac university nntp host posting
  comp.windows.x: like article use reply thanks window com nntp posting host
  misc.forsale: com usa mail new distribution nntp host posting university sale
  rec.autos: distribution like just university nntp host posting car article com
  rec.motorcycles: don just like bike nntp host posting dod com article
  rec.sport.baseball: think just year com baseball host nntp posting university article
  rec.sport.hockey: nhl game hockey team nntp host article ca posting university
  sci.crypt: just nntp host encryption posting article chip key clipper com
  sci.electronics: does like know university use article com nntp host posting
  sci.med: like reply university nntp don host know posting com article
  sci.space: nasa like university just com nntp host posting space article
  soc.religion.christian: just don like rutgers know university article think people god
  talk.politics.guns: like university don gun people nntp host posting article com
  talk.politics.mideast: like just nntp host israeli university posting israel people article
  talk.politics.misc: like university nntp host just don posting people com article
  talk.religion.misc: think know christian posting god people just don article com

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
  ...                                      remove=('headers', 'footers', 'quotes'))
  >>> vectors_test = vectorizer.transform(newsgroups_test.data)
  >>> pred = clf.predict(vectors_test)
  >>> metrics.f1_score(pred, newsgroups_test.target)
  0.51830104911679742

This classifier lost over a third of its F-score, just because we removed
metadata that has little to do with topic classification. It recovers only a
bit if we also strip this metadata from the training data:

  >>> newsgroups_train = fetch_20newsgroups(subset='train',
                                            remove=('headers', 'footers', 'quotes'))
  >>> vectors = vectorizer.fit_transform(newsgroups_train.data)
  >>> clf = BernoulliNB(alpha=.01)
  >>> clf.fit(vectors, newsgroups_train.target)
  >>> vectors_test = vectorizer.transform(newsgroups_test.data)
  >>> pred = clf.predict(vectors_test)
  >>> metrics.f1_score(pred, newsgroups_test.target)
  0.56907392353755404

Some other classifiers cope better with this harder version of the task. Try
running :ref:`example_grid_search_text_feature_extraction.py` with and without
the ``--filter`` option to compare the results.

.. topic:: Recommendation

  When evaluating natural language classifiers on the 20 Newsgroups data, you
  should strip newsgroup-related metadata. In scikit-learn, you can do this by
  setting ``remove=('headers', 'footers', 'quotes')``. The F-score will be
  lower because it is more realistic.

.. topic:: Examples

   * :ref:`example_grid_search_text_feature_extraction.py`

   * :ref:`example_document_classification_20newsgroups.py`

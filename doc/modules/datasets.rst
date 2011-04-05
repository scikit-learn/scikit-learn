=========================
Dataset loading utilities
=========================

.. currentmodule:: scikits.learn.datasets

The ``scikits.learn.datasets`` package embeds some small toy datasets
as introduced in the "Getting Started" section.

To evaluate the impact of the scale of the dataset (``n_features`` and
``n_samples``) while controlling the statistical properties of the data
(typically the correlation and informativeness of the features), it is
also possible to generate synthetic data data

This package also features helpers to fetch larger datasets commonly
used by the machine learning community to benchmark algorithm on data
that comes from the 'real world'.


Dataset generators
==================

TODO



The Labeled Faces in the Wild face recognition dataset
======================================================

This dataset is a collection of JPEG pictures of famous people collected
over the internet, all details are available on the official website:

    http://vis-www.cs.umass.edu/lfw/

Each picture is centered on a single face. The typical task is called
Face Verification: given a pair of two pictures, a binary classifier
must predict whether the two images are from the same person.

An alternative task, Face Recognition or Face Identification is:
given the picture of the face of an unknown person, identify the name
of the person by refering to a gallery of previously seen pictures of
identified persons.

Both Face Verification and Face Recognition are tasks that are typically
performed on the output of a model trained to perform Face Detection. The
most popular model for Face Detection is called Viola-Johns and is
implemented in the OpenCV library. The LFW faces were extracted by this
face detector from various online websites.


Usage
-----

``scikit-learn`` provides two loaders that will automatically download,
cache, parse the metadata files, decode the jpeg and convert the
interesting slices into memmaped numpy arrays. This dataset size if more
than 200 MB. The first load typically takes more than a couple of minutes
to fully decode the relevant part of the JPEG files into numpy arrays. If
the dataset has  been loaded once, the following times the loading times
less than 200ms by using a memmaped version memoized on the disk in the
``~/scikit_learn_data/lfw_home/`` folder using ``joblib``.

The first loader is used for the Face Identification task: a multi-class
classification task (hence supervised learning)::

  >>> from scikits.learn.datasets import fetch_lfw_people
  >>> lfw_people = fetch_lfw_people(min_faces_per_person=70, resize=0.4)

  >>> for name in lfw_people.target_names:
  ...     print name
  ...
  Ariel Sharon
  Colin Powell
  Donald Rumsfeld
  George W Bush
  Gerhard Schroeder
  Hugo Chavez
  Tony Blair

The default slice is a rectangular shape around the face, removing
most of the background::

  >>> lfw_people.data.dtype
  dtype('float32')

  >>> lfw_people.data.shape
  (1288, 50, 37)

Each of the ``1140`` faces is assigned to a single person id in the ``target``
array::

  >>> lfw_people.target.shape
  (1288,)

  >>> list(lfw_people.target[:10])
  [5, 6, 3, 1, 0, 1, 3, 4, 3, 0]

The second loader is typically used for the face verification task: each sample
is a pair of two picture belonging or not to the same person::

  >>> from scikits.learn.datasets import fetch_lfw_pairs
  >>> lfw_pairs_train = fetch_lfw_pairs(subset='train')

  >>> list(lfw_pairs_train.target_names)
  ['Different persons', 'Same person']

  >>> lfw_pairs_train.data.shape
  (2200, 2, 62, 47)

  >>> lfw_pairs_train.target.shape
  (2200,)

Both for the ``fetch_lfw_people`` and ``fetch_lfw_pairs`` function it is
possible to get an additional dimension with the RGB color channels by
passing ``color=True``, in that case the shape will be
``(2200, 2, 62, 47, 3)``.

The ``fetch_lfw_pairs`` datasets is subdived in 3 subsets: the development
``train`` set, the development ``test`` set and an evaluation ``10_folds``
set meant to compute performance metrics using a 10-folds cross
validation scheme.

.. topic:: References:

 * `Labeled Faces in the Wild: A Database for Studying Face Recognition
   in Unconstrained Environments.
   <http://vis-www.cs.umass.edu/lfw/lfw.pdf>`_
   Gary B. Huang, Manu Ramesh, Tamara Berg, and Erik Learned-Miller.
   University of Massachusetts, Amherst, Technical Report 07-49, October, 2007.


Examples
--------

:ref:`example_applications_plot_face_recognition.py`



The 20 newsgroups text dataset
==============================

The 20 newsgroups dataset comprises around 18000 newsgroups posts on
20 topics splitted in two subsets: one for training (or development)
and the other one for testing (or for performance evaluation). The split
between the train and test set is based upon a messages posted before
and after a specific date.


Usage
-----

The ``scikits.learn.datasets.fetch_20newsgroups`` function is a data
fetching / caching functions that downloads the data archive from
the original `20 newsgroups website`_, extracts the archive contents
in the ``~/scikit_learn_data/20news_home`` folder and calls the
``scikits.learn.datasets.load_filenames`` on either the training or
testing set folder::

  >>> from scikits.learn.datasets import fetch_20newsgroups
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
``scikits.learn.feature_extraction.text`` as demonstrated in the following
example that extract `TF-IDF`_ vectors of unigram tokens::


  >>> from scikits.learn.feature_extraction.text import Vectorizer
  >>> documents = [open(f).read() for f in newsgroups_train.filenames]
  >>> vectorizer = Vectorizer()
  >>> vectors = vectorizer.fit_transform(documents)
  >>> vectors.shape
  (1073, 21107)

The extracted TF-IDF vectors are very sparse with an average of 118 non zero
components by sample in a more than 20000 dimensional space (less than 1% non
zero features)::

  >>> vectors.nnz / vectors.shape[0]
  118

.. _`20 newsgroups website`: http://people.csail.mit.edu/jrennie/20Newsgroups/
.. _`TF-IDF`: http://en.wikipedia.org/wiki/Tf-idf


Examples
--------

:ref:`example_grid_search_text_feature_extraction.py`

:ref:`example_document_classification_20newsgroups.py`


